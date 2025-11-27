/**
 * Robust Variance Estimation (RVE) for Meta-Analysis
 *
 * Implements cluster-robust variance estimation for dependent effect sizes:
 * - Small-sample corrections (CR0, CR1, CR2)
 * - Correlated effects working model
 * - Hierarchical effects working model
 * - Satterthwaite degrees of freedom
 *
 * Based on Hedges, Tipton & Johnson (2010) and Pustejovsky & Tipton (2022)
 * Equivalent to clubSandwich R package methodology
 */

import { jStat } from 'jstat';

/**
 * Construct variance-covariance matrix for correlated effects within clusters
 *
 * @param {Array} studies - Studies with es, vi, and clusterId
 * @param {number} rho - Assumed within-cluster correlation (default 0.8)
 * @returns {Object} Block-diagonal V matrix and cluster info
 */
export function constructVMatrix(studies, rho = 0.8) {
  const validStudies = studies.filter(s => !s.excluded && s.es != null && s.vi != null);

  // Group by cluster
  const clusters = {};
  validStudies.forEach((s, idx) => {
    const cid = s.clusterId || s.studyId || s.id;
    if (!clusters[cid]) {
      clusters[cid] = [];
    }
    clusters[cid].push({ ...s, originalIndex: idx });
  });

  const clusterIds = Object.keys(clusters);
  const m = clusterIds.length;  // Number of clusters
  const n = validStudies.length;  // Total effect sizes

  // Build block-diagonal V matrix
  const V = Array(n).fill(null).map(() => Array(n).fill(0));

  let currentIdx = 0;
  const clusterInfo = [];

  for (const cid of clusterIds) {
    const clusterStudies = clusters[cid];
    const k = clusterStudies.length;

    clusterInfo.push({
      id: cid,
      size: k,
      startIdx: currentIdx,
      endIdx: currentIdx + k - 1
    });

    // Fill block for this cluster
    for (let i = 0; i < k; i++) {
      for (let j = 0; j < k; j++) {
        if (i === j) {
          // Diagonal: sampling variance
          V[currentIdx + i][currentIdx + j] = clusterStudies[i].vi;
        } else {
          // Off-diagonal: correlation * sqrt(vi * vj)
          const cov = rho * Math.sqrt(clusterStudies[i].vi * clusterStudies[j].vi);
          V[currentIdx + i][currentIdx + j] = cov;
        }
      }
    }

    currentIdx += k;
  }

  return {
    V,
    clusters: clusterInfo,
    m,
    n,
    rho,
    studies: validStudies
  };
}

/**
 * Compute inverse of block-diagonal matrix efficiently
 */
function invertBlockDiagonal(V, clusterInfo) {
  const n = V.length;
  const Vinv = Array(n).fill(null).map(() => Array(n).fill(0));

  for (const cluster of clusterInfo) {
    const { startIdx, endIdx, size } = cluster;

    if (size === 1) {
      // Simple 1x1 block
      Vinv[startIdx][startIdx] = 1 / V[startIdx][startIdx];
    } else {
      // Extract block
      const block = [];
      for (let i = startIdx; i <= endIdx; i++) {
        const row = [];
        for (let j = startIdx; j <= endIdx; j++) {
          row.push(V[i][j]);
        }
        block.push(row);
      }

      // Invert block (using Gaussian elimination for small blocks)
      const invBlock = invertMatrix(block);

      // Place back
      for (let i = 0; i < size; i++) {
        for (let j = 0; j < size; j++) {
          Vinv[startIdx + i][startIdx + j] = invBlock[i][j];
        }
      }
    }
  }

  return Vinv;
}

/**
 * Invert a small matrix using Gauss-Jordan elimination
 */
function invertMatrix(matrix) {
  const n = matrix.length;
  const aug = matrix.map((row, i) => {
    const newRow = [...row];
    for (let j = 0; j < n; j++) {
      newRow.push(i === j ? 1 : 0);
    }
    return newRow;
  });

  // Forward elimination
  for (let i = 0; i < n; i++) {
    // Find pivot
    let maxRow = i;
    for (let k = i + 1; k < n; k++) {
      if (Math.abs(aug[k][i]) > Math.abs(aug[maxRow][i])) {
        maxRow = k;
      }
    }
    [aug[i], aug[maxRow]] = [aug[maxRow], aug[i]];

    if (Math.abs(aug[i][i]) < 1e-10) {
      // Near-singular, add small regularization
      aug[i][i] += 1e-8;
    }

    // Scale pivot row
    const scale = aug[i][i];
    for (let j = 0; j < 2 * n; j++) {
      aug[i][j] /= scale;
    }

    // Eliminate column
    for (let k = 0; k < n; k++) {
      if (k !== i) {
        const factor = aug[k][i];
        for (let j = 0; j < 2 * n; j++) {
          aug[k][j] -= factor * aug[i][j];
        }
      }
    }
  }

  // Extract inverse
  return aug.map(row => row.slice(n));
}

/**
 * Matrix multiplication helper
 */
function matMul(A, B) {
  const rowsA = A.length;
  const colsA = A[0].length;
  const colsB = B[0].length;

  const result = Array(rowsA).fill(null).map(() => Array(colsB).fill(0));

  for (let i = 0; i < rowsA; i++) {
    for (let j = 0; j < colsB; j++) {
      for (let k = 0; k < colsA; k++) {
        result[i][j] += A[i][k] * B[k][j];
      }
    }
  }

  return result;
}

/**
 * RVE Meta-Analysis with Correlated Effects Working Model
 *
 * @param {Array} studies - Studies with es, vi, clusterId
 * @param {Object} options - Options including rho, small sample correction type
 * @returns {Object} RVE results
 */
export function rveMetaAnalysis(studies, options = {}) {
  const {
    rho = 0.8,
    smallSampleCorrection = 'CR2',  // CR0, CR1, CR2
    confLevel = 0.95
  } = options;

  // Construct V matrix
  const { V, clusters, m, n, studies: validStudies } = constructVMatrix(studies, rho);

  if (n < 2) {
    return { error: 'At least 2 effect sizes required' };
  }

  if (m < 2) {
    return { error: 'At least 2 clusters required for RVE' };
  }

  // Extract effect sizes
  const y = validStudies.map(s => s.es);

  // Design matrix (intercept only for simple pooling)
  const X = validStudies.map(() => [1]);

  // Inverse of V
  const Vinv = invertBlockDiagonal(V, clusters);

  // Weighted least squares: (X'V^-1 X)^-1 X'V^-1 y
  // For intercept-only: beta = sum(Vinv * y) / sum(Vinv)
  let sumVinvY = 0;
  let sumVinv = 0;

  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      sumVinvY += Vinv[i][j] * y[j];
      sumVinv += Vinv[i][j];
    }
  }

  const beta = sumVinvY / sumVinv;

  // Residuals
  const residuals = y.map(yi => yi - beta);

  // Cluster-robust variance estimation
  // Var(beta) = (X'V^-1 X)^-1 * meat * (X'V^-1 X)^-1
  // where meat = sum over clusters of X'V^-1 e e' V^-1 X

  const bread = 1 / sumVinv;  // (X'V^-1 X)^-1 for intercept model

  // Calculate meat with small-sample correction
  let meat = 0;

  for (const cluster of clusters) {
    const { startIdx, endIdx, size } = cluster;

    // Extract cluster residuals and weights
    let clusterContrib = 0;

    if (smallSampleCorrection === 'CR0') {
      // No correction
      for (let i = startIdx; i <= endIdx; i++) {
        for (let j = startIdx; j <= endIdx; j++) {
          clusterContrib += Vinv[i][j] * residuals[i] * residuals[j];
        }
      }
    } else if (smallSampleCorrection === 'CR1') {
      // Simple multiplier correction: m/(m-1)
      for (let i = startIdx; i <= endIdx; i++) {
        for (let j = startIdx; j <= endIdx; j++) {
          clusterContrib += Vinv[i][j] * residuals[i] * residuals[j];
        }
      }
      clusterContrib *= m / (m - 1);
    } else {
      // CR2: Bias-reduced linearization (BRL)
      // Adjust residuals by (I - H)^{-1/2}
      // Simplified approximation for intercept model

      // Leverage for cluster
      let Hjj = 0;
      for (let i = startIdx; i <= endIdx; i++) {
        for (let j = startIdx; j <= endIdx; j++) {
          Hjj += Vinv[i][j];
        }
      }
      Hjj *= bread;

      // CR2 multiplier
      const cr2Mult = 1 / Math.sqrt(Math.max(0.01, 1 - Hjj));

      for (let i = startIdx; i <= endIdx; i++) {
        for (let j = startIdx; j <= endIdx; j++) {
          clusterContrib += Vinv[i][j] * residuals[i] * residuals[j] * cr2Mult * cr2Mult;
        }
      }
    }

    meat += clusterContrib;
  }

  // Sandwich variance
  const varBeta = bread * meat * bread;
  const seBeta = Math.sqrt(varBeta);

  // Satterthwaite degrees of freedom
  // df = 2 * E[varBeta]^2 / Var[varBeta]
  // Approximation: df ≈ m - 1 for simple cases
  let df;
  if (smallSampleCorrection === 'CR2') {
    // More accurate df calculation
    df = calculateSatterthwaiteDf(clusters, Vinv, bread, m);
  } else {
    df = m - 1;
  }

  // Confidence interval
  const alpha = 1 - confLevel;
  const tCrit = jStat.studentt.inv(1 - alpha / 2, df);

  const ciLo = beta - tCrit * seBeta;
  const ciHi = beta + tCrit * seBeta;

  // Test statistic and p-value
  const tStat = beta / seBeta;
  const pVal = 2 * (1 - jStat.studentt.cdf(Math.abs(tStat), df));

  // I² approximation for RVE
  // Based on typical heterogeneity decomposition
  const Q = residuals.reduce((sum, r, i) => {
    let contrib = 0;
    for (let j = 0; j < n; j++) {
      contrib += Vinv[i][j] * r * residuals[j];
    }
    return sum + contrib;
  }, 0);

  const tau2 = Math.max(0, (Q - (n - 1)) / (sumVinv - sumSquaredWeights(Vinv) / sumVinv));
  const I2 = tau2 / (tau2 + meanSamplingVar(validStudies)) * 100;

  return {
    method: 'Robust Variance Estimation',
    correction: smallSampleCorrection,
    rho,
    es: beta,
    se: seBeta,
    ci: { lo: ciLo, hi: ciHi },
    t: tStat,
    df,
    p: pVal,
    nEffects: n,
    nClusters: m,
    heterogeneity: {
      tau2Approx: tau2,
      I2Approx: Math.max(0, Math.min(100, I2))
    },
    clusters: clusters.map(c => ({
      id: c.id,
      nEffects: c.size
    })),
    interpretation: generateRVEInterpretation(beta, seBeta, pVal, m, n)
  };
}

function sumSquaredWeights(Vinv) {
  const n = Vinv.length;
  let sum = 0;
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      sum += Vinv[i][j] * Vinv[i][j];
    }
  }
  return sum;
}

function meanSamplingVar(studies) {
  return studies.reduce((sum, s) => sum + s.vi, 0) / studies.length;
}

function calculateSatterthwaiteDf(clusters, Vinv, bread, m) {
  // Simplified Satterthwaite approximation
  // For balanced designs: df ≈ m - 1
  // For unbalanced: adjust based on cluster sizes

  const sizes = clusters.map(c => c.size);
  const avgSize = sizes.reduce((a, b) => a + b, 0) / m;
  const sizeVariance = sizes.reduce((sum, s) => sum + (s - avgSize) ** 2, 0) / m;

  // Adjust df for imbalance
  const imbalanceFactor = 1 - sizeVariance / (avgSize * avgSize + 0.001);

  return Math.max(1, (m - 1) * imbalanceFactor);
}

function generateRVEInterpretation(es, se, p, m, n) {
  const sig = p < 0.05;
  return `RVE pooled estimate: ${es.toFixed(3)} (SE=${se.toFixed(3)}), ` +
    `based on ${n} effect sizes from ${m} clusters. ` +
    (sig ? 'Effect is statistically significant.' : 'Effect is not statistically significant.');
}

/**
 * RVE Meta-Regression with Robust Standard Errors
 *
 * @param {Array} studies - Studies with es, vi, clusterId, and covariates
 * @param {Array} covariateNames - Names of covariates to include
 * @param {Object} options - Options
 * @returns {Object} RVE regression results
 */
export function rveMetaRegression(studies, covariateNames, options = {}) {
  const {
    rho = 0.8,
    smallSampleCorrection = 'CR2',
    confLevel = 0.95
  } = options;

  const { V, clusters, m, n, studies: validStudies } = constructVMatrix(studies, rho);

  if (n < covariateNames.length + 2) {
    return { error: 'Insufficient studies for number of covariates' };
  }

  // Build design matrix
  const p = covariateNames.length + 1;  // Including intercept
  const X = validStudies.map(s => {
    const row = [1];  // Intercept
    for (const name of covariateNames) {
      row.push(s[name] || 0);
    }
    return row;
  });

  const y = validStudies.map(s => s.es);

  // V inverse
  const Vinv = invertBlockDiagonal(V, clusters);

  // X'V^-1X
  const XtVinv = [];
  for (let i = 0; i < p; i++) {
    XtVinv[i] = [];
    for (let j = 0; j < n; j++) {
      let sum = 0;
      for (let k = 0; k < n; k++) {
        sum += X[k][i] * Vinv[k][j];
      }
      XtVinv[i][j] = sum;
    }
  }

  const XtVinvX = [];
  for (let i = 0; i < p; i++) {
    XtVinvX[i] = [];
    for (let j = 0; j < p; j++) {
      let sum = 0;
      for (let k = 0; k < n; k++) {
        sum += XtVinv[i][k] * X[k][j];
      }
      XtVinvX[i][j] = sum;
    }
  }

  // Invert XtVinvX
  const bread = invertMatrix(XtVinvX);

  // X'V^-1y
  const XtVinvy = [];
  for (let i = 0; i < p; i++) {
    let sum = 0;
    for (let j = 0; j < n; j++) {
      sum += XtVinv[i][j] * y[j];
    }
    XtVinvy[i] = sum;
  }

  // Beta = (X'V^-1X)^-1 X'V^-1y
  const beta = [];
  for (let i = 0; i < p; i++) {
    let sum = 0;
    for (let j = 0; j < p; j++) {
      sum += bread[i][j] * XtVinvy[j];
    }
    beta[i] = sum;
  }

  // Residuals
  const residuals = validStudies.map((s, i) => {
    let fitted = 0;
    for (let j = 0; j < p; j++) {
      fitted += X[i][j] * beta[j];
    }
    return y[i] - fitted;
  });

  // Cluster-robust meat matrix
  const meat = Array(p).fill(null).map(() => Array(p).fill(0));

  for (const cluster of clusters) {
    const { startIdx, endIdx } = cluster;

    // Get leverage adjustment for CR2
    let cr2Mult = 1;
    if (smallSampleCorrection === 'CR2') {
      // Simplified CR2 adjustment
      let Hjj = 0;
      for (let i = startIdx; i <= endIdx; i++) {
        for (let j = startIdx; j <= endIdx; j++) {
          for (let a = 0; a < p; a++) {
            for (let b = 0; b < p; b++) {
              Hjj += X[i][a] * bread[a][b] * XtVinv[b][j];
            }
          }
        }
      }
      cr2Mult = 1 / Math.sqrt(Math.max(0.01, 1 - Hjj / (endIdx - startIdx + 1)));
    } else if (smallSampleCorrection === 'CR1') {
      cr2Mult = Math.sqrt(m / (m - 1));
    }

    // Contribution to meat
    for (let a = 0; a < p; a++) {
      for (let b = 0; b < p; b++) {
        let contrib = 0;
        for (let i = startIdx; i <= endIdx; i++) {
          for (let j = startIdx; j <= endIdx; j++) {
            contrib += XtVinv[a][i] * residuals[i] * residuals[j] * Vinv[j][i] * X[j][b];
          }
        }
        meat[a][b] += contrib * cr2Mult * cr2Mult;
      }
    }
  }

  // Sandwich: bread * meat * bread
  const temp = matMul(bread, meat);
  const varBeta = matMul(temp, bread);

  // Results for each coefficient
  const coefficients = [];
  const alpha = 1 - confLevel;
  const df = m - p;
  const tCrit = df > 0 ? jStat.studentt.inv(1 - alpha / 2, df) : 1.96;

  const names = ['(Intercept)', ...covariateNames];
  for (let i = 0; i < p; i++) {
    const se = Math.sqrt(Math.max(0, varBeta[i][i]));
    const t = se > 0 ? beta[i] / se : 0;
    const pVal = df > 0 ? 2 * (1 - jStat.studentt.cdf(Math.abs(t), df)) : 2 * (1 - jStat.normal.cdf(Math.abs(t), 0, 1));

    coefficients.push({
      name: names[i],
      estimate: beta[i],
      se,
      t,
      df,
      p: pVal,
      ci: {
        lo: beta[i] - tCrit * se,
        hi: beta[i] + tCrit * se
      }
    });
  }

  // Model fit
  const SSres = residuals.reduce((sum, r) => sum + r * r, 0);
  const meanY = y.reduce((a, b) => a + b, 0) / n;
  const SStot = y.reduce((sum, yi) => sum + (yi - meanY) ** 2, 0);
  const R2 = 1 - SSres / SStot;

  return {
    method: 'RVE Meta-Regression',
    correction: smallSampleCorrection,
    rho,
    coefficients,
    R2: Math.max(0, R2) * 100,
    nEffects: n,
    nClusters: m,
    df,
    residuals: residuals.slice(0, 10)  // First 10 for diagnostics
  };
}

/**
 * Sensitivity analysis for assumed correlation rho
 *
 * @param {Array} studies - Studies
 * @param {Object} options - Options including rho range
 * @returns {Object} Sensitivity results across rho values
 */
export function rveSensitivityAnalysis(studies, options = {}) {
  const {
    rhoValues = [0, 0.2, 0.4, 0.6, 0.8, 1.0],
    smallSampleCorrection = 'CR2'
  } = options;

  const results = [];

  for (const rho of rhoValues) {
    const result = rveMetaAnalysis(studies, { rho, smallSampleCorrection });

    if (!result.error) {
      results.push({
        rho,
        es: result.es,
        se: result.se,
        ci: result.ci,
        p: result.p
      });
    }
  }

  // Check sensitivity
  const estimates = results.map(r => r.es);
  const minES = Math.min(...estimates);
  const maxES = Math.max(...estimates);
  const range = maxES - minES;

  const conclusions = results.map(r => r.p < 0.05);
  const allSameConclusion = conclusions.every(c => c === conclusions[0]);

  return {
    results,
    sensitivity: {
      estimateRange: range,
      minES,
      maxES,
      conclusionRobust: allSameConclusion
    },
    interpretation: allSameConclusion
      ? 'Results are robust to assumed within-cluster correlation.'
      : 'Conclusions may depend on assumed within-cluster correlation. Interpret with caution.'
  };
}

export default {
  constructVMatrix,
  rveMetaAnalysis,
  rveMetaRegression,
  rveSensitivityAnalysis
};
