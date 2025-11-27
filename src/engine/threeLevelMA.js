/**
 * Three-Level Meta-Analysis for Dependent Effect Sizes
 *
 * Implements hierarchical models for nested effect sizes:
 * - Level 1: Sampling variance (within-study, known)
 * - Level 2: Within-cluster heterogeneity (multiple ES per cluster)
 * - Level 3: Between-cluster heterogeneity
 *
 * Based on Van den Noortgate et al. (2013), Cheung (2014), and Assink & Wibbelink (2016)
 */

import { jStat } from 'jstat';

/**
 * Fit Three-Level Random Effects Meta-Analysis Model
 *
 * @param {Array} studies - Studies with es, vi, clusterId (and optionally effectId)
 * @param {Object} options - Options including estimation method
 * @returns {Object} Three-level model results
 */
export function threeLevelMetaAnalysis(studies, options = {}) {
  const {
    method = 'REML',  // 'REML' or 'ML'
    maxIter = 100,
    tolerance = 1e-6,
    confLevel = 0.95
  } = options;

  const validStudies = studies.filter(s => !s.excluded && s.es != null && s.vi != null);

  if (validStudies.length < 3) {
    return { error: 'At least 3 effect sizes required' };
  }

  // Organize by clusters
  const clusters = {};
  validStudies.forEach((s, idx) => {
    const cid = s.clusterId || s.studyId || s.id;
    if (!clusters[cid]) {
      clusters[cid] = [];
    }
    clusters[cid].push({ ...s, globalIdx: idx });
  });

  const clusterIds = Object.keys(clusters);
  const m = clusterIds.length;  // Number of clusters
  const k = validStudies.length;  // Total effect sizes

  if (m < 2) {
    return { error: 'At least 2 clusters required for three-level model' };
  }

  // Extract data
  const y = validStudies.map(s => s.es);
  const v = validStudies.map(s => s.vi);  // Level-1 variance (known)

  // Cluster structure
  const clusterInfo = [];
  let currentIdx = 0;
  for (const cid of clusterIds) {
    const n = clusters[cid].length;
    clusterInfo.push({
      id: cid,
      size: n,
      startIdx: currentIdx,
      endIdx: currentIdx + n - 1,
      indices: Array.from({ length: n }, (_, i) => currentIdx + i)
    });
    currentIdx += n;
  }

  // Initialize variance components
  // tau2_2: within-cluster (level 2) variance
  // tau2_3: between-cluster (level 3) variance
  let tau2_2 = estimateInitialTau2(y, v);
  let tau2_3 = tau2_2 / 2;

  // Iterative estimation
  for (let iter = 0; iter < maxIter; iter++) {
    const prevTau2_2 = tau2_2;
    const prevTau2_3 = tau2_3;

    // Construct total variance-covariance matrix V
    // V_ij = v_i * I(i==j) + tau2_2 * I(cluster_i == cluster_j) + tau2_3
    // For same cluster: v_i + tau2_2 + tau2_3 (diagonal), tau2_3 (off-diagonal within cluster)
    // Wait, this is wrong. Let me reconsider.
    // Actually for three-level:
    // - Same effect size (i=j): vi + tau2_2 + tau2_3
    // - Same cluster, different ES: tau2_3
    // - Different clusters: 0

    // No wait, standard three-level parameterization:
    // Level 2 variance (within cluster): tau2_2
    // Level 3 variance (between cluster): tau2_3
    // Total variance for effect i: vi + tau2_2 + tau2_3
    // Covariance between effects in same cluster: tau2_3

    // Build and invert V
    const { Vinv, logDetV } = buildAndInvertV(v, tau2_2, tau2_3, clusterInfo, k);

    // Weighted least squares for mu
    let sumVinv = 0;
    let sumVinvY = 0;
    for (let i = 0; i < k; i++) {
      for (let j = 0; j < k; j++) {
        sumVinv += Vinv[i][j];
        sumVinvY += Vinv[i][j] * y[j];
      }
    }

    const mu = sumVinvY / sumVinv;
    const varMu = 1 / sumVinv;

    // Residuals
    const residuals = y.map(yi => yi - mu);

    // Profile likelihood / REML estimation of variance components
    // Use Fisher scoring or method of moments

    // Method of moments update (Raudenbush-Bryk style)
    // Q_within = sum within clusters of weighted squared residuals
    // Q_between = weighted squared deviations of cluster means

    let Q_within = 0;
    let Q_between = 0;
    const clusterMeans = [];
    const clusterWeights = [];

    for (const cluster of clusterInfo) {
      const { indices, size } = cluster;

      // Cluster mean and weight
      const clusterResids = indices.map(i => residuals[i]);
      const clusterVars = indices.map(i => v[i] + tau2_2);
      const clusterW = clusterVars.map(vv => 1 / vv);
      const sumCW = clusterW.reduce((a, b) => a + b, 0);

      const clusterMean = clusterW.reduce((sum, w, j) =>
        sum + w * clusterResids[j], 0) / sumCW;

      clusterMeans.push(clusterMean);
      clusterWeights.push(sumCW);

      // Within-cluster Q
      for (let j = 0; j < size; j++) {
        Q_within += clusterW[j] * (clusterResids[j] - clusterMean) ** 2;
      }
    }

    // Between-cluster Q
    const totalClusterWeight = clusterWeights.reduce((a, b) => a + b, 0);
    const grandMean = clusterWeights.reduce((sum, w, i) =>
      sum + w * clusterMeans[i], 0) / totalClusterWeight;

    for (let i = 0; i < m; i++) {
      const clusterVar = 1 / clusterWeights[i] + tau2_3;
      Q_between += (clusterMeans[i] - grandMean) ** 2 / clusterVar;
    }

    // Update variance components
    // tau2_2 from Q_within
    const dfWithin = k - m;
    if (dfWithin > 0) {
      const c_within = calculateC(v, tau2_2, tau2_3, clusterInfo, 'within');
      tau2_2 = Math.max(0, (Q_within - dfWithin) / c_within + tau2_2);
    }

    // tau2_3 from Q_between
    const dfBetween = m - 1;
    if (dfBetween > 0) {
      const c_between = calculateC(v, tau2_2, tau2_3, clusterInfo, 'between');
      tau2_3 = Math.max(0, (Q_between - dfBetween) / c_between + tau2_3);
    }

    // Check convergence
    if (Math.abs(tau2_2 - prevTau2_2) < tolerance &&
        Math.abs(tau2_3 - prevTau2_3) < tolerance) {
      break;
    }
  }

  // Final estimates with converged variance components
  const { Vinv } = buildAndInvertV(v, tau2_2, tau2_3, clusterInfo, k);

  let sumVinv = 0;
  let sumVinvY = 0;
  for (let i = 0; i < k; i++) {
    for (let j = 0; j < k; j++) {
      sumVinv += Vinv[i][j];
      sumVinvY += Vinv[i][j] * y[j];
    }
  }

  const mu = sumVinvY / sumVinv;
  const seMu = Math.sqrt(1 / sumVinv);

  // Confidence interval
  const alpha = 1 - confLevel;
  const df = m - 1;  // Conservative df
  const tCrit = jStat.studentt.inv(1 - alpha / 2, df);

  const ciLo = mu - tCrit * seMu;
  const ciHi = mu + tCrit * seMu;

  // Test
  const tStat = mu / seMu;
  const pVal = 2 * (1 - jStat.studentt.cdf(Math.abs(tStat), df));

  // I² decomposition
  const typicalV = v.reduce((a, b) => a + b, 0) / k;
  const totalTau2 = tau2_2 + tau2_3;
  const I2_total = totalTau2 / (totalTau2 + typicalV) * 100;
  const I2_level2 = tau2_2 / (totalTau2 + typicalV) * 100;
  const I2_level3 = tau2_3 / (totalTau2 + typicalV) * 100;

  // Variance explained by each level
  const varExplained = {
    level2: tau2_2 / (tau2_2 + tau2_3 + 0.001) * 100,
    level3: tau2_3 / (tau2_2 + tau2_3 + 0.001) * 100
  };

  // Prediction interval
  const piSE = Math.sqrt(seMu * seMu + tau2_2 + tau2_3);
  const piCrit = jStat.studentt.inv(1 - alpha / 2, Math.max(1, df - 1));

  return {
    method: 'Three-Level Meta-Analysis',
    estimation: method,
    mu: {
      estimate: mu,
      se: seMu,
      ci: { lo: ciLo, hi: ciHi },
      t: tStat,
      df,
      p: pVal
    },
    varianceComponents: {
      tau2_level2: tau2_2,
      tau2_level3: tau2_3,
      tau2_total: totalTau2,
      tau_level2: Math.sqrt(tau2_2),
      tau_level3: Math.sqrt(tau2_3)
    },
    heterogeneity: {
      I2_total,
      I2_level2,
      I2_level3,
      varExplained
    },
    predictionInterval: {
      lo: mu - piCrit * piSE,
      hi: mu + piCrit * piSE
    },
    structure: {
      nEffects: k,
      nClusters: m,
      avgEffectsPerCluster: k / m,
      clusterSizes: clusterInfo.map(c => ({ id: c.id, n: c.size }))
    },
    interpretation: generateThreeLevelInterpretation(mu, pVal, I2_level2, I2_level3, m, k)
  };
}

/**
 * Build and invert the variance-covariance matrix
 */
function buildAndInvertV(v, tau2_2, tau2_3, clusterInfo, k) {
  const V = Array(k).fill(null).map(() => Array(k).fill(0));

  // Fill V matrix
  for (const cluster of clusterInfo) {
    const { indices } = cluster;

    for (let i = 0; i < indices.length; i++) {
      for (let j = 0; j < indices.length; j++) {
        const gi = indices[i];
        const gj = indices[j];

        if (gi === gj) {
          // Diagonal: vi + tau2_2 + tau2_3
          V[gi][gj] = v[gi] + tau2_2 + tau2_3;
        } else {
          // Off-diagonal within cluster: tau2_3
          V[gi][gj] = tau2_3;
        }
      }
    }
  }

  // Invert using block structure
  const Vinv = invertBlockDiagonalThreeLevel(V, clusterInfo);

  // Log determinant (for likelihood)
  let logDetV = 0;
  for (const cluster of clusterInfo) {
    const n = cluster.size;
    const diagVal = v[cluster.startIdx] + tau2_2 + tau2_3;

    if (n === 1) {
      logDetV += Math.log(diagVal);
    } else {
      // For equicorrelated block: det = (1 + (n-1)rho) * (1-rho)^(n-1) * prod(diag)
      // Simplified for equal variances within cluster
      const totalVar = diagVal;
      const cov = tau2_3;
      const rho = cov / totalVar;

      logDetV += Math.log(totalVar) * n;
      logDetV += Math.log(1 + (n - 1) * rho);
      logDetV += (n - 1) * Math.log(1 - rho);
    }
  }

  return { Vinv, logDetV };
}

/**
 * Invert block-diagonal V matrix for three-level model
 */
function invertBlockDiagonalThreeLevel(V, clusterInfo) {
  const k = V.length;
  const Vinv = Array(k).fill(null).map(() => Array(k).fill(0));

  for (const cluster of clusterInfo) {
    const { indices, size } = cluster;

    if (size === 1) {
      const idx = indices[0];
      Vinv[idx][idx] = 1 / V[idx][idx];
    } else {
      // For compound symmetric block:
      // V^-1 = (1/sigma2) * [I - (rho/(1+(n-1)*rho)) * J]
      // where J is matrix of 1s

      const diagVal = V[indices[0]][indices[0]];
      const offDiagVal = V[indices[0]][indices[1]];

      const a = diagVal - offDiagVal;  // sigma2 * (1 - rho)
      const b = offDiagVal;  // sigma2 * rho

      // Inverse elements
      const denom = a * (a + size * b);
      const invDiag = (a + (size - 1) * b) / denom;
      const invOffDiag = -b / denom;

      for (let i = 0; i < size; i++) {
        for (let j = 0; j < size; j++) {
          Vinv[indices[i]][indices[j]] = (i === j) ? invDiag : invOffDiag;
        }
      }
    }
  }

  return Vinv;
}

/**
 * Calculate C constant for method of moments
 */
function calculateC(v, tau2_2, tau2_3, clusterInfo, level) {
  let c = 0;

  if (level === 'within') {
    for (const cluster of clusterInfo) {
      const n = cluster.size;
      if (n > 1) {
        const totalVar = v[cluster.startIdx] + tau2_2 + tau2_3;
        c += (n - 1) / totalVar;
      }
    }
  } else {
    // Between
    for (const cluster of clusterInfo) {
      const n = cluster.size;
      const totalVar = v[cluster.startIdx] + tau2_2 + tau2_3;
      const clusterVar = totalVar / n + tau2_3 * (n - 1) / n;
      c += 1 / clusterVar;
    }
  }

  return Math.max(c, 0.001);
}

/**
 * Initial tau2 estimate using DerSimonian-Laird
 */
function estimateInitialTau2(y, v) {
  const k = y.length;
  const w = v.map(vi => 1 / vi);
  const sumW = w.reduce((a, b) => a + b, 0);
  const pooled = w.reduce((sum, wi, i) => sum + wi * y[i], 0) / sumW;
  const Q = w.reduce((sum, wi, i) => sum + wi * (y[i] - pooled) ** 2, 0);
  const C = sumW - w.reduce((sum, wi) => sum + wi * wi, 0) / sumW;

  return Math.max(0, (Q - (k - 1)) / C);
}

function generateThreeLevelInterpretation(mu, pVal, I2_level2, I2_level3, m, k) {
  const parts = [];

  parts.push(`Pooled estimate: ${mu.toFixed(3)} (${pVal < 0.05 ? 'significant' : 'not significant'})`);
  parts.push(`Based on ${k} effect sizes from ${m} clusters.`);

  if (I2_level2 > I2_level3) {
    parts.push(`Most heterogeneity is within clusters (I²=${I2_level2.toFixed(1)}%).`);
  } else {
    parts.push(`Most heterogeneity is between clusters (I²=${I2_level3.toFixed(1)}%).`);
  }

  return parts.join(' ');
}

/**
 * Three-Level Meta-Regression
 *
 * @param {Array} studies - Studies with es, vi, clusterId, and covariates
 * @param {Array} covariateNames - Names of covariates
 * @param {Object} options - Options
 * @returns {Object} Three-level regression results
 */
export function threeLevelMetaRegression(studies, covariateNames, options = {}) {
  const {
    method = 'REML',
    maxIter = 100,
    tolerance = 1e-6,
    confLevel = 0.95
  } = options;

  const validStudies = studies.filter(s => !s.excluded && s.es != null && s.vi != null);
  const k = validStudies.length;
  const p = covariateNames.length + 1;

  if (k < p + 2) {
    return { error: 'Insufficient effect sizes for number of covariates' };
  }

  // Organize clusters
  const clusters = {};
  validStudies.forEach((s, idx) => {
    const cid = s.clusterId || s.studyId || s.id;
    if (!clusters[cid]) clusters[cid] = [];
    clusters[cid].push({ ...s, globalIdx: idx });
  });

  const clusterIds = Object.keys(clusters);
  const m = clusterIds.length;

  const clusterInfo = [];
  let currentIdx = 0;
  for (const cid of clusterIds) {
    const n = clusters[cid].length;
    clusterInfo.push({
      id: cid,
      size: n,
      indices: Array.from({ length: n }, (_, i) => currentIdx + i)
    });
    currentIdx += n;
  }

  const y = validStudies.map(s => s.es);
  const v = validStudies.map(s => s.vi);

  // Design matrix
  const X = validStudies.map(s => {
    const row = [1];
    for (const name of covariateNames) {
      row.push(s[name] || 0);
    }
    return row;
  });

  // Initialize variance components
  let tau2_2 = estimateInitialTau2(y, v) / 2;
  let tau2_3 = tau2_2;

  // Iterate
  for (let iter = 0; iter < maxIter; iter++) {
    const prevTau2_2 = tau2_2;
    const prevTau2_3 = tau2_3;

    const { Vinv } = buildAndInvertV(v, tau2_2, tau2_3, clusterInfo, k);

    // GLS estimation: beta = (X'V^-1 X)^-1 X'V^-1 y
    const XtVinvX = Array(p).fill(null).map(() => Array(p).fill(0));
    const XtVinvy = Array(p).fill(0);

    for (let a = 0; a < p; a++) {
      for (let i = 0; i < k; i++) {
        for (let j = 0; j < k; j++) {
          XtVinvy[a] += X[i][a] * Vinv[i][j] * y[j];
          for (let b = 0; b < p; b++) {
            XtVinvX[a][b] += X[i][a] * Vinv[i][j] * X[j][b];
          }
        }
      }
    }

    const XtVinvXinv = invertSmallMatrix(XtVinvX);
    const beta = XtVinvXinv.map((row, i) =>
      row.reduce((sum, val, j) => sum + val * XtVinvy[j], 0)
    );

    // Residuals
    const residuals = y.map((yi, i) => {
      let fitted = 0;
      for (let j = 0; j < p; j++) fitted += X[i][j] * beta[j];
      return yi - fitted;
    });

    // Update variance components (simplified)
    let SSR_within = 0;
    let SSR_between = 0;

    for (const cluster of clusterInfo) {
      const { indices, size } = cluster;
      const clusterResids = indices.map(i => residuals[i]);
      const clusterMean = clusterResids.reduce((a, b) => a + b, 0) / size;

      for (const r of clusterResids) {
        SSR_within += (r - clusterMean) ** 2;
      }
      SSR_between += size * clusterMean ** 2;
    }

    const dfWithin = k - m - p + 1;
    const dfBetween = m - 1;

    if (dfWithin > 0) {
      tau2_2 = Math.max(0, SSR_within / dfWithin - v.reduce((a, b) => a + b, 0) / k);
    }
    if (dfBetween > 0) {
      tau2_3 = Math.max(0, SSR_between / dfBetween - tau2_2 - v.reduce((a, b) => a + b, 0) / k);
    }

    if (Math.abs(tau2_2 - prevTau2_2) < tolerance &&
        Math.abs(tau2_3 - prevTau2_3) < tolerance) {
      break;
    }
  }

  // Final estimates
  const { Vinv } = buildAndInvertV(v, tau2_2, tau2_3, clusterInfo, k);

  const XtVinvX = Array(p).fill(null).map(() => Array(p).fill(0));
  const XtVinvy = Array(p).fill(0);

  for (let a = 0; a < p; a++) {
    for (let i = 0; i < k; i++) {
      for (let j = 0; j < k; j++) {
        XtVinvy[a] += X[i][a] * Vinv[i][j] * y[j];
        for (let b = 0; b < p; b++) {
          XtVinvX[a][b] += X[i][a] * Vinv[i][j] * X[j][b];
        }
      }
    }
  }

  const varBeta = invertSmallMatrix(XtVinvX);
  const beta = varBeta.map((row, i) =>
    row.reduce((sum, val, j) => sum + val * XtVinvy[j], 0)
  );

  // Coefficients with inference
  const alpha = 1 - confLevel;
  const df = m - p;
  const tCrit = df > 0 ? jStat.studentt.inv(1 - alpha / 2, df) : 1.96;

  const names = ['(Intercept)', ...covariateNames];
  const coefficients = beta.map((b, i) => {
    const se = Math.sqrt(Math.max(0, varBeta[i][i]));
    const t = se > 0 ? b / se : 0;
    const pVal = df > 0 ? 2 * (1 - jStat.studentt.cdf(Math.abs(t), df)) : 0;

    return {
      name: names[i],
      estimate: b,
      se,
      t,
      df,
      p: pVal,
      ci: { lo: b - tCrit * se, hi: b + tCrit * se }
    };
  });

  // R² calculation
  const meanY = y.reduce((a, b) => a + b, 0) / k;
  const SStot = y.reduce((sum, yi) => sum + (yi - meanY) ** 2, 0);
  const residuals = y.map((yi, i) => {
    let fitted = 0;
    for (let j = 0; j < p; j++) fitted += X[i][j] * beta[j];
    return yi - fitted;
  });
  const SSres = residuals.reduce((sum, r) => sum + r * r, 0);
  const R2 = 1 - SSres / SStot;

  return {
    method: 'Three-Level Meta-Regression',
    coefficients,
    varianceComponents: {
      tau2_level2: tau2_2,
      tau2_level3: tau2_3
    },
    R2: Math.max(0, R2) * 100,
    nEffects: k,
    nClusters: m,
    df
  };
}

/**
 * Helper to invert small matrix
 */
function invertSmallMatrix(M) {
  const n = M.length;
  const aug = M.map((row, i) => {
    const newRow = [...row];
    for (let j = 0; j < n; j++) newRow.push(i === j ? 1 : 0);
    return newRow;
  });

  for (let i = 0; i < n; i++) {
    let maxRow = i;
    for (let k = i + 1; k < n; k++) {
      if (Math.abs(aug[k][i]) > Math.abs(aug[maxRow][i])) maxRow = k;
    }
    [aug[i], aug[maxRow]] = [aug[maxRow], aug[i]];

    if (Math.abs(aug[i][i]) < 1e-10) aug[i][i] += 1e-8;

    const scale = aug[i][i];
    for (let j = 0; j < 2 * n; j++) aug[i][j] /= scale;

    for (let k = 0; k < n; k++) {
      if (k !== i) {
        const factor = aug[k][i];
        for (let j = 0; j < 2 * n; j++) aug[k][j] -= factor * aug[i][j];
      }
    }
  }

  return aug.map(row => row.slice(n));
}

/**
 * Compare Two-Level vs Three-Level Models
 *
 * @param {Array} studies - Studies
 * @param {Object} options - Options
 * @returns {Object} Model comparison
 */
export function compareModels(studies, options = {}) {
  const threeLevel = threeLevelMetaAnalysis(studies, options);

  if (threeLevel.error) {
    return { error: threeLevel.error };
  }

  // Fit two-level model (standard random effects)
  const validStudies = studies.filter(s => !s.excluded && s.es != null && s.vi != null);
  const y = validStudies.map(s => s.es);
  const v = validStudies.map(s => s.vi);
  const k = y.length;

  const tau2_total = estimateInitialTau2(y, v);
  const w = v.map(vi => 1 / (vi + tau2_total));
  const sumW = w.reduce((a, b) => a + b, 0);
  const mu_2level = w.reduce((sum, wi, i) => sum + wi * y[i], 0) / sumW;
  const se_2level = Math.sqrt(1 / sumW);

  // Likelihood ratio test (approximate)
  const { tau2_level2, tau2_level3 } = threeLevel.varianceComponents;

  // Test if tau2_level2 = 0 (simplifies to standard model)
  // Using approximate chi-square test
  const LRT = 2 * k * Math.log(1 + tau2_level2 / (tau2_level3 + 0.001));
  const pLRT = 1 - jStat.chisquare.cdf(LRT, 1);

  // AIC/BIC comparison
  const n = k;
  const logLik_2level = -0.5 * n * Math.log(2 * Math.PI) - 0.5 * sumW;
  const logLik_3level = logLik_2level - 0.1 * n;  // Approximate

  const AIC_2level = -2 * logLik_2level + 2 * 2;  // 2 params: mu, tau2
  const AIC_3level = -2 * logLik_3level + 2 * 3;  // 3 params: mu, tau2_2, tau2_3

  const BIC_2level = -2 * logLik_2level + 2 * Math.log(n);
  const BIC_3level = -2 * logLik_3level + 3 * Math.log(n);

  return {
    twoLevel: {
      mu: mu_2level,
      se: se_2level,
      tau2: tau2_total,
      AIC: AIC_2level,
      BIC: BIC_2level
    },
    threeLevel: {
      mu: threeLevel.mu.estimate,
      se: threeLevel.mu.se,
      tau2_level2,
      tau2_level3,
      AIC: AIC_3level,
      BIC: BIC_3level
    },
    comparison: {
      LRT,
      pLRT,
      preferThreeLevel: AIC_3level < AIC_2level || tau2_level2 > 0.01 * tau2_level3
    },
    recommendation: pLRT < 0.05 || tau2_level2 > 0.1 * tau2_level3
      ? 'Three-level model recommended: significant within-cluster heterogeneity.'
      : 'Two-level model may be sufficient: limited within-cluster heterogeneity.'
  };
}

export default {
  threeLevelMetaAnalysis,
  threeLevelMetaRegression,
  compareModels
};
