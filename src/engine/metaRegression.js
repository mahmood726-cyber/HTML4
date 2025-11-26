/**
 * Meta-regression analysis
 * Explores sources of heterogeneity through covariate modeling
 * Reference: Thompson & Sharp (1999), Higgins & Thompson (2004)
 */

import { tCritical, pFromT, chiSqP } from './distributions.js';
import { tauDL, tauREML } from './tau.js';

/**
 * Univariate meta-regression with a single covariate
 * Uses weighted least squares with random effects
 *
 * Model: θᵢ = β₀ + β₁xᵢ + ζᵢ + εᵢ
 * where ζᵢ ~ N(0, τ²) and εᵢ ~ N(0, σᵢ²)
 *
 * @param {Array} effects - Array of effect objects with es, vi
 * @param {Array} covariate - Covariate values (same length as effects)
 * @param {Object} options - Analysis options
 * @returns {Object} Meta-regression results
 */
export function univariateMetaRegression(effects, covariate, options = {}) {
  const {
    tauMethod = 'REML',  // REML is preferred for meta-regression
    confLevel = 0.95,
    knha = true,  // Knapp-Hartung adjustment
    permutations = 0  // Number of permutation tests (0 = none)
  } = options;

  // Filter valid studies with covariate values
  const valid = effects
    .map((e, i) => ({ ...e, x: covariate[i], index: i }))
    .filter(e => !e.excluded && e.x !== null && e.x !== undefined && !isNaN(e.x));

  if (valid.length < 3) {
    return { error: 'Need at least 3 studies with valid covariate values' };
  }

  const k = valid.length;
  const y = valid.map(e => e.es);
  const x = valid.map(e => e.x);
  const vi = valid.map(e => e.vi);

  // Step 1: Estimate τ² using method of moments or REML
  // Initial estimate using DL on residuals
  let tau2;

  if (tauMethod === 'REML') {
    tau2 = estimateTau2REML(y, x, vi);
  } else {
    tau2 = estimateTau2MM(y, x, vi);
  }

  // Step 2: Calculate weights w = 1/(vi + τ²)
  const w = vi.map(v => 1 / (v + tau2));
  const sumW = w.reduce((a, b) => a + b, 0);

  // Step 3: Weighted least squares
  // Design matrix X = [1, x] for intercept and slope
  const sumWX = w.reduce((acc, wi, i) => acc + wi * x[i], 0);
  const sumWY = w.reduce((acc, wi, i) => acc + wi * y[i], 0);
  const sumWXY = w.reduce((acc, wi, i) => acc + wi * x[i] * y[i], 0);
  const sumWX2 = w.reduce((acc, wi, i) => acc + wi * x[i] * x[i], 0);

  // Normal equations: X'WX β = X'Wy
  // [sumW    sumWX  ] [β₀]   [sumWY ]
  // [sumWX   sumWX2 ] [β₁] = [sumWXY]

  const det = sumW * sumWX2 - sumWX * sumWX;
  if (Math.abs(det) < 1e-10) {
    return { error: 'Singular design matrix - check covariate for constant values' };
  }

  // Solve for coefficients
  const beta0 = (sumWX2 * sumWY - sumWX * sumWXY) / det;  // Intercept
  const beta1 = (sumW * sumWXY - sumWX * sumWY) / det;    // Slope

  // Step 4: Calculate fitted values and residuals
  const fitted = x.map(xi => beta0 + beta1 * xi);
  const residuals = y.map((yi, i) => yi - fitted[i]);

  // Step 5: Residual heterogeneity (QE)
  const QE = w.reduce((acc, wi, i) => acc + wi * residuals[i] ** 2, 0);
  const dfE = k - 2;  // df for residual heterogeneity
  const pQE = dfE > 0 ? chiSqP(QE, dfE) : 1;

  // I² residual
  const I2res = dfE > 0 ? Math.max(0, (QE - dfE) / QE * 100) : 0;

  // Step 6: Model test (QM) - tests if β₁ = 0
  // QM = β₁² / Var(β₁)
  const varBeta1 = sumW / det;
  const QM = beta1 ** 2 / varBeta1;
  const pQM = chiSqP(QM, 1);

  // Step 7: Standard errors and CIs
  let seBeta0 = Math.sqrt(sumWX2 / det);
  let seBeta1 = Math.sqrt(sumW / det);

  // Knapp-Hartung adjustment (recommended)
  if (knha && dfE > 0) {
    const khFactor = Math.max(1, QE / dfE);
    seBeta0 *= Math.sqrt(khFactor);
    seBeta1 *= Math.sqrt(khFactor);
  }

  const alpha = 1 - confLevel;
  const df = knha ? dfE : Infinity;
  const critVal = df > 0 && df < Infinity ? tCritical(df, alpha) : 1.96;

  // Test statistics and p-values
  const tBeta0 = beta0 / seBeta0;
  const tBeta1 = beta1 / seBeta1;
  const pBeta0 = df > 0 && df < Infinity ? pFromT(tBeta0, df) : 2 * (1 - normalCDF(Math.abs(tBeta0)));
  const pBeta1 = df > 0 && df < Infinity ? pFromT(tBeta1, df) : 2 * (1 - normalCDF(Math.abs(tBeta1)));

  // Confidence intervals
  const ci0 = {
    lo: beta0 - critVal * seBeta0,
    hi: beta0 + critVal * seBeta0
  };
  const ci1 = {
    lo: beta1 - critVal * seBeta1,
    hi: beta1 + critVal * seBeta1
  };

  // Step 8: R² - proportion of heterogeneity explained
  // R² = 1 - τ²_model / τ²_null
  const tau2Null = tauMethod === 'REML' ? tauREML(valid) : tauDL(valid);
  const R2 = tau2Null > 0 ? Math.max(0, 1 - tau2 / tau2Null) * 100 : 0;

  // Step 9: Permutation test (if requested)
  let permP = null;
  if (permutations > 0) {
    permP = permutationTest(y, x, vi, beta1, permutations, tauMethod);
  }

  // Prepare study-level results
  const studies = valid.map((e, i) => ({
    id: e.id,
    x: x[i],
    y: y[i],
    fitted: fitted[i],
    residual: residuals[i],
    weight: (w[i] / sumW) * 100
  }));

  return {
    // Coefficients
    intercept: {
      estimate: beta0,
      se: seBeta0,
      t: tBeta0,
      p: pBeta0,
      ci: ci0
    },
    slope: {
      estimate: beta1,
      se: seBeta1,
      t: tBeta1,
      p: pBeta1,
      ci: ci1,
      permP
    },

    // Model statistics
    k,
    tau2,
    tau2Null,
    R2,

    // Heterogeneity tests
    QM,       // Model test (covariate effect)
    pQM,
    dfM: 1,
    QE,       // Residual heterogeneity
    pQE,
    dfE,
    I2res,

    // Settings
    method: tauMethod,
    knha,
    confLevel,

    // Study data
    studies
  };
}

/**
 * Multivariate meta-regression with multiple covariates
 *
 * @param {Array} effects - Array of effect objects
 * @param {Array<Array>} covariates - 2D array of covariate values
 * @param {Array<string>} names - Covariate names
 * @param {Object} options - Analysis options
 * @returns {Object} Meta-regression results
 */
export function multivariateMetaRegression(effects, covariates, names = [], options = {}) {
  const {
    tauMethod = 'REML',
    confLevel = 0.95,
    knha = true
  } = options;

  // Number of covariates
  const p = covariates.length;
  if (p === 0) {
    return { error: 'No covariates provided' };
  }

  // Filter valid studies
  const valid = effects
    .map((e, i) => ({
      ...e,
      x: covariates.map(cov => cov[i]),
      index: i
    }))
    .filter(e => !e.excluded && e.x.every(v => v !== null && v !== undefined && !isNaN(v)));

  if (valid.length < p + 2) {
    return { error: `Need at least ${p + 2} studies for ${p} covariate(s)` };
  }

  const k = valid.length;
  const n = p + 1;  // Number of parameters (intercept + covariates)

  // Build design matrix X (k x n)
  // First column is 1 for intercept
  const X = valid.map(e => [1, ...e.x]);
  const y = valid.map(e => e.es);
  const vi = valid.map(e => e.vi);

  // Estimate τ² using iterative method
  let tau2 = estimateTau2MultiREML(y, X, vi, 100);

  // Calculate weights
  const w = vi.map(v => 1 / (v + tau2));
  const W = w;  // Diagonal weight matrix (as array)

  // WLS solution: β = (X'WX)⁻¹ X'Wy
  // Compute X'WX (n x n)
  const XtWX = [];
  for (let i = 0; i < n; i++) {
    XtWX[i] = [];
    for (let j = 0; j < n; j++) {
      let sum = 0;
      for (let s = 0; s < k; s++) {
        sum += X[s][i] * W[s] * X[s][j];
      }
      XtWX[i][j] = sum;
    }
  }

  // Compute X'Wy (n x 1)
  const XtWy = [];
  for (let i = 0; i < n; i++) {
    let sum = 0;
    for (let s = 0; s < k; s++) {
      sum += X[s][i] * W[s] * y[s];
    }
    XtWy[i] = sum;
  }

  // Invert X'WX
  const XtWXinv = invertMatrix(XtWX);
  if (!XtWXinv) {
    return { error: 'Singular design matrix - check for multicollinearity' };
  }

  // Calculate coefficients β = (X'WX)⁻¹ X'Wy
  const beta = [];
  for (let i = 0; i < n; i++) {
    let sum = 0;
    for (let j = 0; j < n; j++) {
      sum += XtWXinv[i][j] * XtWy[j];
    }
    beta[i] = sum;
  }

  // Fitted values and residuals
  const fitted = valid.map((e, i) => {
    return beta.reduce((sum, b, j) => sum + b * X[i][j], 0);
  });
  const residuals = y.map((yi, i) => yi - fitted[i]);

  // Residual heterogeneity QE
  const QE = W.reduce((acc, wi, i) => acc + wi * residuals[i] ** 2, 0);
  const dfE = k - n;

  // Knapp-Hartung adjustment
  let varScale = 1;
  if (knha && dfE > 0) {
    varScale = Math.max(1, QE / dfE);
  }

  // Standard errors from diagonal of (X'WX)⁻¹
  const se = XtWXinv.map((row, i) => Math.sqrt(row[i] * varScale));

  const alpha = 1 - confLevel;
  const critVal = knha && dfE > 0 ? tCritical(dfE, alpha) : 1.96;

  // Build coefficient results
  const coefficients = beta.map((b, i) => {
    const coefSE = se[i];
    const t = b / coefSE;
    const pVal = knha && dfE > 0 ? pFromT(t, dfE) : 2 * (1 - normalCDF(Math.abs(t)));

    return {
      name: i === 0 ? 'Intercept' : (names[i - 1] || `X${i}`),
      estimate: b,
      se: coefSE,
      t,
      p: pVal,
      ci: {
        lo: b - critVal * coefSE,
        hi: b + critVal * coefSE
      }
    };
  });

  // Model test QM - tests if all slopes = 0
  // QM = β'₁ Var(β₁)⁻¹ β₁ where β₁ excludes intercept
  let QM = 0;
  if (p > 0) {
    const beta1 = beta.slice(1);
    const cov1 = XtWXinv.slice(1).map(row => row.slice(1).map(v => v * varScale));
    const cov1inv = invertMatrix(cov1);
    if (cov1inv) {
      for (let i = 0; i < p; i++) {
        for (let j = 0; j < p; j++) {
          QM += beta1[i] * cov1inv[i][j] * beta1[j];
        }
      }
    }
  }
  const pQM = chiSqP(QM, p);

  const pQE = dfE > 0 ? chiSqP(QE, dfE) : 1;
  const I2res = dfE > 0 ? Math.max(0, (QE - dfE) / QE * 100) : 0;

  // R² calculation
  const tau2Null = tauMethod === 'REML' ? tauREML(valid) : tauDL(valid);
  const R2 = tau2Null > 0 ? Math.max(0, 1 - tau2 / tau2Null) * 100 : 0;

  return {
    coefficients,
    k,
    p,
    tau2,
    tau2Null,
    R2,
    QM,
    pQM,
    dfM: p,
    QE,
    pQE,
    dfE,
    I2res,
    method: tauMethod,
    knha,
    confLevel,
    studies: valid.map((e, i) => ({
      id: e.id,
      x: e.x,
      y: y[i],
      fitted: fitted[i],
      residual: residuals[i],
      weight: (W[i] / W.reduce((a, b) => a + b, 0)) * 100
    }))
  };
}

/**
 * Categorical subgroup analysis as meta-regression
 *
 * @param {Array} effects - Effect objects
 * @param {Array} groups - Group assignments (categorical)
 * @param {Object} options - Analysis options
 * @returns {Object} Subgroup comparison results
 */
export function subgroupMetaRegression(effects, groups, options = {}) {
  const {
    confLevel = 0.95,
    tauMethod = 'REML'
  } = options;

  // Get unique groups
  const uniqueGroups = [...new Set(groups.filter(g => g !== null && g !== undefined && g !== ''))];

  if (uniqueGroups.length < 2) {
    return { error: 'Need at least 2 groups for subgroup analysis' };
  }

  // Create dummy variables (reference = first group)
  const dummies = uniqueGroups.slice(1).map(g => groups.map(gi => gi === g ? 1 : 0));

  const result = multivariateMetaRegression(
    effects,
    dummies,
    uniqueGroups.slice(1),
    { confLevel, tauMethod, knha: true }
  );

  if (result.error) return result;

  // Add group-specific pooled estimates
  const groupResults = {};

  uniqueGroups.forEach(g => {
    const groupEffects = effects.filter((e, i) => !e.excluded && groups[i] === g);
    if (groupEffects.length > 0) {
      const tau2g = tauMethod === 'REML' ? tauREML(groupEffects) : tauDL(groupEffects);
      const wg = groupEffects.map(e => 1 / (e.vi + tau2g));
      const sumWg = wg.reduce((a, b) => a + b, 0);
      const esg = groupEffects.reduce((acc, e, i) => acc + wg[i] * e.es, 0) / sumWg;
      const seg = Math.sqrt(1 / sumWg);

      groupResults[g] = {
        k: groupEffects.length,
        es: esg,
        se: seg,
        ci: {
          lo: esg - 1.96 * seg,
          hi: esg + 1.96 * seg
        },
        tau2: tau2g
      };
    }
  });

  return {
    ...result,
    type: 'subgroup',
    groups: uniqueGroups,
    groupResults,
    referenceGroup: uniqueGroups[0]
  };
}

// ============ Helper Functions ============

/**
 * Estimate τ² for meta-regression using method of moments
 */
function estimateTau2MM(y, x, vi) {
  const k = y.length;
  const w = vi.map(v => 1 / v);
  const sumW = w.reduce((a, b) => a + b, 0);

  // Calculate WLS residuals assuming τ² = 0
  const sumWX = w.reduce((acc, wi, i) => acc + wi * x[i], 0);
  const sumWY = w.reduce((acc, wi, i) => acc + wi * y[i], 0);
  const sumWXY = w.reduce((acc, wi, i) => acc + wi * x[i] * y[i], 0);
  const sumWX2 = w.reduce((acc, wi, i) => acc + wi * x[i] * x[i], 0);

  const det = sumW * sumWX2 - sumWX * sumWX;
  if (Math.abs(det) < 1e-10) return 0;

  const beta0 = (sumWX2 * sumWY - sumWX * sumWXY) / det;
  const beta1 = (sumW * sumWXY - sumWX * sumWY) / det;

  // Q statistic for residuals
  const Q = w.reduce((acc, wi, i) => {
    const resid = y[i] - beta0 - beta1 * x[i];
    return acc + wi * resid * resid;
  }, 0);

  // Moment estimate
  const C = sumW - (sumWX2 * sumW * sumW + sumW * sumWX * sumWX - 2 * sumWX * sumW * sumWX) / (sumW * det);

  const tau2 = Math.max(0, (Q - (k - 2)) / C);
  return tau2;
}

/**
 * Estimate τ² for meta-regression using REML
 */
function estimateTau2REML(y, x, vi, maxIter = 50) {
  let tau2 = estimateTau2MM(y, x, vi);
  const k = y.length;

  for (let iter = 0; iter < maxIter; iter++) {
    const w = vi.map(v => 1 / (v + tau2));
    const sumW = w.reduce((a, b) => a + b, 0);
    const sumWX = w.reduce((acc, wi, i) => acc + wi * x[i], 0);
    const sumWY = w.reduce((acc, wi, i) => acc + wi * y[i], 0);
    const sumWXY = w.reduce((acc, wi, i) => acc + wi * x[i] * y[i], 0);
    const sumWX2 = w.reduce((acc, wi, i) => acc + wi * x[i] * x[i], 0);

    const det = sumW * sumWX2 - sumWX * sumWX;
    if (Math.abs(det) < 1e-10) break;

    const beta0 = (sumWX2 * sumWY - sumWX * sumWXY) / det;
    const beta1 = (sumW * sumWXY - sumWX * sumWY) / det;

    // REML score and Fisher information
    const w2 = w.map(wi => wi * wi);
    const sumW2 = w2.reduce((a, b) => a + b, 0);

    // Simplified REML update using Q statistic
    let sumW2Resid2 = 0;
    for (let i = 0; i < k; i++) {
      const resid = y[i] - beta0 - beta1 * x[i];
      sumW2Resid2 += w2[i] * resid * resid;
    }

    // Simplified update
    const Q = w.reduce((acc, wi, i) => {
      const resid = y[i] - beta0 - beta1 * x[i];
      return acc + wi * resid * resid;
    }, 0);

    const tau2New = Math.max(0, tau2 + (Q - (k - 2)) / sumW2);

    if (Math.abs(tau2New - tau2) < 1e-6) break;
    tau2 = tau2New;
  }

  return tau2;
}

/**
 * Estimate τ² for multivariate meta-regression using REML
 */
function estimateTau2MultiREML(y, X, vi, maxIter = 50) {
  const k = y.length;
  const p = X[0].length;
  let tau2 = tauDL(y.map((yi, i) => ({ es: yi, vi: vi[i] })));

  for (let iter = 0; iter < maxIter; iter++) {
    const w = vi.map(v => 1 / (v + tau2));

    // Compute X'WX and X'Wy
    const XtWX = [];
    for (let i = 0; i < p; i++) {
      XtWX[i] = [];
      for (let j = 0; j < p; j++) {
        let sum = 0;
        for (let s = 0; s < k; s++) {
          sum += X[s][i] * w[s] * X[s][j];
        }
        XtWX[i][j] = sum;
      }
    }

    const XtWy = [];
    for (let i = 0; i < p; i++) {
      let sum = 0;
      for (let s = 0; s < k; s++) {
        sum += X[s][i] * w[s] * y[s];
      }
      XtWy[i] = sum;
    }

    const XtWXinv = invertMatrix(XtWX);
    if (!XtWXinv) break;

    const beta = [];
    for (let i = 0; i < p; i++) {
      let sum = 0;
      for (let j = 0; j < p; j++) {
        sum += XtWXinv[i][j] * XtWy[j];
      }
      beta[i] = sum;
    }

    // Residual Q
    const Q = w.reduce((acc, wi, i) => {
      const fitted = beta.reduce((sum, b, j) => sum + b * X[i][j], 0);
      return acc + wi * (y[i] - fitted) ** 2;
    }, 0);

    const sumW2 = w.reduce((a, b) => a + b * b, 0);
    const tau2New = Math.max(0, tau2 + (Q - (k - p)) / sumW2);

    if (Math.abs(tau2New - tau2) < 1e-6) break;
    tau2 = tau2New;
  }

  return tau2;
}

/**
 * Matrix inversion using Gauss-Jordan elimination
 */
function invertMatrix(matrix) {
  const n = matrix.length;
  const aug = matrix.map((row, i) => [...row, ...Array(n).fill(0).map((_, j) => i === j ? 1 : 0)]);

  for (let i = 0; i < n; i++) {
    // Find pivot
    let maxRow = i;
    for (let k = i + 1; k < n; k++) {
      if (Math.abs(aug[k][i]) > Math.abs(aug[maxRow][i])) maxRow = k;
    }
    [aug[i], aug[maxRow]] = [aug[maxRow], aug[i]];

    if (Math.abs(aug[i][i]) < 1e-10) return null;  // Singular

    // Scale pivot row
    const pivot = aug[i][i];
    for (let j = 0; j < 2 * n; j++) aug[i][j] /= pivot;

    // Eliminate column
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
 * Permutation test for meta-regression slope
 */
function permutationTest(y, x, vi, observedSlope, nPerm, tauMethod) {
  let count = 0;
  const absObserved = Math.abs(observedSlope);

  for (let p = 0; p < nPerm; p++) {
    // Shuffle y values
    const shuffledY = [...y].sort(() => Math.random() - 0.5);

    // Fit model with shuffled data
    const tau2 = tauMethod === 'REML'
      ? estimateTau2REML(shuffledY, x, vi, 20)
      : estimateTau2MM(shuffledY, x, vi);

    const w = vi.map(v => 1 / (v + tau2));
    const sumW = w.reduce((a, b) => a + b, 0);
    const sumWX = w.reduce((acc, wi, i) => acc + wi * x[i], 0);
    const sumWY = w.reduce((acc, wi, i) => acc + wi * shuffledY[i], 0);
    const sumWXY = w.reduce((acc, wi, i) => acc + wi * x[i] * shuffledY[i], 0);
    const sumWX2 = w.reduce((acc, wi, i) => acc + wi * x[i] * x[i], 0);

    const det = sumW * sumWX2 - sumWX * sumWX;
    if (Math.abs(det) < 1e-10) continue;

    const permSlope = (sumW * sumWXY - sumWX * sumWY) / det;

    if (Math.abs(permSlope) >= absObserved) count++;
  }

  return (count + 1) / (nPerm + 1);
}

/**
 * Normal CDF
 */
function normalCDF(x) {
  const a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741;
  const a4 = -1.453152027, a5 = 1.061405429, p = 0.3275911;

  const sign = x < 0 ? -1 : 1;
  x = Math.abs(x) / Math.sqrt(2);
  const t = 1 / (1 + p * x);
  const y = 1 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-x * x);

  return 0.5 * (1 + sign * y);
}

/**
 * Create bubble plot data for meta-regression
 * @param {Object} regResult - Meta-regression result
 * @returns {Object} Plot data for visualization
 */
export function createBubblePlotData(regResult) {
  if (!regResult || !regResult.studies) return null;

  const { studies, intercept, slope } = regResult;

  // Calculate regression line
  const xValues = studies.map(s => s.x);
  const xMin = Math.min(...xValues);
  const xMax = Math.max(...xValues);
  const xRange = xMax - xMin;

  const lineX = [xMin - 0.1 * xRange, xMax + 0.1 * xRange];
  const lineY = lineX.map(x => intercept.estimate + slope.estimate * x);

  // CI band (approximate)
  const ciY = lineX.map(x => {
    const pred = intercept.estimate + slope.estimate * x;
    const xDev = x - xValues.reduce((a, b) => a + b, 0) / xValues.length;
    const se = Math.sqrt(intercept.se ** 2 + (slope.se * xDev) ** 2);
    return {
      lo: pred - 1.96 * se,
      hi: pred + 1.96 * se
    };
  });

  return {
    points: studies.map(s => ({
      x: s.x,
      y: s.y,
      size: s.weight,
      id: s.id
    })),
    line: { x: lineX, y: lineY },
    ci: { x: lineX, lo: ciY.map(c => c.lo), hi: ciY.map(c => c.hi) }
  };
}
