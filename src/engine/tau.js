/**
 * Tau² (between-study variance) estimators
 */

import { REML_CONFIG, PM_CONFIG } from '../constants/index.js';

/**
 * DerSimonian-Laird tau² estimator
 * @param {Array} effects - Array of effect objects with es and vi properties
 * @returns {number} Tau² estimate
 */
export function tauDL(effects) {
  const k = effects.length;
  if (k < 2) return 0;

  const w = effects.map(e => 1 / e.vi);
  const sumW = w.reduce((a, b) => a + b, 0);
  const sumW2 = w.reduce((a, b) => a + b * b, 0);

  // Fixed-effect weighted mean
  const mu = effects.reduce((acc, e, i) => acc + w[i] * e.es, 0) / sumW;

  // Cochran's Q
  const Q = effects.reduce((acc, e, i) => acc + w[i] * (e.es - mu) ** 2, 0);

  // C constant
  const C = sumW - sumW2 / sumW;

  // DL estimator
  return Math.max(0, (Q - (k - 1)) / C);
}

/**
 * REML (Restricted Maximum Likelihood) tau² estimator
 * Uses Fisher scoring algorithm with proper REML score function
 * @param {Array} effects - Array of effect objects
 * @param {number} maxIter - Maximum iterations (default from config)
 * @param {number} tol - Convergence tolerance (default from config)
 * @returns {number} Tau² estimate
 */
export function tauREML(effects, maxIter = REML_CONFIG.maxIterations, tol = REML_CONFIG.tolerance) {
  const k = effects.length;
  if (k < 2) return 0;

  // Initialize with DL estimate
  let tau2 = tauDL(effects);

  // Upper bound to prevent divergence
  const maxTau2 = effects.reduce((acc, e) => acc + e.es ** 2, 0);

  for (let iter = 0; iter < maxIter; iter++) {
    const w = effects.map(e => 1 / (e.vi + tau2));
    const sumW = w.reduce((a, b) => a + b, 0);
    const sumW2 = w.reduce((a, b) => a + b * b, 0);

    // Weighted mean
    const mu = effects.reduce((acc, e, j) => acc + w[j] * e.es, 0) / sumW;

    // REML score function: S(τ²) = -0.5 * Σwi + 0.5 * Σwi²(yi-μ)²
    // Fisher information: I(τ²) = 0.5 * Σwi²
    // Newton-Raphson: τ² += S / I = [Σwi²(yi-μ)² - Σwi] / Σwi²
    const sumW2Resid2 = effects.reduce((acc, e, j) => acc + w[j] * w[j] * (e.es - mu) ** 2, 0);

    const numerator = sumW2Resid2 - sumW;
    const denominator = sumW2;

    if (Math.abs(denominator) < 1e-10) {
      return tau2;
    }

    const tau2New = Math.max(0, Math.min(maxTau2, tau2 + numerator / denominator));

    // Check convergence
    if (Math.abs(tau2New - tau2) < tol) {
      return tau2New;
    }

    tau2 = tau2New;
  }

  return tau2;
}

/**
 * Paule-Mandel tau² estimator
 * @param {Array} effects - Array of effect objects
 * @param {number} maxIter - Maximum iterations
 * @param {number} tol - Convergence tolerance
 * @returns {number} Tau² estimate
 */
export function tauPM(effects, maxIter = PM_CONFIG.maxIterations, tol = PM_CONFIG.tolerance) {
  const k = effects.length;
  if (k < 2) return 0;

  // Initialize with DL estimate
  let tau2 = tauDL(effects);

  // If DL gives 0, check if data is homogeneous
  if (tau2 === 0) {
    // Check Q with fixed-effect weights
    const wFE = effects.map(e => 1 / e.vi);
    const sumWFE = wFE.reduce((a, b) => a + b, 0);
    const muFE = effects.reduce((acc, e, j) => acc + wFE[j] * e.es, 0) / sumWFE;
    const QFE = effects.reduce((acc, e, j) => acc + wFE[j] * (e.es - muFE) ** 2, 0);

    // If Q < k-1, data is homogeneous, return 0
    if (QFE <= k - 1) {
      return 0;
    }
  }

  // Upper bound for tau² to prevent divergence
  const esRange = Math.max(...effects.map(e => e.es)) - Math.min(...effects.map(e => e.es));
  const maxTau2 = esRange * esRange * 10; // Reasonable upper bound

  for (let i = 0; i < maxIter; i++) {
    const w = effects.map(e => 1 / (e.vi + tau2));
    const sumW = w.reduce((a, b) => a + b, 0);

    // Weighted mean
    const mu = effects.reduce((acc, e, j) => acc + w[j] * e.es, 0) / sumW;

    // Q statistic with current weights
    const Q = effects.reduce((acc, e, j) => acc + w[j] * (e.es - mu) ** 2, 0);

    // Check if Q ≈ k-1 (convergence criterion)
    if (Math.abs(Q - (k - 1)) < tol) {
      return tau2;
    }

    // Derivative of Q with respect to tau²
    const dQ = -effects.reduce((acc, e, j) => {
      return acc + w[j] ** 2 * ((e.es - mu) ** 2 - (1 - w[j] / sumW) / w[j]);
    }, 0);

    if (Math.abs(dQ) < tol) {
      return tau2;
    }

    // Newton-Raphson update with bounds
    const delta = (Q - (k - 1)) / dQ;
    tau2 = Math.max(0, Math.min(maxTau2, tau2 - delta));
  }

  return tau2;
}

/**
 * Sidik-Jonkman tau² estimator
 * @param {Array} effects - Array of effect objects
 * @returns {number} Tau² estimate
 */
export function tauSJ(effects) {
  const k = effects.length;
  if (k < 2) return 0;

  // Unweighted mean
  const mu0 = effects.reduce((a, e) => a + e.es, 0) / k;

  // Initial variance estimate (sample variance of effects)
  const tau2_0 = effects.reduce((a, e) => a + (e.es - mu0) ** 2, 0) / (k - 1);

  // Weights using initial variance
  const w = effects.map(e => 1 / (e.vi + tau2_0));
  const sumW = w.reduce((a, b) => a + b, 0);

  // Weighted mean
  const mu = effects.reduce((acc, e, j) => acc + w[j] * e.es, 0) / sumW;

  // Q statistic
  const Q = effects.reduce((acc, e, j) => acc + w[j] * (e.es - mu) ** 2, 0);

  // SJ estimator
  return Math.max(0, (Q - (k - 1)) * tau2_0 / Q);
}

/**
 * Hedges tau² estimator
 * @param {Array} effects - Array of effect objects
 * @returns {number} Tau² estimate
 */
export function tauHE(effects) {
  const k = effects.length;
  if (k < 2) return 0;

  // Unweighted mean
  const mu = effects.reduce((a, e) => a + e.es, 0) / k;

  // Unweighted variance of effects
  const S2 = effects.reduce((a, e) => a + (e.es - mu) ** 2, 0) / (k - 1);

  // Average sampling variance
  const avgVi = effects.reduce((a, e) => a + e.vi, 0) / k;

  return Math.max(0, S2 - avgVi);
}

/**
 * Hunter-Schmidt tau² estimator
 * Uses sample sizes as weights (not inverse variance)
 * @param {Array} effects - Array of effect objects
 * @returns {number} Tau² estimate
 */
export function tauHS(effects) {
  const k = effects.length;
  if (k < 2) return 0;

  // Get sample sizes - use raw data if available, otherwise estimate from variance
  const n = effects.map(e => {
    if (e.raw) {
      // Binary/continuous data
      if (e.raw.n1 !== undefined && e.raw.n2 !== undefined) {
        return e.raw.n1 + e.raw.n2;
      }
      // Single group data
      if (e.raw.n !== undefined) {
        return e.raw.n;
      }
    }
    // Fallback: estimate N from variance (rough approximation)
    // For SMD, vi ≈ 2/N, so N ≈ 2/vi
    return Math.max(10, Math.round(2 / e.vi));
  });

  const sumN = n.reduce((a, b) => a + b, 0);

  // Sample-size weighted mean
  const mu = effects.reduce((acc, e, i) => acc + n[i] * e.es, 0) / sumN;

  // Weighted variance of effects (between-study variance + sampling variance)
  const weightedVar = effects.reduce((acc, e, i) => acc + n[i] * (e.es - mu) ** 2, 0) / sumN;

  // Sample-size weighted average sampling variance
  const avgVi = effects.reduce((acc, e, i) => acc + n[i] * e.vi, 0) / sumN;

  return Math.max(0, weightedVar - avgVi);
}

/**
 * Empirical Bayes tau² estimator
 * @param {Array} effects - Array of effect objects
 * @param {number} maxIter - Maximum iterations
 * @param {number} tol - Convergence tolerance
 * @returns {number} Tau² estimate
 */
export function tauEB(effects, maxIter = REML_CONFIG.maxIterations, tol = REML_CONFIG.tolerance) {
  const k = effects.length;
  if (k < 2) return 0;

  // Initialize with method of moments
  let tau2 = tauHE(effects);

  for (let i = 0; i < maxIter; i++) {
    const w = effects.map(e => 1 / (e.vi + tau2));
    const sumW = w.reduce((a, b) => a + b, 0);

    // Weighted mean
    const mu = effects.reduce((acc, e, j) => acc + w[j] * e.es, 0) / sumW;

    // EB update
    const num = effects.reduce((acc, e, j) => {
      return acc + (e.es - mu) ** 2 - e.vi;
    }, 0);

    const tau2New = Math.max(0, num / k);

    if (Math.abs(tau2New - tau2) < tol) {
      return tau2New;
    }

    tau2 = tau2New;
  }

  return tau2;
}

/**
 * Get tau² estimator function by method name
 * @param {string} method - Estimator method name
 * @returns {Function} Tau² estimator function
 */
export function getTauEstimator(method) {
  const estimators = {
    dl: tauDL,
    reml: tauREML,
    pm: tauPM,
    sj: tauSJ,
    he: tauHE,
    hs: tauHS,
    eb: tauEB
  };

  return estimators[method] || tauDL;
}
