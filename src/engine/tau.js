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

  for (let i = 0; i < maxIter; i++) {
    const w = effects.map(e => 1 / (e.vi + tau2));
    const sumW = w.reduce((a, b) => a + b, 0);

    // Weighted mean
    const mu = effects.reduce((acc, e, j) => acc + w[j] * e.es, 0) / sumW;

    // Fisher scoring update
    const num = effects.reduce((acc, e, j) => acc + w[j] ** 2 * ((e.es - mu) ** 2 - e.vi), 0);
    const den = effects.reduce((acc, e, j) => acc + w[j] ** 2, 0);

    const tau2New = Math.max(0, tau2 + num / den);

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

    // Newton-Raphson update
    tau2 = Math.max(0, tau2 - (Q - (k - 1)) / dQ);
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
 * @param {Array} effects - Array of effect objects
 * @returns {number} Tau² estimate
 */
export function tauHS(effects) {
  const k = effects.length;
  if (k < 2) return 0;

  // Sample-size weighted mean (approximated by inverse-variance for simplicity)
  const w = effects.map(e => 1 / e.vi);
  const sumW = w.reduce((a, b) => a + b, 0);
  const mu = effects.reduce((acc, e, i) => acc + w[i] * e.es, 0) / sumW;

  // Weighted variance
  const weightedVar = effects.reduce((acc, e, i) => acc + w[i] * (e.es - mu) ** 2, 0) / sumW;

  // Average sampling variance (weighted)
  const avgVi = effects.reduce((acc, e, i) => acc + w[i] * e.vi, 0) / sumW;

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
