/**
 * Meta-analysis pooling functions
 */

import { tCritical, pFromZ, pFromT, chiSqP } from './distributions.js';
import { DEFAULT_CONF_LEVEL } from '../constants/index.js';

/**
 * Pool effects using inverse variance method
 * @param {Array} effects - Array of effect objects
 * @param {number} tau2 - Between-study variance
 * @param {Object} options - Pooling options
 * @returns {Object|null} Pooled results
 */
export function poolInverseVariance(effects, tau2, options = {}) {
  const {
    useHKSJ = true,
    confLevel = DEFAULT_CONF_LEVEL,
    showPredictionInterval = true
  } = options;

  const active = effects.filter(e => !e.excluded);
  if (active.length === 0) return null;

  const k = active.length;
  const df = k - 1;
  const alpha = 1 - confLevel;

  // Calculate weights
  const w = active.map(e => 1 / (e.vi + tau2));
  const sumW = w.reduce((a, b) => a + b, 0);

  // Pooled effect size
  const es = active.reduce((acc, e, i) => acc + w[i] * e.es, 0) / sumW;

  // Standard error
  let se = Math.sqrt(1 / sumW);

  // Cochran's Q
  const Q = active.reduce((acc, e, i) => acc + w[i] * (e.es - es) ** 2, 0);

  // I² statistic
  const I2 = df > 0 ? Math.max(0, (Q - df) / Q * 100) : 0;

  // H² statistic
  const H2 = df > 0 ? Q / df : 1;

  // Critical value for CI
  let critVal = 1.96;

  // Apply HKSJ adjustment if requested and df > 0
  if (useHKSJ && df > 0) {
    const qStar = Q / df;
    se = se * Math.sqrt(Math.max(1, qStar));
    critVal = tCritical(df, alpha);
  }

  // Test statistic and p-value
  const z = es / se;
  const pVal = useHKSJ && df > 0 ? pFromT(z, df) : pFromZ(z);

  // Confidence interval
  const ciLo = es - critVal * se;
  const ciHi = es + critVal * se;

  // Prediction interval (uses k-2 df per Higgins et al. 2009)
  let piLo = null;
  let piHi = null;
  if (showPredictionInterval && k > 2) {
    const piSe = Math.sqrt(se * se + tau2);
    const piCrit = tCritical(k - 2, alpha);
    piLo = es - piCrit * piSe;
    piHi = es + piCrit * piSe;
  }

  // Calculate study-level statistics
  const studyStats = calculateStudyStats(effects, es, tau2, sumW);

  // Calculate total sample size
  const totalN = active.reduce((acc, s) => {
    if (s.raw) {
      return acc + (s.raw.n1 || 0) + (s.raw.n2 || 0) + (s.raw.n || 0);
    }
    return acc;
  }, 0);

  // Q test p-value
  const qPVal = chiSqP(Q, df);

  // Tau and its CI (using Q-profile method)
  const tau = Math.sqrt(tau2);
  const tauCI = calculateTauCI(tau2, Q, df, k, confLevel, active);

  // H² CI
  const H2CI = calculateH2CI(H2, df, confLevel);

  return {
    es,
    se,
    z,
    pVal,
    ciLo,
    ciHi,
    piLo,
    piHi,
    Q,
    qPVal,
    df,
    I2,
    H2,
    H2CI,
    tau2,
    tau,
    tauCI,
    k,
    totalN,
    studies: studyStats
  };
}

/**
 * Calculate study-level statistics (weights, contributions)
 * @param {Array} effects - All effects
 * @param {number} pooledES - Pooled effect size
 * @param {number} tau2 - Between-study variance
 * @param {number} sumW - Sum of weights
 * @returns {Array} Effects with added statistics
 */
function calculateStudyStats(effects, pooledES, tau2, sumW) {
  return effects.map(e => {
    if (e.excluded) {
      return { ...e, w: 0, q: 0, imp: 0 };
    }

    const wi = 1 / (e.vi + tau2);
    const w = (wi / sumW) * 100; // Weight as percentage
    const q = wi * (e.es - pooledES) ** 2; // Contribution to Q

    // Influence on pooled estimate (leave-one-out difference)
    const sumWLoo = sumW - wi;
    const esLoo = sumWLoo > 0
      ? (pooledES * sumW - wi * e.es) / sumWLoo
      : pooledES;
    const imp = Math.abs(pooledES - esLoo);

    return { ...e, w, q, imp };
  });
}

/**
 * Calculate CI for tau² using Q-profile method
 * Finds tau² values where Q(tau²) = chi-square critical values
 * @param {number} tau2 - Tau² estimate
 * @param {number} Q - Q statistic
 * @param {number} df - Degrees of freedom
 * @param {number} k - Number of studies
 * @param {number} confLevel - Confidence level
 * @param {Array} effects - Effect objects for Q-profile calculation
 * @returns {Object} CI bounds { lo, hi }
 */
function calculateTauCI(tau2, Q, df, k, confLevel, effects = null) {
  if (df <= 0) {
    return { lo: 0, hi: 0 };
  }

  const alpha = 1 - confLevel;

  // Chi-square critical values for Q-profile bounds
  // Q ~ chi-square(k-1) under null, find tau² where Q equals critical values
  const chiLo = chiSquareQuantile(alpha / 2, df);
  const chiHi = chiSquareQuantile(1 - alpha / 2, df);

  // If we have effects, use iterative Q-profile method
  if (effects && effects.length >= 2) {
    const active = effects.filter(e => !e.excluded);
    if (active.length >= 2) {
      const lo = qProfileSearch(active, chiHi, 0, tau2 * 10 + 1);
      const hi = qProfileSearch(active, chiLo, tau2, tau2 * 20 + 10);
      return { lo: Math.sqrt(Math.max(0, lo)), hi: Math.sqrt(hi) };
    }
  }

  // Fallback: Jackson's approximation when effects not available
  // Based on Q ~ chi-square relationship
  if (Q <= df) {
    return { lo: 0, hi: 0 };
  }

  // Approximate using scaled relationship
  const C = df; // Approximation of the C constant
  const tau2Lo = Math.max(0, (Q - chiHi) / C);
  const tau2Hi = Math.max(tau2, (Q - chiLo) / C);

  return { lo: Math.sqrt(tau2Lo), hi: Math.sqrt(tau2Hi) };
}

/**
 * Q-profile search: find tau² where Q(tau²) = target
 * @param {Array} effects - Effect objects
 * @param {number} target - Target Q value
 * @param {number} lower - Lower search bound
 * @param {number} upper - Upper search bound
 * @returns {number} Tau² value
 */
function qProfileSearch(effects, target, lower, upper) {
  const maxIter = 50;
  const tol = 1e-6;

  for (let i = 0; i < maxIter; i++) {
    const mid = (lower + upper) / 2;
    const Q = calculateQ(effects, mid);

    if (Math.abs(Q - target) < tol || upper - lower < tol) {
      return mid;
    }

    // Q decreases as tau² increases
    if (Q > target) {
      lower = mid;
    } else {
      upper = mid;
    }
  }

  return (lower + upper) / 2;
}

/**
 * Calculate Q statistic for given tau²
 * @param {Array} effects - Effect objects
 * @param {number} tau2 - Between-study variance
 * @returns {number} Q statistic
 */
function calculateQ(effects, tau2) {
  const w = effects.map(e => 1 / (e.vi + tau2));
  const sumW = w.reduce((a, b) => a + b, 0);
  const mu = effects.reduce((acc, e, i) => acc + w[i] * e.es, 0) / sumW;
  return effects.reduce((acc, e, i) => acc + w[i] * (e.es - mu) ** 2, 0);
}

/**
 * Chi-square quantile function (Wilson-Hilferty approximation inverse)
 * @param {number} p - Probability
 * @param {number} df - Degrees of freedom
 * @returns {number} Chi-square quantile
 */
function chiSquareQuantile(p, df) {
  if (df <= 0 || p <= 0 || p >= 1) return df;

  // Wilson-Hilferty approximation inverse
  const z = normQuantileApprox(p);
  const term = 1 - 2 / (9 * df) + z * Math.sqrt(2 / (9 * df));
  return df * Math.pow(Math.max(0.001, term), 3);
}

/**
 * Simple normal quantile approximation
 * @param {number} p - Probability
 * @returns {number} Z-score
 */
function normQuantileApprox(p) {
  if (p <= 0) return -8;
  if (p >= 1) return 8;
  if (p === 0.5) return 0;

  // Rational approximation (Abramowitz & Stegun 26.2.23)
  const t = p < 0.5 ? Math.sqrt(-2 * Math.log(p)) : Math.sqrt(-2 * Math.log(1 - p));
  const c0 = 2.515517, c1 = 0.802853, c2 = 0.010328;
  const d1 = 1.432788, d2 = 0.189269, d3 = 0.001308;

  const z = t - (c0 + c1 * t + c2 * t * t) / (1 + d1 * t + d2 * t * t + d3 * t * t * t);
  return p < 0.5 ? -z : z;
}

/**
 * Calculate CI for H² statistic
 * @param {number} H2 - H² statistic
 * @param {number} df - Degrees of freedom
 * @param {number} confLevel - Confidence level
 * @returns {Object} CI bounds { lo, hi }
 */
function calculateH2CI(H2, df, confLevel) {
  if (df <= 0) {
    return { lo: 1, hi: 1 };
  }

  // Use log transformation for CI
  const lnH = Math.log(Math.sqrt(Math.max(1, H2)));
  const seLnH = 0.5 / Math.sqrt(df);
  const z = 1.96; // For 95% CI

  const lo = Math.exp(lnH - z * seLnH) ** 2;
  const hi = Math.exp(lnH + z * seLnH) ** 2;

  return { lo: Math.max(1, lo), hi };
}

/**
 * Perform leave-one-out sensitivity analysis
 * @param {Array} effects - Array of effect objects
 * @param {Function} tauEstimator - Tau² estimator function
 * @param {Object} options - Pooling options
 * @returns {Array} Leave-one-out results
 */
export function leaveOneOut(effects, tauEstimator, options = {}) {
  const active = effects.filter(e => !e.excluded);
  if (active.length < 3) return [];

  return active.map((excluded, i) => {
    const subset = active.filter((_, j) => j !== i);
    const tau2 = tauEstimator(subset);
    const pooled = poolInverseVariance(
      subset.map(s => ({ ...s })),
      tau2,
      { ...options, useHKSJ: false }
    );

    return {
      excludedStudy: excluded.id,
      es: pooled?.es,
      ciLo: pooled?.ciLo,
      ciHi: pooled?.ciHi,
      I2: pooled?.I2
    };
  });
}

/**
 * Perform cumulative meta-analysis
 * @param {Array} effects - Array of effect objects (will be sorted by year)
 * @param {Function} tauEstimator - Tau² estimator function
 * @param {Object} options - Pooling options
 * @returns {Array} Cumulative results
 */
export function cumulativeMetaAnalysis(effects, tauEstimator, options = {}) {
  const active = effects.filter(e => !e.excluded);
  if (active.length < 2) return [];

  // Sort by year (if available)
  const sorted = [...active].sort((a, b) => {
    const yearA = parseInt(a.year) || 9999;
    const yearB = parseInt(b.year) || 9999;
    return yearA - yearB;
  });

  const results = [];

  for (let i = 1; i <= sorted.length; i++) {
    const subset = sorted.slice(0, i);
    const tau2 = tauEstimator(subset);
    const pooled = poolInverseVariance(
      subset.map(s => ({ ...s })),
      tau2,
      { ...options, useHKSJ: false }
    );

    if (pooled) {
      results.push({
        k: i,
        study: sorted[i - 1].id,
        year: sorted[i - 1].year,
        es: pooled.es,
        ciLo: pooled.ciLo,
        ciHi: pooled.ciHi,
        I2: pooled.I2
      });
    }
  }

  return results;
}

/**
 * Calculate influence diagnostics
 * @param {Object} pooledResult - Pooled meta-analysis result
 * @returns {Array} Influence diagnostics for each study
 */
export function calculateInfluenceDiagnostics(pooledResult) {
  const active = pooledResult.studies.filter(s => !s.excluded);
  if (active.length < 3) return [];

  return active.map(s => {
    // Standardized residual
    const resid = (s.es - pooledResult.es) / Math.sqrt(s.vi + pooledResult.tau2);

    // DFFITS approximation
    const dffits = resid * Math.sqrt(s.w / 100);

    // Cook's D approximation
    const leverage = s.w / 100;
    const cookD = leverage > 1 ? 0 : resid ** 2 * leverage / (1 - leverage);

    // Hat value (leverage)
    const hat = leverage;

    return {
      id: s.id,
      residual: resid,
      dffits,
      cookD,
      hat,
      contribution: s.q,
      weight: s.w
    };
  });
}
