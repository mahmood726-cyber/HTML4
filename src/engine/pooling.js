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

  // Tau and its CI (using Q-profile method approximation)
  const tau = Math.sqrt(tau2);
  const tauCI = calculateTauCI(tau2, Q, df, k, confLevel);

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
 * Calculate CI for tau using Q-profile method approximation
 * @param {number} tau2 - Tau² estimate
 * @param {number} Q - Q statistic
 * @param {number} df - Degrees of freedom
 * @param {number} k - Number of studies
 * @param {number} confLevel - Confidence level
 * @returns {Object} CI bounds { lo, hi }
 */
function calculateTauCI(tau2, Q, df, k, confLevel) {
  if (df <= 0 || tau2 === 0) {
    return { lo: 0, hi: 0 };
  }

  // Simplified approximation using chi-square distribution
  // This is a rough approximation; full Q-profile requires iteration
  const alpha = 1 - confLevel;
  const mult = 1.5; // Rough multiplier for CI width

  const lo = Math.max(0, Math.sqrt(tau2) * (1 - mult / Math.sqrt(k)));
  const hi = Math.sqrt(tau2) * (1 + mult / Math.sqrt(k));

  return { lo, hi };
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
