/**
 * Effect size calculations for different data types
 */

import { parseNumber } from '../utils/index.js';
import { DEFAULT_CONTINUITY_CORRECTION } from '../constants/index.js';

/**
 * Calculate effect size from study data
 * @param {Object} row - Study data row
 * @param {string} type - Data type (binary, continuous, etc.)
 * @param {string} metric - Effect measure (OR, RR, RD, etc.)
 * @param {number} cc - Continuity correction (default: 0.5)
 * @returns {Object|null} Effect size object or null if invalid
 */
export function calculateEffect(row, type, metric, cc = DEFAULT_CONTINUITY_CORRECTION) {
  switch (type) {
    case 'binary':
      return calculateBinaryEffect(row, metric, cc);
    case 'continuous':
      return calculateContinuousEffect(row);
    case 'proportion':
      return calculateProportionEffect(row);
    case 'survival':
      return calculateSurvivalEffect(row);
    case 'correlation':
      return calculateCorrelationEffect(row);
    case 'generic':
      return calculateGenericEffect(row);
    default:
      return null;
  }
}

/**
 * Calculate effect size for binary outcomes
 * @param {Object} row - Study data with e1, n1, e2, n2
 * @param {string} metric - OR, RR, or RD
 * @param {number} cc - Continuity correction
 * @returns {Object|null} Effect size object
 */
function calculateBinaryEffect(row, metric, cc) {
  let e1 = parseNumber(row.e1);
  let n1 = parseNumber(row.n1);
  let e2 = parseNumber(row.e2);
  let n2 = parseNumber(row.n2);

  // Validate inputs
  if ([e1, n1, e2, n2].some(v => v === null || v < 0)) return null;
  if (n1 === 0 || n2 === 0) return null;
  if (e1 > n1 || e2 > n2) return null;

  // Apply continuity correction if needed
  const needsCC = e1 === 0 || e2 === 0 || e1 === n1 || e2 === n2;
  const e1c = needsCC ? e1 + cc : e1;
  const n1c = needsCC ? n1 + 2 * cc : n1;
  const e2c = needsCC ? e2 + cc : e2;
  const n2c = needsCC ? n2 + 2 * cc : n2;

  const p1 = e1c / n1c;
  const p2 = e2c / n2c;

  let es, vi, display;

  switch (metric) {
    case 'OR':
      es = Math.log((e1c * (n2c - e2c)) / (e2c * (n1c - e1c)));
      vi = 1 / e1c + 1 / (n1c - e1c) + 1 / e2c + 1 / (n2c - e2c);
      display = Math.exp(es);
      break;

    case 'RR':
      es = Math.log(p1 / p2);
      // Correct variance formula: (n-e)/(e*n) for each group
      vi = (n1c - e1c) / (e1c * n1c) + (n2c - e2c) / (e2c * n2c);
      display = Math.exp(es);
      break;

    case 'RD':
      es = p1 - p2;
      vi = p1 * (1 - p1) / n1c + p2 * (1 - p2) / n2c;
      display = es;
      break;

    default:
      return null;
  }

  // Validate result
  if (!Number.isFinite(es) || !Number.isFinite(vi) || vi <= 0) {
    return null;
  }

  return {
    es,
    vi,
    se: Math.sqrt(vi),
    display,
    raw: {
      e1,
      n1,
      e2,
      n2,
      p1: e1 / n1,
      p2: e2 / n2
    }
  };
}

/**
 * Calculate standardized mean difference (Hedges' g)
 * @param {Object} row - Study data with m1, s1, n1, m2, s2, n2
 * @returns {Object|null} Effect size object
 */
function calculateContinuousEffect(row) {
  const m1 = parseNumber(row.m1);
  const s1 = parseNumber(row.s1);
  const n1 = parseNumber(row.n1);
  const m2 = parseNumber(row.m2);
  const s2 = parseNumber(row.s2);
  const n2 = parseNumber(row.n2);

  // Validate inputs
  if ([m1, s1, n1, m2, s2, n2].some(v => v === null)) return null;
  if (s1 <= 0 || s2 <= 0 || n1 <= 0 || n2 <= 0) return null;

  // Calculate pooled SD
  const pooledSD = Math.sqrt(((n1 - 1) * s1 * s1 + (n2 - 1) * s2 * s2) / (n1 + n2 - 2));

  // Cohen's d
  const d = (m1 - m2) / pooledSD;

  // Hedges' g correction factor
  const j = 1 - 3 / (4 * (n1 + n2 - 2) - 1);

  // Hedges' g (bias-corrected SMD)
  const es = d * j;

  // Variance of Hedges' g
  const vi = (n1 + n2) / (n1 * n2) + (es * es) / (2 * (n1 + n2));

  if (!Number.isFinite(es) || !Number.isFinite(vi) || vi <= 0) {
    return null;
  }

  return {
    es,
    vi,
    se: Math.sqrt(vi),
    display: es,
    raw: { m1, s1, n1, m2, s2, n2 }
  };
}

/**
 * Calculate effect size for single proportion (logit transformation)
 * @param {Object} row - Study data with e, n
 * @returns {Object|null} Effect size object
 */
function calculateProportionEffect(row) {
  const e = parseNumber(row.e);
  const n = parseNumber(row.n);

  if (e === null || n === null || n <= 0 || e < 0 || e > n) return null;

  // Apply continuity correction
  const pAdj = (e + 0.5) / (n + 1);

  // Logit transformation
  const es = Math.log(pAdj / (1 - pAdj));
  const vi = 1 / (n * pAdj * (1 - pAdj));

  if (!Number.isFinite(es) || !Number.isFinite(vi) || vi <= 0) {
    return null;
  }

  return {
    es,
    vi,
    se: Math.sqrt(vi),
    display: e / n,
    raw: { e, n }
  };
}

/**
 * Calculate effect size for survival data (hazard ratio)
 * @param {Object} row - Study data with hr, ll, ul (or se)
 * @returns {Object|null} Effect size object
 */
function calculateSurvivalEffect(row) {
  const hr = parseNumber(row.hr);
  const ll = parseNumber(row.ll);
  const ul = parseNumber(row.ul);

  if (hr === null || hr <= 0) return null;

  const es = Math.log(hr);
  let vi;

  if (ll && ul && ll > 0 && ul > 0 && ll < ul) {
    // Calculate SE from CI
    const se = (Math.log(ul) - Math.log(ll)) / (2 * 1.96);
    vi = se * se;
  } else {
    // Use provided SE if available
    const se = parseNumber(row.se);
    if (!se || se <= 0) return null;
    vi = se * se;
  }

  if (!Number.isFinite(es) || !Number.isFinite(vi) || vi <= 0) {
    return null;
  }

  return {
    es,
    vi,
    se: Math.sqrt(vi),
    display: hr,
    raw: { hr, ll, ul }
  };
}

/**
 * Calculate effect size for correlation (Fisher's z transformation)
 * @param {Object} row - Study data with r, n
 * @returns {Object|null} Effect size object
 */
function calculateCorrelationEffect(row) {
  const r = parseNumber(row.r);
  const n = parseNumber(row.n);

  if (r === null || n === null || n <= 3) return null;
  if (Math.abs(r) >= 1) return null;

  // Fisher's z transformation
  const es = 0.5 * Math.log((1 + r) / (1 - r));
  const vi = 1 / (n - 3);

  if (!Number.isFinite(es) || !Number.isFinite(vi) || vi <= 0) {
    return null;
  }

  return {
    es,
    vi,
    se: Math.sqrt(vi),
    display: r,
    raw: { r, n }
  };
}

/**
 * Use pre-computed effect size and SE
 * @param {Object} row - Study data with es, se
 * @returns {Object|null} Effect size object
 */
function calculateGenericEffect(row) {
  const es = parseNumber(row.es);
  const se = parseNumber(row.se);

  if (es === null || se === null || se <= 0) return null;

  return {
    es,
    vi: se * se,
    se,
    display: es
  };
}

/**
 * Transform effect back from log scale for ratios
 * @param {number} logES - Effect size on log scale
 * @param {string} metric - Effect measure
 * @returns {number} Effect on natural scale
 */
export function transformEffect(logES, metric) {
  if (['OR', 'RR', 'HR'].includes(metric)) {
    return Math.exp(logES);
  }
  return logES;
}

/**
 * Check if metric is ratio-based
 * @param {string} metric - Effect measure
 * @returns {boolean} True if ratio-based
 */
export function isRatioMetric(metric) {
  return ['OR', 'RR', 'HR'].includes(metric);
}
