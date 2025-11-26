/**
 * Mantel-Haenszel and Peto methods for binary outcome meta-analysis
 * Optimal for sparse data and rare events (Bradburn et al. 2007)
 */

import { pFromZ, chiSqP } from './distributions.js';

/**
 * Mantel-Haenszel meta-analysis for binary outcomes
 * Does not require continuity correction, handles zero cells naturally
 * Reference: Mantel & Haenszel (1959), Greenland & Robins (1985)
 *
 * @param {Array} studies - Array of study objects with e1, n1, e2, n2
 * @param {string} metric - Effect measure: 'OR', 'RR', or 'RD'
 * @param {Object} options - Analysis options
 * @returns {Object} MH pooled result
 */
export function mantelHaenszel(studies, metric = 'OR', options = {}) {
  const {
    confLevel = 0.95,
    includeZeroZero = false  // Include double-zero studies
  } = options;

  // Filter to valid binary studies
  const valid = studies.filter(s => {
    if (!s.raw || s.excluded) return false;
    const { e1, n1, e2, n2 } = s.raw;
    if ([e1, n1, e2, n2].some(v => v == null || v < 0)) return false;
    if (n1 === 0 || n2 === 0) return false;

    // Handle double-zero studies
    const isDoubleZero = (e1 === 0 && e2 === 0) || (e1 === n1 && e2 === n2);
    if (isDoubleZero && !includeZeroZero) return false;

    return true;
  });

  if (valid.length === 0) return null;

  const k = valid.length;
  const alpha = 1 - confLevel;
  const z = Math.abs(normalQuantile(alpha / 2));

  let result;

  switch (metric) {
    case 'OR':
      result = mhOddsRatio(valid, z);
      break;
    case 'RR':
      result = mhRiskRatio(valid, z);
      break;
    case 'RD':
      result = mhRiskDifference(valid, z);
      break;
    default:
      return null;
  }

  if (!result) return null;

  // Add heterogeneity statistics
  const { Q, I2, pHet } = calculateMHHeterogeneity(valid, result.es, metric);

  return {
    ...result,
    method: 'Mantel-Haenszel',
    metric,
    k,
    Q,
    I2,
    pHet,
    df: k - 1,
    studies: valid.map(s => ({
      id: s.id,
      es: s.es,
      se: s.se,
      raw: s.raw
    }))
  };
}

/**
 * Mantel-Haenszel Odds Ratio
 * Reference: Robins, Breslow & Greenland (1986)
 */
function mhOddsRatio(studies, z) {
  let R = 0, S = 0;
  let varR = 0, varS = 0, covRS = 0;

  studies.forEach(s => {
    const { e1, n1, e2, n2 } = s.raw;
    const a = e1, b = n1 - e1;  // Treatment: events, non-events
    const c = e2, d = n2 - e2;  // Control: events, non-events
    const N = n1 + n2;

    // MH components
    R += (a * d) / N;
    S += (b * c) / N;

    // Variance components (Robins-Breslow-Greenland)
    const P = (a + d) / N;
    const Q = (b + c) / N;

    varR += (a * d * (a + d)) / (N * N);
    varS += (b * c * (b + c)) / (N * N);
    covRS += (a * d * Q + b * c * P) / (N * N);
  });

  if (S === 0) return null;

  const OR = R / S;
  const lnOR = Math.log(OR);

  // Variance of ln(OR) using RBG formula
  const varLnOR = varR / (2 * R * R) + covRS / (2 * R * S) + varS / (2 * S * S);
  const seLnOR = Math.sqrt(varLnOR);

  // Test statistic
  const zStat = lnOR / seLnOR;
  const pVal = 2 * (1 - normalCDF(Math.abs(zStat)));

  return {
    es: lnOR,
    se: seLnOR,
    display: OR,
    ciLo: Math.exp(lnOR - z * seLnOR),
    ciHi: Math.exp(lnOR + z * seLnOR),
    z: zStat,
    pVal,
    isRatio: true
  };
}

/**
 * Mantel-Haenszel Risk Ratio
 * Reference: Greenland & Robins (1985)
 */
function mhRiskRatio(studies, z) {
  let R = 0, S = 0;
  let varLnRR = 0;

  studies.forEach(s => {
    const { e1, n1, e2, n2 } = s.raw;
    const N = n1 + n2;

    // MH components for RR
    R += (e1 * n2) / N;
    S += (e2 * n1) / N;
  });

  if (S === 0 || R === 0) return null;

  const RR = R / S;
  const lnRR = Math.log(RR);

  // Greenland-Robins variance estimator
  let P = 0;
  studies.forEach(s => {
    const { e1, n1, e2, n2 } = s.raw;
    const N = n1 + n2;
    P += ((n1 * n2 * (e1 + e2)) / (N * N) - (e1 * e2) / N);
  });

  varLnRR = P / (R * S);
  const seLnRR = Math.sqrt(varLnRR);

  const zStat = lnRR / seLnRR;
  const pVal = 2 * (1 - normalCDF(Math.abs(zStat)));

  return {
    es: lnRR,
    se: seLnRR,
    display: RR,
    ciLo: Math.exp(lnRR - z * seLnRR),
    ciHi: Math.exp(lnRR + z * seLnRR),
    z: zStat,
    pVal,
    isRatio: true
  };
}

/**
 * Mantel-Haenszel Risk Difference
 * Reference: Greenland & Robins (1985)
 */
function mhRiskDifference(studies, z) {
  let sumW = 0, sumWRD = 0;
  let sumW2V = 0;

  studies.forEach(s => {
    const { e1, n1, e2, n2 } = s.raw;
    const N = n1 + n2;

    const p1 = e1 / n1;
    const p2 = e2 / n2;
    const rd = p1 - p2;

    // MH weight for RD
    const w = (n1 * n2) / N;

    // Variance within study
    const v = (p1 * (1 - p1) / n1) + (p2 * (1 - p2) / n2);

    sumW += w;
    sumWRD += w * rd;
    sumW2V += w * w * v;
  });

  if (sumW === 0) return null;

  const RD = sumWRD / sumW;
  const varRD = sumW2V / (sumW * sumW);
  const seRD = Math.sqrt(varRD);

  const zStat = RD / seRD;
  const pVal = 2 * (1 - normalCDF(Math.abs(zStat)));

  return {
    es: RD,
    se: seRD,
    display: RD,
    ciLo: RD - z * seRD,
    ciHi: RD + z * seRD,
    z: zStat,
    pVal,
    isRatio: false
  };
}

/**
 * Peto Odds Ratio method
 * Optimal for rare events with balanced groups
 * Reference: Yusuf et al. (1985), Peto (1987)
 *
 * @param {Array} studies - Array of study objects with e1, n1, e2, n2
 * @param {Object} options - Analysis options
 * @returns {Object} Peto pooled result
 */
export function petoOddsRatio(studies, options = {}) {
  const {
    confLevel = 0.95,
    excludeZeroZero = true
  } = options;

  // Filter valid studies
  const valid = studies.filter(s => {
    if (!s.raw || s.excluded) return false;
    const { e1, n1, e2, n2 } = s.raw;
    if ([e1, n1, e2, n2].some(v => v == null || v < 0)) return false;
    if (n1 === 0 || n2 === 0) return false;

    // Exclude double-zero studies for Peto
    if (excludeZeroZero && e1 === 0 && e2 === 0) return false;
    if (excludeZeroZero && e1 === n1 && e2 === n2) return false;

    return true;
  });

  if (valid.length === 0) return null;

  const k = valid.length;
  const alpha = 1 - confLevel;
  const z = Math.abs(normalQuantile(alpha / 2));

  let O_E = 0;  // Sum of (Observed - Expected)
  let V = 0;    // Sum of hypergeometric variances

  valid.forEach(s => {
    const { e1, n1, e2, n2 } = s.raw;
    const N = n1 + n2;
    const E = e1 + e2;  // Total events

    // Expected events in treatment group under null
    const expected = (n1 * E) / N;

    // Hypergeometric variance
    const variance = (n1 * n2 * E * (N - E)) / (N * N * (N - 1));

    O_E += e1 - expected;
    V += variance;
  });

  if (V === 0) return null;

  // Peto log odds ratio
  const lnOR = O_E / V;
  const seLnOR = 1 / Math.sqrt(V);

  const OR = Math.exp(lnOR);
  const zStat = lnOR / seLnOR;
  const pVal = 2 * (1 - normalCDF(Math.abs(zStat)));

  // Heterogeneity test (Peto method)
  let Q = 0;
  valid.forEach(s => {
    const { e1, n1, e2, n2 } = s.raw;
    const N = n1 + n2;
    const E = e1 + e2;
    const expected = (n1 * E) / N;
    const variance = (n1 * n2 * E * (N - E)) / (N * N * (N - 1));

    if (variance > 0) {
      const contrib = ((e1 - expected) - lnOR * variance) ** 2 / variance;
      Q += contrib;
    }
  });

  const df = k - 1;
  const I2 = df > 0 ? Math.max(0, (Q - df) / Q * 100) : 0;
  const pHet = chiSqP(Q, df);

  return {
    method: 'Peto',
    metric: 'OR',
    es: lnOR,
    se: seLnOR,
    display: OR,
    ciLo: Math.exp(lnOR - z * seLnOR),
    ciHi: Math.exp(lnOR + z * seLnOR),
    z: zStat,
    pVal,
    isRatio: true,
    k,
    Q,
    I2,
    pHet,
    df,
    O_E,
    V,
    studies: valid.map(s => ({
      id: s.id,
      raw: s.raw
    }))
  };
}

/**
 * Calculate heterogeneity statistics for MH analysis
 */
function calculateMHHeterogeneity(studies, pooledES, metric) {
  const k = studies.length;
  const df = k - 1;

  if (df <= 0) {
    return { Q: 0, I2: 0, pHet: 1 };
  }

  // Cochran's Q for MH
  let Q = 0;

  studies.forEach(s => {
    const { e1, n1, e2, n2 } = s.raw;
    let studyES, studyVar;

    if (metric === 'OR') {
      const a = e1, b = n1 - e1, c = e2, d = n2 - e2;
      if (a > 0 && b > 0 && c > 0 && d > 0) {
        studyES = Math.log((a * d) / (b * c));
        studyVar = 1/a + 1/b + 1/c + 1/d;
      } else return;
    } else if (metric === 'RR') {
      if (e1 > 0 && e2 > 0) {
        studyES = Math.log((e1/n1) / (e2/n2));
        studyVar = (n1 - e1)/(e1 * n1) + (n2 - e2)/(e2 * n2);
      } else return;
    } else {
      const p1 = e1/n1, p2 = e2/n2;
      studyES = p1 - p2;
      studyVar = (p1*(1-p1)/n1) + (p2*(1-p2)/n2);
    }

    if (studyVar > 0) {
      Q += (studyES - pooledES) ** 2 / studyVar;
    }
  });

  const I2 = Math.max(0, (Q - df) / Q * 100);
  const pHet = chiSqP(Q, df);

  return { Q, I2, pHet };
}

/**
 * Normal CDF approximation
 */
function normalCDF(x) {
  const a1 = 0.254829592;
  const a2 = -0.284496736;
  const a3 = 1.421413741;
  const a4 = -1.453152027;
  const a5 = 1.061405429;
  const p = 0.3275911;

  const sign = x < 0 ? -1 : 1;
  x = Math.abs(x) / Math.sqrt(2);

  const t = 1.0 / (1.0 + p * x);
  const y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-x * x);

  return 0.5 * (1.0 + sign * y);
}

/**
 * Normal quantile (inverse CDF) approximation
 */
function normalQuantile(p) {
  if (p <= 0) return -Infinity;
  if (p >= 1) return Infinity;
  if (p === 0.5) return 0;

  const a = [
    -3.969683028665376e+01, 2.209460984245205e+02,
    -2.759285104469687e+02, 1.383577518672690e+02,
    -3.066479806614716e+01, 2.506628277459239e+00
  ];
  const b = [
    -5.447609879822406e+01, 1.615858368580409e+02,
    -1.556989798598866e+02, 6.680131188771972e+01,
    -1.328068155288572e+01
  ];
  const c = [
    -7.784894002430293e-03, -3.223964580411365e-01,
    -2.400758277161838e+00, -2.549732539343734e+00,
    4.374664141464968e+00, 2.938163982698783e+00
  ];
  const d = [
    7.784695709041462e-03, 3.224671290700398e-01,
    2.445134137142996e+00, 3.754408661907416e+00
  ];

  const pLow = 0.02425;
  const pHigh = 1 - pLow;
  let q, r;

  if (p < pLow) {
    q = Math.sqrt(-2 * Math.log(p));
    return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
           ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
  } else if (p <= pHigh) {
    q = p - 0.5;
    r = q * q;
    return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
           (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
  } else {
    q = Math.sqrt(-2 * Math.log(1 - p));
    return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
            ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
  }
}

/**
 * Select optimal method for binary meta-analysis
 * Based on Bradburn et al. (2007) recommendations
 *
 * @param {Array} studies - Study data
 * @param {string} metric - Effect measure
 * @returns {Object} Recommendation with method and rationale
 */
export function recommendBinaryMethod(studies) {
  const valid = studies.filter(s => s.raw && !s.excluded);

  if (valid.length === 0) {
    return { method: 'IV', rationale: 'No valid studies' };
  }

  // Calculate characteristics
  let totalEvents = 0;
  let totalN = 0;
  let zeroStudies = 0;
  let doubleZero = 0;
  let unbalanced = 0;

  valid.forEach(s => {
    const { e1, n1, e2, n2 } = s.raw;
    totalEvents += e1 + e2;
    totalN += n1 + n2;

    if (e1 === 0 || e2 === 0 || e1 === n1 || e2 === n2) zeroStudies++;
    if ((e1 === 0 && e2 === 0) || (e1 === n1 && e2 === n2)) doubleZero++;
    if (Math.abs(n1 - n2) / Math.max(n1, n2) > 0.3) unbalanced++;
  });

  const eventRate = totalEvents / totalN;
  const propZero = zeroStudies / valid.length;
  const propUnbalanced = unbalanced / valid.length;

  // Decision logic based on Bradburn et al.
  if (eventRate < 0.01) {
    // Very rare events
    if (propUnbalanced < 0.3) {
      return {
        method: 'Peto',
        rationale: 'Rare events with balanced groups - Peto method optimal'
      };
    } else {
      return {
        method: 'MH',
        rationale: 'Rare events with unbalanced groups - MH method preferred'
      };
    }
  } else if (eventRate < 0.1 || propZero > 0.2) {
    // Sparse data
    return {
      method: 'MH',
      rationale: 'Sparse data with zero cells - MH method avoids continuity correction bias'
    };
  } else {
    return {
      method: 'IV',
      rationale: 'Adequate events - Inverse variance method appropriate'
    };
  }
}
