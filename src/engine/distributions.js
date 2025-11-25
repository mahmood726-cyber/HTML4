/**
 * Statistical distribution functions
 */

import { NORM_CDF_COEFFICIENTS, T_CRIT_TABLE } from '../constants/index.js';

/**
 * Standard normal cumulative distribution function (Abramowitz & Stegun approximation)
 * @param {number} z - Z-score
 * @returns {number} Cumulative probability
 */
export function normCDF(z) {
  // Handle edge cases
  if (!Number.isFinite(z)) {
    return z < 0 ? 0 : 1;
  }

  const { a1, a2, a3, a4, a5, p } = NORM_CDF_COEFFICIENTS;
  const sign = z < 0 ? -1 : 1;
  const absZ = Math.abs(z) / Math.SQRT2;
  const t = 1 / (1 + p * absZ);
  const y = 1 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-absZ * absZ);

  return 0.5 * (1 + sign * y);
}

/**
 * Two-tailed p-value from z-score
 * @param {number} z - Z-score
 * @returns {number} Two-tailed p-value
 */
export function pFromZ(z) {
  if (!Number.isFinite(z)) return 1;
  return 2 * (1 - normCDF(Math.abs(z)));
}

/**
 * Get t-distribution critical value
 * @param {number} df - Degrees of freedom
 * @param {number} alpha - Significance level (default: 0.05)
 * @returns {number} Critical value
 */
export function tCritical(df, alpha = 0.05) {
  if (!Number.isFinite(df) || df <= 0) return 1.96;

  // Use lookup table for common values
  if (T_CRIT_TABLE[df] !== undefined) {
    return T_CRIT_TABLE[df];
  }

  // Find closest value in table or use approximation
  const knownDf = Object.keys(T_CRIT_TABLE)
    .map(Number)
    .filter(d => Number.isFinite(d))
    .sort((a, b) => a - b);

  // If df is large enough, use normal approximation
  if (df > 100) {
    return T_CRIT_TABLE[Infinity];
  }

  // Cornish-Fisher expansion for intermediate values
  const z = 1.96; // For alpha = 0.05
  return z + (z ** 3 + z) / (4 * df) + (5 * z ** 5 + 16 * z ** 3 + 3 * z) / (96 * df ** 2);
}

/**
 * T-distribution CDF approximation
 * @param {number} t - T-statistic
 * @param {number} df - Degrees of freedom
 * @returns {number} Cumulative probability
 */
export function tCDF(t, df) {
  if (!Number.isFinite(t) || !Number.isFinite(df) || df <= 0) {
    return 0.5;
  }

  // For large df, use normal approximation
  if (df > 100) {
    return normCDF(t);
  }

  // Approximation using adjusted z-score
  const z = t * Math.sqrt(1 - 1 / (4 * df) - 7 / (120 * df * df));
  return normCDF(z);
}

/**
 * Two-tailed p-value from t-statistic
 * @param {number} t - T-statistic
 * @param {number} df - Degrees of freedom
 * @returns {number} Two-tailed p-value
 */
export function pFromT(t, df) {
  if (!Number.isFinite(t) || !Number.isFinite(df) || df <= 0) {
    return 1;
  }
  return 2 * (1 - tCDF(Math.abs(t), df));
}

/**
 * Chi-squared p-value (Wilson-Hilferty approximation)
 * @param {number} x - Chi-squared statistic
 * @param {number} df - Degrees of freedom
 * @returns {number} P-value
 */
export function chiSqP(x, df) {
  if (!Number.isFinite(x) || !Number.isFinite(df) || x <= 0 || df <= 0) {
    return 1;
  }

  const z = Math.pow(x / df, 1 / 3) - (1 - 2 / (9 * df));
  const se = Math.sqrt(2 / (9 * df));

  return 1 - normCDF(z / se);
}

/**
 * Inverse normal CDF (quantile function) - Rational approximation
 * @param {number} p - Probability (0 < p < 1)
 * @returns {number} Z-score
 */
export function normQuantile(p) {
  if (p <= 0) return -Infinity;
  if (p >= 1) return Infinity;
  if (p === 0.5) return 0;

  // Rational approximation coefficients
  const a = [
    -3.969683028665376e1,
    2.209460984245205e2,
    -2.759285104469687e2,
    1.383577518672690e2,
    -3.066479806614716e1,
    2.506628277459239e0
  ];
  const b = [
    -5.447609879822406e1,
    1.615858368580409e2,
    -1.556989798598866e2,
    6.680131188771972e1,
    -1.328068155288572e1
  ];
  const c = [
    -7.784894002430293e-3,
    -3.223964580411365e-1,
    -2.400758277161838e0,
    -2.549732539343734e0,
    4.374664141464968e0,
    2.938163982698783e0
  ];
  const d = [
    7.784695709041462e-3,
    3.224671290700398e-1,
    2.445134137142996e0,
    3.754408661907416e0
  ];

  const pLow = 0.02425;
  const pHigh = 1 - pLow;

  let q, r;

  if (p < pLow) {
    q = Math.sqrt(-2 * Math.log(p));
    return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
           ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
  } else if (p <= pHigh) {
    q = p - 0.5;
    r = q * q;
    return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q /
           (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1);
  } else {
    q = Math.sqrt(-2 * Math.log(1 - p));
    return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) /
            ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1);
  }
}
