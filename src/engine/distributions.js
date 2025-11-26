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
 * T-distribution CDF using regularized incomplete beta function
 * Accurate for all df including small values
 * @param {number} t - T-statistic
 * @param {number} df - Degrees of freedom
 * @returns {number} Cumulative probability
 */
export function tCDF(t, df) {
  if (!Number.isFinite(t) || !Number.isFinite(df) || df <= 0) {
    return 0.5;
  }

  // For very large df, use normal approximation
  if (df > 1000) {
    return normCDF(t);
  }

  // t-CDF via incomplete beta function:
  // F(t) = 1 - 0.5 * I_x(df/2, 1/2) where x = df/(df + t²)
  // For t > 0: F(t) = 1 - 0.5 * I_x(df/2, 1/2)
  // For t < 0: F(t) = 0.5 * I_x(df/2, 1/2)
  const x = df / (df + t * t);
  const beta = incompleteBeta(x, df / 2, 0.5);

  if (t >= 0) {
    return 1 - 0.5 * beta;
  } else {
    return 0.5 * beta;
  }
}

/**
 * Regularized incomplete beta function I_x(a, b)
 * Uses continued fraction expansion for accuracy
 * @param {number} x - Upper limit (0 <= x <= 1)
 * @param {number} a - Shape parameter a > 0
 * @param {number} b - Shape parameter b > 0
 * @returns {number} I_x(a, b)
 */
function incompleteBeta(x, a, b) {
  if (x <= 0) return 0;
  if (x >= 1) return 1;

  // Use symmetry relation if x > (a+1)/(a+b+2)
  if (x > (a + 1) / (a + b + 2)) {
    return 1 - incompleteBeta(1 - x, b, a);
  }

  // Compute log of beta function B(a,b) = Γ(a)Γ(b)/Γ(a+b)
  const lbeta = logGamma(a) + logGamma(b) - logGamma(a + b);

  // Front factor: x^a * (1-x)^b / (a * B(a,b))
  const front = Math.exp(a * Math.log(x) + b * Math.log(1 - x) - lbeta) / a;

  // Continued fraction (Lentz's algorithm)
  const maxIter = 200;
  const eps = 1e-14;

  let f = 1;
  let C = 1;
  let D = 0;

  for (let m = 0; m <= maxIter; m++) {
    let num;
    if (m === 0) {
      num = 1;
    } else if (m % 2 === 1) {
      // Odd terms
      const k = (m - 1) / 2;
      num = -(a + k) * (a + b + k) * x / ((a + 2 * k) * (a + 2 * k + 1));
    } else {
      // Even terms
      const k = m / 2;
      num = k * (b - k) * x / ((a + 2 * k - 1) * (a + 2 * k));
    }

    D = 1 + num * D;
    if (Math.abs(D) < 1e-30) D = 1e-30;
    D = 1 / D;

    C = 1 + num / C;
    if (Math.abs(C) < 1e-30) C = 1e-30;

    const delta = C * D;
    f *= delta;

    if (Math.abs(delta - 1) < eps) {
      break;
    }
  }

  return front * (f - 1);
}

/**
 * Log-gamma function using Lanczos approximation
 * @param {number} x - Input value
 * @returns {number} ln(Γ(x))
 */
function logGamma(x) {
  if (x <= 0) return Infinity;

  // Lanczos coefficients for g=7
  const c = [
    0.99999999999980993,
    676.5203681218851,
    -1259.1392167224028,
    771.32342877765313,
    -176.61502916214059,
    12.507343278686905,
    -0.13857109526572012,
    9.9843695780195716e-6,
    1.5056327351493116e-7
  ];

  if (x < 0.5) {
    // Reflection formula
    return Math.log(Math.PI / Math.sin(Math.PI * x)) - logGamma(1 - x);
  }

  x -= 1;
  let sum = c[0];
  for (let i = 1; i < c.length; i++) {
    sum += c[i] / (x + i);
  }

  const t = x + 7.5;
  return 0.5 * Math.log(2 * Math.PI) + (x + 0.5) * Math.log(t) - t + Math.log(sum);
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
