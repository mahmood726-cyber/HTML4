/**
 * Publication bias tests and adjustments
 */

import { pFromT, pFromZ } from './distributions.js';
import { tauDL } from './tau.js';
import { poolInverseVariance } from './pooling.js';

/**
 * Egger's regression test for funnel plot asymmetry
 * Uses weighted least squares as recommended (Sterne & Egger 2005)
 * @param {Array} effects - Array of effect objects
 * @returns {Object} Test results { p, intercept, t, slope }
 */
export function eggerTest(effects) {
  const active = effects.filter(e => !e.excluded);
  const n = active.length;

  if (n < 3) {
    return { p: null, intercept: null, t: null, slope: null };
  }

  // Standardized effect (y) vs precision (x), weighted by inverse variance
  const x = active.map(e => 1 / e.se);
  const y = active.map(e => e.es / e.se);
  const w = active.map(e => 1 / e.vi); // Inverse variance weights

  // Weighted linear regression
  const sumW = w.reduce((a, b) => a + b, 0);
  const sumWX = w.reduce((acc, wi, i) => acc + wi * x[i], 0);
  const sumWY = w.reduce((acc, wi, i) => acc + wi * y[i], 0);
  const sumWXY = w.reduce((acc, wi, i) => acc + wi * x[i] * y[i], 0);
  const sumWX2 = w.reduce((acc, wi, i) => acc + wi * x[i] * x[i], 0);

  const denominator = sumW * sumWX2 - sumWX * sumWX;
  if (Math.abs(denominator) < 1e-10) {
    return { p: null, intercept: null, t: null, slope: null };
  }

  const slope = (sumW * sumWXY - sumWX * sumWY) / denominator;
  const intercept = (sumWY - slope * sumWX) / sumW;

  // Calculate weighted residual sum of squares
  const yPred = x.map(xi => intercept + slope * xi);
  const wsse = w.reduce((acc, wi, i) => acc + wi * (y[i] - yPred[i]) ** 2, 0);
  const mse = wsse / (n - 2);

  // Standard error of intercept
  const seInt = Math.sqrt(mse * sumWX2 / denominator);

  // T-test for intercept (tests for asymmetry)
  const t = intercept / seInt;
  const p = pFromT(t, n - 2);

  return { p, intercept, t, slope };
}

/**
 * Begg's rank correlation test (adjusted for ties)
 * Uses Kendall's tau-b with tie corrections
 * @param {Array} effects - Array of effect objects
 * @returns {Object} Test results { p, tau, z }
 */
export function beggTest(effects) {
  const active = effects.filter(e => !e.excluded);
  const n = active.length;

  if (n < 3) {
    return { p: null, tau: null, z: null };
  }

  // Get standardized effect sizes (adjusted for variance)
  // Begg & Mazumdar recommend using variance-standardized residuals
  const tau2 = tauDL(active);
  const w = active.map(e => 1 / (e.vi + tau2));
  const sumW = w.reduce((a, b) => a + b, 0);
  const mu = active.reduce((acc, e, i) => acc + w[i] * e.es, 0) / sumW;

  // Standardized effect: (es - mu) / se
  const data = active.map(e => ({
    standardizedES: (e.es - mu) / e.se,
    variance: e.vi
  }));

  // Calculate Kendall's tau-b between standardized ES and variance
  let concordant = 0;
  let discordant = 0;
  let tiesX = 0; // Ties in standardized ES
  let tiesY = 0; // Ties in variance
  let tiesXY = 0; // Joint ties

  for (let i = 0; i < n - 1; i++) {
    for (let j = i + 1; j < n; j++) {
      const xDiff = data[i].standardizedES - data[j].standardizedES;
      const yDiff = data[i].variance - data[j].variance;

      const xTied = Math.abs(xDiff) < 1e-10;
      const yTied = Math.abs(yDiff) < 1e-10;

      if (xTied && yTied) {
        tiesXY++;
      } else if (xTied) {
        tiesX++;
      } else if (yTied) {
        tiesY++;
      } else if (xDiff * yDiff > 0) {
        concordant++;
      } else {
        discordant++;
      }
    }
  }

  const nPairs = n * (n - 1) / 2;

  // Tau-b formula with tie correction
  const denominator = Math.sqrt((nPairs - tiesX - tiesXY) * (nPairs - tiesY - tiesXY));

  if (denominator < 1e-10) {
    return { p: 1, tau: 0, z: 0 };
  }

  const tau = (concordant - discordant) / denominator;

  // Variance of tau-b with tie correction (Kendall 1970)
  const v0 = n * (n - 1) * (2 * n + 5);
  const vt = tiesX * (tiesX - 1) * (2 * tiesX + 5);
  const vu = tiesY * (tiesY - 1) * (2 * tiesY + 5);
  const v1 = (tiesX * (tiesX - 1)) * (tiesY * (tiesY - 1)) / (2 * n * (n - 1));
  const v2 = (tiesX * (tiesX - 1) * (tiesX - 2)) * (tiesY * (tiesY - 1) * (tiesY - 2)) /
             (9 * n * (n - 1) * (n - 2));

  const variance = (v0 - vt - vu) / 18 + v1 + v2;
  const se = Math.sqrt(Math.max(variance, 1e-10)) / denominator;

  const z = tau / se;
  const p = pFromZ(z);

  return { p, tau, z };
}

/**
 * Peters' test for binary outcomes (recommended over Egger for OR)
 * @param {Array} effects - Array of effect objects with raw binary data
 * @returns {Object} Test results { p, intercept, t }
 */
export function petersTest(effects) {
  const active = effects.filter(e => !e.excluded && e.raw?.n1 && e.raw?.n2);
  const n = active.length;

  if (n < 3) {
    return { p: null, intercept: null, t: null };
  }

  // Peters' regression: ES against 1/sqrt(total_n)
  const x = active.map(e => 1 / Math.sqrt(e.raw.n1 + e.raw.n2));
  const y = active.map(e => e.es);
  const w = active.map(e => 1 / e.vi);

  // Weighted least squares
  const sumW = w.reduce((a, b) => a + b, 0);
  const sumWX = w.reduce((acc, wi, i) => acc + wi * x[i], 0);
  const sumWY = w.reduce((acc, wi, i) => acc + wi * y[i], 0);
  const sumWXY = w.reduce((acc, wi, i) => acc + wi * x[i] * y[i], 0);
  const sumWX2 = w.reduce((acc, wi, i) => acc + wi * x[i] * x[i], 0);

  const denominator = sumW * sumWX2 - sumWX * sumWX;
  if (denominator === 0) {
    return { p: null, intercept: null, t: null };
  }

  const slope = (sumW * sumWXY - sumWX * sumWY) / denominator;
  const intercept = (sumWY - slope * sumWX) / sumW;

  // Calculate SE and t-test
  const yPred = x.map(xi => intercept + slope * xi);
  const sse = w.reduce((acc, wi, i) => acc + wi * (y[i] - yPred[i]) ** 2, 0);
  const mse = sse / (n - 2);
  const seSlope = Math.sqrt(mse * sumW / denominator);
  const t = slope / seSlope;
  const p = pFromT(t, n - 2);

  return { p, intercept, t };
}

/**
 * Rosenthal's fail-safe N
 * @param {Array} effects - Array of effect objects
 * @param {number} targetAlpha - Target significance level (default: 0.05)
 * @returns {number} Number of null studies needed
 */
export function failSafeN(effects, targetAlpha = 0.05) {
  const active = effects.filter(e => !e.excluded);

  if (active.length < 2) return 0;

  // Critical z for target alpha (one-tailed)
  const zCrit = 1.645; // For alpha = 0.05

  // Sum of z-scores
  const zSum = active.reduce((acc, e) => acc + e.es / e.se, 0);

  // Rosenthal formula
  const fsn = Math.max(0, Math.floor((zSum * zSum) / (zCrit * zCrit) - active.length));

  return fsn;
}

/**
 * Orwin's fail-safe N (for standardized mean differences)
 * @param {Array} effects - Array of effect objects
 * @param {number} criterion - Criterion effect size (default: 0.1)
 * @param {number} nullEffect - Assumed effect of null studies (default: 0)
 * @returns {number} Number of null studies needed
 */
export function orwinFailSafeN(effects, criterion = 0.1, nullEffect = 0) {
  const active = effects.filter(e => !e.excluded);

  if (active.length < 2) return 0;

  const meanES = active.reduce((acc, e) => acc + e.es, 0) / active.length;

  if (Math.abs(meanES - criterion) < 0.0001) return Infinity;
  if (Math.abs(criterion - nullEffect) < 0.0001) return Infinity;

  const fsn = Math.ceil(active.length * (meanES - criterion) / (criterion - nullEffect));

  return Math.max(0, fsn);
}

/**
 * Trim and fill analysis using R₀ estimator (Duval & Tweedie 2000)
 * @param {Array} effects - Array of effect objects
 * @param {string} estimator - Estimator type: 'R0' or 'L0' (default: 'R0')
 * @param {string} side - Side to impute: 'right', 'left', or 'auto' (default: 'auto')
 * @param {number} maxIter - Maximum iterations (default: 10)
 * @returns {Object} Results { k0, imputedStudies, adjusted, side }
 */
export function trimAndFill(effects, estimator = 'R0', side = 'auto', maxIter = 10) {
  const active = effects.filter(e => !e.excluded).map(e => ({ ...e }));
  const n = active.length;

  if (n < 3) {
    return { k0: 0, imputedStudies: [], adjusted: null, side: null };
  }

  // Initial pooled estimate to determine center
  let tau2 = tauDL(active);
  let pooled = poolInverseVariance(active, tau2, { useHKSJ: false });
  let center = pooled.es;

  // Determine side if auto
  if (side === 'auto') {
    // Check which side has more extreme studies
    const deviations = active.map(e => e.es - center);
    const rightSkew = deviations.filter(d => d > 0).reduce((a, b) => a + b, 0);
    const leftSkew = Math.abs(deviations.filter(d => d < 0).reduce((a, b) => a + b, 0));
    side = rightSkew > leftSkew ? 'right' : 'left';
  }

  // Iterative trim and fill
  let k0 = 0;
  let trimmed = [...active];

  for (let iter = 0; iter < maxIter; iter++) {
    // Sort by distance from center
    const sorted = [...trimmed].sort((a, b) => {
      const distA = side === 'right' ? a.es - center : center - a.es;
      const distB = side === 'right' ? b.es - center : center - b.es;
      return distB - distA; // Most extreme first
    });

    // Calculate ranks (1 = most extreme on asymmetric side)
    const ranks = sorted.map((_, i) => i + 1);

    // R₀ estimator: k0 = max(0, round((4*T - n) / (2*n + 1)))
    // where T = sum of ranks for studies on asymmetric side
    let T = 0;
    sorted.forEach((e, i) => {
      const deviation = side === 'right' ? e.es - center : center - e.es;
      if (deviation > 0) {
        T += ranks[i];
      }
    });

    let k0New;
    if (estimator === 'L0') {
      // L₀ estimator: k0 = round((4*S² - n) / (2*n - 1))
      // where S is count of positive deviations
      const S = sorted.filter(e => (side === 'right' ? e.es - center : center - e.es) > 0).length;
      k0New = Math.max(0, Math.round((4 * S * S - n) / (2 * n - 1)));
    } else {
      // R₀ estimator
      k0New = Math.max(0, Math.round((4 * T - n * (n + 1) / 2) / (2 * n - 1)));
    }

    // Check convergence
    if (k0New === k0) {
      break;
    }
    k0 = k0New;

    // Trim k0 most extreme studies and recalculate center
    if (k0 > 0 && k0 < n) {
      trimmed = sorted.slice(k0);
      tau2 = tauDL(trimmed);
      pooled = poolInverseVariance(trimmed, tau2, { useHKSJ: false });
      if (pooled) {
        center = pooled.es;
      }
    }
  }

  if (k0 === 0) {
    return { k0: 0, imputedStudies: [], adjusted: pooled, side };
  }

  // Identify studies to mirror (k0 most extreme on asymmetric side)
  const sortedByExtreme = [...active].sort((a, b) => {
    const distA = side === 'right' ? a.es - center : center - a.es;
    const distB = side === 'right' ? b.es - center : center - b.es;
    return distB - distA;
  });

  const toMirror = sortedByExtreme.slice(0, Math.min(k0, n - 1));

  // Create imputed studies (mirror around center)
  const imputedStudies = toMirror.map(e => ({
    id: `Imputed (${e.id})`,
    es: 2 * center - e.es,
    vi: e.vi,
    se: e.se,
    excluded: false,
    imputed: true
  }));

  // Pool with imputed studies
  const augmented = [...active, ...imputedStudies];
  tau2 = tauDL(augmented);
  const adjusted = poolInverseVariance(augmented, tau2, { useHKSJ: false });

  return { k0, imputedStudies, adjusted, side };
}

/**
 * Calculate funnel plot data for visualization
 * @param {Array} effects - Array of effect objects
 * @param {Object} pooledResult - Pooled meta-analysis result
 * @param {boolean} isRatio - Whether effect is ratio-based
 * @returns {Object} Funnel plot data
 */
export function getFunnelPlotData(effects, pooledResult, isRatio = false) {
  const active = effects.filter(e => !e.excluded);

  if (active.length < 3) return null;

  const pooledES = pooledResult.es;
  const maxSE = Math.max(...active.map(e => e.se));

  // Generate pseudo-confidence interval lines
  const seRange = Array.from({ length: 50 }, (_, i) => (i * maxSE / 49) + 0.01);

  const transform = v => isRatio ? Math.exp(v) : v;

  return {
    studies: active.map(e => ({
      id: e.id,
      x: transform(e.es),
      y: e.se
    })),
    pooledES: transform(pooledES),
    ciLines: {
      lower: seRange.map(se => transform(pooledES - 1.96 * se)),
      upper: seRange.map(se => transform(pooledES + 1.96 * se)),
      se: seRange
    },
    maxSE
  };
}

/**
 * Calculate contour-enhanced funnel plot regions
 * @param {Object} pooledResult - Pooled result
 * @param {number} maxSE - Maximum SE for plot bounds
 * @returns {Object} Contour data for significance regions
 */
export function getContourFunnelData(pooledResult, maxSE) {
  const seRange = Array.from({ length: 50 }, (_, i) => (i * maxSE / 49) + 0.01);

  // Significance contours at p = 0.01, 0.05, 0.10
  const contours = [
    { p: 0.01, z: 2.576, color: 'rgba(220, 38, 38, 0.1)' },
    { p: 0.05, z: 1.96, color: 'rgba(234, 179, 8, 0.1)' },
    { p: 0.10, z: 1.645, color: 'rgba(34, 197, 94, 0.1)' }
  ];

  return contours.map(({ p, z, color }) => ({
    p,
    color,
    lower: seRange.map(se => -z * se),
    upper: seRange.map(se => z * se),
    se: seRange
  }));
}
