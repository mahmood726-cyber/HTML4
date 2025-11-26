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
 * Begg's rank correlation test
 * @param {Array} effects - Array of effect objects
 * @returns {Object} Test results { p, tau, z }
 */
export function beggTest(effects) {
  const active = effects.filter(e => !e.excluded);
  const n = active.length;

  if (n < 3) {
    return { p: null, tau: null, z: null };
  }

  // Sort by effect size and assign ranks
  const sorted = [...active].sort((a, b) => a.es - b.es);
  const withRanks = sorted.map((s, i) => ({ ...s, rank: i + 1 }));

  // Calculate Kendall's tau between effect size rank and SE
  let concordant = 0;
  let discordant = 0;

  for (let i = 0; i < n - 1; i++) {
    for (let j = i + 1; j < n; j++) {
      const seComp = withRanks[i].se - withRanks[j].se;
      const rankComp = withRanks[i].rank - withRanks[j].rank;

      if (seComp * rankComp > 0) {
        concordant++;
      } else if (seComp * rankComp < 0) {
        discordant++;
      }
      // Ties (seComp * rankComp === 0) are ignored
    }
  }

  const tau = (concordant - discordant) / (n * (n - 1) / 2);
  const z = 3 * tau * Math.sqrt(n * (n - 1)) / Math.sqrt(2 * (2 * n + 5));
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
 * Trim and fill analysis
 * @param {Array} effects - Array of effect objects
 * @param {string} side - Side to trim ('right' or 'left', default: 'right')
 * @returns {Object} Results { k0, imputedStudies, adjusted }
 */
export function trimAndFill(effects, side = 'right') {
  const active = effects.filter(e => !e.excluded).map(e => ({ ...e }));
  const n = active.length;

  if (n < 3) {
    return { k0: 0, imputedStudies: [], adjusted: null };
  }

  // Sort by effect size
  active.sort((a, b) => a.es - b.es);

  // Find median
  const median = active[Math.floor(n / 2)].es;

  // Count studies on each side
  const left = active.filter(e => e.es < median).length;
  const right = active.filter(e => e.es > median).length;

  // Estimate number of missing studies (simplified R0 estimator)
  const k0 = Math.abs(right - left);

  if (k0 === 0) {
    return { k0: 0, imputedStudies: [], adjusted: null };
  }

  // Identify studies to mirror
  const toMirror = side === 'right' || right > left
    ? active.slice(-k0)
    : active.slice(0, k0);

  // Create imputed studies
  const imputedStudies = toMirror.map(e => ({
    id: `Imputed (${e.id})`,
    es: 2 * median - e.es,
    vi: e.vi,
    se: e.se,
    excluded: false,
    imputed: true
  }));

  // Pool with imputed studies
  const augmented = [...active, ...imputedStudies];
  const tau2 = tauDL(augmented);
  const adjusted = poolInverseVariance(augmented, tau2, { useHKSJ: false });

  return { k0, imputedStudies, adjusted };
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
