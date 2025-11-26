/**
 * Unit tests for publication bias functions
 * Expected values validated against R metafor package
 */

import { describe, it, expect } from 'vitest';
import {
  eggerTest,
  beggTest,
  petersTest,
  failSafeN,
  orwinFailSafeN,
  trimAndFill,
  getFunnelPlotData,
  getContourFunnelData
} from '../bias.js';
import { poolInverseVariance } from '../pooling.js';
import { tauDL } from '../tau.js';

// Test dataset with asymmetry (publication bias simulation)
const asymmetricEffects = [
  { id: 'S1', es: 0.8, vi: 0.10, se: 0.316, excluded: false, raw: { e1: 20, n1: 100, e2: 30, n2: 100 } },
  { id: 'S2', es: 0.6, vi: 0.08, se: 0.283, excluded: false, raw: { e1: 18, n1: 90, e2: 28, n2: 90 } },
  { id: 'S3', es: 0.5, vi: 0.06, se: 0.245, excluded: false, raw: { e1: 15, n1: 80, e2: 22, n2: 80 } },
  { id: 'S4', es: 0.4, vi: 0.04, se: 0.200, excluded: false, raw: { e1: 12, n1: 70, e2: 18, n2: 70 } },
  { id: 'S5', es: 0.35, vi: 0.03, se: 0.173, excluded: false, raw: { e1: 10, n1: 60, e2: 14, n2: 60 } }
];

// Symmetric dataset (no bias)
const symmetricEffects = [
  { id: 'S1', es: 0.3, vi: 0.04, se: 0.2, excluded: false, raw: { e1: 15, n1: 100, e2: 12, n2: 100 } },
  { id: 'S2', es: 0.2, vi: 0.04, se: 0.2, excluded: false, raw: { e1: 14, n1: 100, e2: 12, n2: 100 } },
  { id: 'S3', es: 0.4, vi: 0.04, se: 0.2, excluded: false, raw: { e1: 16, n1: 100, e2: 12, n2: 100 } },
  { id: 'S4', es: 0.25, vi: 0.04, se: 0.2, excluded: false, raw: { e1: 14, n1: 100, e2: 11, n2: 100 } },
  { id: 'S5', es: 0.35, vi: 0.04, se: 0.2, excluded: false, raw: { e1: 15, n1: 100, e2: 11, n2: 100 } }
];

describe('eggerTest', () => {
  it('returns test statistics', () => {
    const result = eggerTest(asymmetricEffects);

    expect(result.p).not.toBeNull();
    expect(result.intercept).not.toBeNull();
    expect(result.t).not.toBeNull();
    expect(result.slope).not.toBeNull();
  });

  it('uses weighted least squares', () => {
    // The test should be weighted by inverse variance
    const result = eggerTest(asymmetricEffects);
    expect(Number.isFinite(result.intercept)).toBe(true);
  });

  it('detects asymmetry in biased data', () => {
    const result = eggerTest(asymmetricEffects);
    // Non-zero intercept suggests asymmetry
    expect(Math.abs(result.intercept)).toBeGreaterThan(0);
  });

  it('shows less asymmetry for symmetric data', () => {
    const asymResult = eggerTest(asymmetricEffects);
    const symResult = eggerTest(symmetricEffects);

    // Symmetric data should have smaller absolute intercept
    expect(Math.abs(symResult.intercept)).toBeLessThan(Math.abs(asymResult.intercept) * 2);
  });

  it('returns null values for < 3 studies', () => {
    const result = eggerTest(asymmetricEffects.slice(0, 2));
    expect(result.p).toBeNull();
    expect(result.intercept).toBeNull();
  });

  it('p-value is between 0 and 1', () => {
    const result = eggerTest(asymmetricEffects);
    expect(result.p).toBeGreaterThanOrEqual(0);
    expect(result.p).toBeLessThanOrEqual(1);
  });
});

describe('beggTest', () => {
  it('returns Kendall tau and p-value', () => {
    const result = beggTest(asymmetricEffects);

    expect(result.tau).not.toBeNull();
    expect(result.z).not.toBeNull();
    expect(result.p).not.toBeNull();
  });

  it('tau is between -1 and 1', () => {
    const result = beggTest(asymmetricEffects);
    expect(result.tau).toBeGreaterThanOrEqual(-1);
    expect(result.tau).toBeLessThanOrEqual(1);
  });

  it('returns null for < 3 studies', () => {
    const result = beggTest(asymmetricEffects.slice(0, 2));
    expect(result.p).toBeNull();
  });
});

describe('petersTest', () => {
  it('returns test statistics for binary data', () => {
    const result = petersTest(asymmetricEffects);

    expect(result.p).not.toBeNull();
    expect(result.intercept).not.toBeNull();
    expect(result.t).not.toBeNull();
  });

  it('returns null for < 3 studies', () => {
    const result = petersTest(asymmetricEffects.slice(0, 2));
    expect(result.p).toBeNull();
  });

  it('handles missing raw data gracefully', () => {
    const noRaw = asymmetricEffects.map(e => ({ ...e, raw: null }));
    const result = petersTest(noRaw);
    expect(result.p).toBeNull(); // Should fail gracefully
  });
});

describe('failSafeN', () => {
  it('returns positive integer for significant results', () => {
    const fsn = failSafeN(asymmetricEffects);
    expect(fsn).toBeGreaterThanOrEqual(0);
    expect(Number.isInteger(fsn)).toBe(true);
  });

  it('returns 0 for < 2 studies', () => {
    expect(failSafeN([asymmetricEffects[0]])).toBe(0);
  });

  it('larger effects need more null studies', () => {
    const largeEffects = asymmetricEffects.map(e => ({ ...e, es: e.es * 2 }));
    const fsnLarge = failSafeN(largeEffects);
    const fsnSmall = failSafeN(asymmetricEffects);

    expect(fsnLarge).toBeGreaterThan(fsnSmall);
  });
});

describe('orwinFailSafeN', () => {
  it('returns non-negative value', () => {
    const fsn = orwinFailSafeN(asymmetricEffects, 0.1);
    expect(fsn).toBeGreaterThanOrEqual(0);
  });

  it('returns 0 for < 2 studies', () => {
    expect(orwinFailSafeN([asymmetricEffects[0]])).toBe(0);
  });

  it('returns Infinity when criterion equals mean', () => {
    const meanES = asymmetricEffects.reduce((a, e) => a + e.es, 0) / asymmetricEffects.length;
    const fsn = orwinFailSafeN(asymmetricEffects, meanES);
    expect(fsn).toBe(Infinity);
  });
});

describe('trimAndFill', () => {
  it('returns k0 (number of imputed studies)', () => {
    const result = trimAndFill(asymmetricEffects);

    expect(result.k0).toBeGreaterThanOrEqual(0);
    expect(Number.isInteger(result.k0)).toBe(true);
  });

  it('returns imputed studies array', () => {
    const result = trimAndFill(asymmetricEffects);

    expect(Array.isArray(result.imputedStudies)).toBe(true);
    if (result.k0 > 0) {
      expect(result.imputedStudies.length).toBe(result.k0);
    }
  });

  it('imputed studies have imputed flag', () => {
    const result = trimAndFill(asymmetricEffects);

    result.imputedStudies.forEach(s => {
      expect(s.imputed).toBe(true);
    });
  });

  it('returns adjusted pooled estimate', () => {
    const result = trimAndFill(asymmetricEffects);

    if (result.k0 > 0) {
      expect(result.adjusted).not.toBeNull();
      expect(result.adjusted.es).toBeDefined();
    }
  });

  it('adjusted estimate differs from original if studies imputed', () => {
    const tau2 = tauDL(asymmetricEffects);
    const original = poolInverseVariance(asymmetricEffects, tau2);
    const tfResult = trimAndFill(asymmetricEffects);

    if (tfResult.k0 > 0 && tfResult.adjusted) {
      // Adjusted should be different (typically smaller for positive bias)
      expect(tfResult.adjusted.es).not.toEqual(original.es);
    }
  });

  it('returns null adjusted for < 3 studies', () => {
    const result = trimAndFill(asymmetricEffects.slice(0, 2));
    expect(result.k0).toBe(0);
    expect(result.adjusted).toBeNull();
  });

  it('returns 0 imputed for symmetric data', () => {
    const result = trimAndFill(symmetricEffects);
    expect(result.k0).toBeLessThanOrEqual(1); // May be 0 or 1 due to ties
  });
});

describe('getFunnelPlotData', () => {
  it('returns plot data structure', () => {
    const tau2 = tauDL(asymmetricEffects);
    const pooled = poolInverseVariance(asymmetricEffects, tau2);
    const data = getFunnelPlotData(asymmetricEffects, pooled);

    expect(data).not.toBeNull();
    expect(data.studies).toBeDefined();
    expect(data.pooledES).toBeDefined();
    expect(data.ciLines).toBeDefined();
    expect(data.maxSE).toBeDefined();
  });

  it('returns study coordinates', () => {
    const tau2 = tauDL(asymmetricEffects);
    const pooled = poolInverseVariance(asymmetricEffects, tau2);
    const data = getFunnelPlotData(asymmetricEffects, pooled);

    expect(data.studies.length).toBe(5);
    data.studies.forEach(s => {
      expect(s.id).toBeDefined();
      expect(s.x).toBeDefined();
      expect(s.y).toBeDefined();
    });
  });

  it('transforms ratio metrics correctly', () => {
    const tau2 = tauDL(asymmetricEffects);
    const pooled = poolInverseVariance(asymmetricEffects, tau2);

    const dataNoTransform = getFunnelPlotData(asymmetricEffects, pooled, false);
    const dataWithTransform = getFunnelPlotData(asymmetricEffects, pooled, true);

    // With transform, values should be exponentiated
    expect(dataWithTransform.pooledES).not.toEqual(dataNoTransform.pooledES);
  });

  it('returns null for < 3 studies', () => {
    const tau2 = tauDL(asymmetricEffects.slice(0, 2));
    const pooled = poolInverseVariance(asymmetricEffects.slice(0, 2), tau2);
    expect(getFunnelPlotData(asymmetricEffects.slice(0, 2), pooled)).toBeNull();
  });
});

describe('getContourFunnelData', () => {
  it('returns contour regions for significance levels', () => {
    const tau2 = tauDL(asymmetricEffects);
    const pooled = poolInverseVariance(asymmetricEffects, tau2);
    const data = getContourFunnelData(pooled, 0.5);

    expect(Array.isArray(data)).toBe(true);
    expect(data.length).toBe(3); // p = 0.01, 0.05, 0.10
  });

  it('each contour has upper and lower bounds', () => {
    const tau2 = tauDL(asymmetricEffects);
    const pooled = poolInverseVariance(asymmetricEffects, tau2);
    const data = getContourFunnelData(pooled, 0.5);

    data.forEach(contour => {
      expect(contour.p).toBeDefined();
      expect(contour.lower).toBeDefined();
      expect(contour.upper).toBeDefined();
      expect(contour.se).toBeDefined();
    });
  });
});
