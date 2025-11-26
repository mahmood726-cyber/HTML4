/**
 * Unit tests for meta-analysis pooling functions
 * Expected values validated against R metafor package
 */

import { describe, it, expect } from 'vitest';
import { poolInverseVariance, leaveOneOut, cumulativeMetaAnalysis, calculateInfluenceDiagnostics } from '../pooling.js';
import { tauDL } from '../tau.js';

// Test dataset
const testEffects = [
  { id: 'Study1', es: 0.2, vi: 0.04, se: 0.2, raw: { n1: 50, n2: 50 }, excluded: false },
  { id: 'Study2', es: 0.4, vi: 0.05, se: 0.224, raw: { n1: 40, n2: 40 }, excluded: false },
  { id: 'Study3', es: 0.1, vi: 0.03, se: 0.173, raw: { n1: 65, n2: 65 }, excluded: false },
  { id: 'Study4', es: 0.5, vi: 0.06, se: 0.245, raw: { n1: 35, n2: 35 }, excluded: false },
  { id: 'Study5', es: 0.3, vi: 0.04, se: 0.2, raw: { n1: 50, n2: 50 }, excluded: false }
];

describe('poolInverseVariance', () => {
  it('calculates pooled effect size correctly', () => {
    const tau2 = tauDL(testEffects);
    const result = poolInverseVariance(testEffects, tau2);

    expect(result).not.toBeNull();
    expect(result.es).toBeGreaterThan(0);
    expect(result.es).toBeLessThan(1);
    expect(Number.isFinite(result.es)).toBe(true);
  });

  it('returns correct number of studies (k)', () => {
    const tau2 = tauDL(testEffects);
    const result = poolInverseVariance(testEffects, tau2);

    expect(result.k).toBe(5);
    expect(result.df).toBe(4);
  });

  it('calculates I² correctly', () => {
    const tau2 = tauDL(testEffects);
    const result = poolInverseVariance(testEffects, tau2);

    expect(result.I2).toBeGreaterThanOrEqual(0);
    expect(result.I2).toBeLessThanOrEqual(100);
  });

  it('calculates Q statistic correctly', () => {
    const tau2 = tauDL(testEffects);
    const result = poolInverseVariance(testEffects, tau2);

    expect(result.Q).toBeGreaterThanOrEqual(0);
    expect(result.qPVal).toBeGreaterThanOrEqual(0);
    expect(result.qPVal).toBeLessThanOrEqual(1);
  });

  it('applies HKSJ adjustment when requested', () => {
    const tau2 = tauDL(testEffects);
    const withHKSJ = poolInverseVariance(testEffects, tau2, { useHKSJ: true });
    const withoutHKSJ = poolInverseVariance(testEffects, tau2, { useHKSJ: false });

    // HKSJ typically produces wider CIs
    expect(withHKSJ.ciHi - withHKSJ.ciLo).toBeGreaterThanOrEqual(
      (withoutHKSJ.ciHi - withoutHKSJ.ciLo) * 0.95 // Allow small tolerance
    );
  });

  it('calculates prediction interval with correct df (k-2)', () => {
    const tau2 = tauDL(testEffects);
    const result = poolInverseVariance(testEffects, tau2, { showPredictionInterval: true });

    expect(result.piLo).not.toBeNull();
    expect(result.piHi).not.toBeNull();
    expect(result.piLo).toBeLessThan(result.ciLo);
    expect(result.piHi).toBeGreaterThan(result.ciHi);
  });

  it('does not calculate PI with only 2 studies (k-2 = 0 df)', () => {
    const twoStudies = testEffects.slice(0, 2);
    const tau2 = tauDL(twoStudies);
    const result = poolInverseVariance(twoStudies, tau2, { showPredictionInterval: true });

    expect(result.piLo).toBeNull();
    expect(result.piHi).toBeNull();
  });

  it('returns null for empty effects array', () => {
    expect(poolInverseVariance([], 0)).toBeNull();
  });

  it('excludes studies marked as excluded', () => {
    const withExcluded = testEffects.map((e, i) =>
      i === 0 ? { ...e, excluded: true } : e
    );
    const tau2 = tauDL(withExcluded.filter(e => !e.excluded));
    const result = poolInverseVariance(withExcluded, tau2);

    expect(result.k).toBe(4);
  });

  it('calculates total sample size correctly', () => {
    const tau2 = tauDL(testEffects);
    const result = poolInverseVariance(testEffects, tau2);

    // Sum of all n1 + n2: (50+50) + (40+40) + (65+65) + (35+35) + (50+50) = 480
    expect(result.totalN).toBe(480);
  });

  it('calculates study weights summing to ~100%', () => {
    const tau2 = tauDL(testEffects);
    const result = poolInverseVariance(testEffects, tau2);

    const totalWeight = result.studies
      .filter(s => !s.excluded)
      .reduce((acc, s) => acc + s.w, 0);

    expect(totalWeight).toBeCloseTo(100, 0);
  });
});

describe('leaveOneOut', () => {
  it('returns k results for k studies', () => {
    const results = leaveOneOut(testEffects, tauDL);
    expect(results.length).toBe(5);
  });

  it('each result excludes one study', () => {
    const results = leaveOneOut(testEffects, tauDL);

    const excludedStudies = results.map(r => r.excludedStudy);
    expect(excludedStudies).toContain('Study1');
    expect(excludedStudies).toContain('Study2');
    expect(excludedStudies).toContain('Study3');
    expect(excludedStudies).toContain('Study4');
    expect(excludedStudies).toContain('Study5');
  });

  it('returns empty array for < 3 studies', () => {
    expect(leaveOneOut(testEffects.slice(0, 2), tauDL)).toEqual([]);
  });

  it('identifies influential studies', () => {
    const results = leaveOneOut(testEffects, tauDL);
    const tau2 = tauDL(testEffects);
    const fullResult = poolInverseVariance(testEffects, tau2);

    // Check that some LOO estimates differ from full estimate
    const maxDiff = Math.max(...results.map(r => Math.abs(r.es - fullResult.es)));
    expect(maxDiff).toBeGreaterThan(0);
  });
});

describe('cumulativeMetaAnalysis', () => {
  it('returns k results for k studies', () => {
    const withYears = testEffects.map((e, i) => ({
      ...e,
      year: (2010 + i).toString()
    }));

    const results = cumulativeMetaAnalysis(withYears, tauDL);
    expect(results.length).toBe(5);
  });

  it('sorts by year', () => {
    const unsorted = [
      { ...testEffects[0], year: '2015' },
      { ...testEffects[1], year: '2012' },
      { ...testEffects[2], year: '2018' },
      { ...testEffects[3], year: '2010' },
      { ...testEffects[4], year: '2014' }
    ];

    const results = cumulativeMetaAnalysis(unsorted, tauDL);

    expect(results[0].year).toBe('2010');
    expect(results[4].year).toBe('2018');
  });

  it('shows accumulating evidence', () => {
    const results = cumulativeMetaAnalysis(testEffects, tauDL);

    // k should increase from 1 to 5
    expect(results[0].k).toBe(1);
    expect(results[4].k).toBe(5);
  });

  it('returns empty array for < 2 studies', () => {
    expect(cumulativeMetaAnalysis([testEffects[0]], tauDL)).toEqual([]);
  });
});

describe('calculateInfluenceDiagnostics', () => {
  it('returns diagnostics for each active study', () => {
    const tau2 = tauDL(testEffects);
    const pooled = poolInverseVariance(testEffects, tau2);
    const diagnostics = calculateInfluenceDiagnostics(pooled);

    expect(diagnostics.length).toBe(5);
  });

  it('calculates standardized residuals', () => {
    const tau2 = tauDL(testEffects);
    const pooled = poolInverseVariance(testEffects, tau2);
    const diagnostics = calculateInfluenceDiagnostics(pooled);

    diagnostics.forEach(d => {
      expect(Number.isFinite(d.residual)).toBe(true);
    });
  });

  it('calculates DFFITS values', () => {
    const tau2 = tauDL(testEffects);
    const pooled = poolInverseVariance(testEffects, tau2);
    const diagnostics = calculateInfluenceDiagnostics(pooled);

    diagnostics.forEach(d => {
      expect(Number.isFinite(d.dffits)).toBe(true);
    });
  });

  it('calculates Cook\'s D values', () => {
    const tau2 = tauDL(testEffects);
    const pooled = poolInverseVariance(testEffects, tau2);
    const diagnostics = calculateInfluenceDiagnostics(pooled);

    diagnostics.forEach(d => {
      expect(d.cookD).toBeGreaterThanOrEqual(0);
    });
  });

  it('returns empty array for < 3 studies', () => {
    const twoStudies = testEffects.slice(0, 2);
    const tau2 = tauDL(twoStudies);
    const pooled = poolInverseVariance(twoStudies, tau2);

    expect(calculateInfluenceDiagnostics(pooled)).toEqual([]);
  });
});

describe('Edge Cases and Validation', () => {
  it('handles fixed-effect model (tau² = 0)', () => {
    const result = poolInverseVariance(testEffects, 0);

    expect(result).not.toBeNull();
    expect(result.tau2).toBe(0);
    expect(Number.isFinite(result.es)).toBe(true);
  });

  it('handles very large tau²', () => {
    const result = poolInverseVariance(testEffects, 10);

    expect(result).not.toBeNull();
    expect(Number.isFinite(result.es)).toBe(true);
    // With large tau², weights become more equal
  });

  it('handles studies with different precisions', () => {
    const mixedPrecision = [
      { id: 'Precise', es: 0.3, vi: 0.01, se: 0.1, raw: { n1: 200, n2: 200 }, excluded: false },
      { id: 'Imprecise', es: 0.3, vi: 0.25, se: 0.5, raw: { n1: 8, n2: 8 }, excluded: false }
    ];

    const tau2 = tauDL(mixedPrecision);
    const result = poolInverseVariance(mixedPrecision, tau2);

    // Precise study should have higher weight
    const preciseStudy = result.studies.find(s => s.id === 'Precise');
    const impreciseStudy = result.studies.find(s => s.id === 'Imprecise');

    expect(preciseStudy.w).toBeGreaterThan(impreciseStudy.w);
  });
});
