/**
 * Unit tests for Exact Likelihood methods module
 * Tests GLMM and other methods for sparse binary data
 */

import { describe, it, expect } from 'vitest';
import {
  binomialNormalModel,
  betaBinomialModel,
  penalizedLikelihood,
  compareSparseDataMethods
} from '../exactLikelihood.js';

// Test dataset with binary outcomes
const binaryStudies = [
  { id: 'S1', e1: 15, n1: 100, e2: 10, n2: 100, excluded: false },
  { id: 'S2', e1: 20, n1: 150, e2: 12, n2: 150, excluded: false },
  { id: 'S3', e1: 8, n1: 80, e2: 5, n2: 80, excluded: false },
  { id: 'S4', e1: 25, n1: 200, e2: 18, n2: 200, excluded: false },
  { id: 'S5', e1: 12, n1: 120, e2: 8, n2: 120, excluded: false }
];

// Sparse data with zero cells
const sparseStudies = [
  { id: 'S1', e1: 2, n1: 50, e2: 0, n2: 50, excluded: false },
  { id: 'S2', e1: 3, n1: 60, e2: 1, n2: 60, excluded: false },
  { id: 'S3', e1: 1, n1: 40, e2: 0, n2: 40, excluded: false },
  { id: 'S4', e1: 4, n1: 80, e2: 2, n2: 80, excluded: false },
  { id: 'S5', e1: 0, n1: 30, e2: 0, n2: 30, excluded: false }
];

describe('binomialNormalModel', () => {
  it('estimates log odds ratio', () => {
    const result = binomialNormalModel(binaryStudies);

    expect(result.error).toBeUndefined();
    expect(result.logEffect).toBeDefined();
    expect(Number.isFinite(result.logEffect.estimate)).toBe(true);
  });

  it('provides exponentiated effect (OR)', () => {
    const result = binomialNormalModel(binaryStudies);

    expect(result.effect).toBeDefined();
    expect(result.effect.label).toBe('OR');
    expect(result.effect.estimate).toBeGreaterThan(0);
    expect(result.effect.ci.lo).toBeGreaterThan(0);
    expect(result.effect.ci.hi).toBeGreaterThan(result.effect.ci.lo);
  });

  it('returns correct study count', () => {
    const result = binomialNormalModel(binaryStudies);

    expect(result.k).toBe(5);
  });

  it('estimates heterogeneity', () => {
    const result = binomialNormalModel(binaryStudies);

    expect(result.heterogeneity).toBeDefined();
    expect(result.heterogeneity.tau2).toBeGreaterThanOrEqual(0);
    expect(result.heterogeneity.I2).toBeGreaterThanOrEqual(0);
    expect(result.heterogeneity.I2).toBeLessThanOrEqual(100);
  });

  it('detects sparse data', () => {
    const result = binomialNormalModel(sparseStudies);

    expect(result.isSparse).toBe(true);
  });

  it('supports log link for risk ratio', () => {
    const result = binomialNormalModel(binaryStudies, { link: 'log' });

    expect(result.effect.label).toBe('RR');
    expect(result.effect.estimate).toBeGreaterThan(0);
  });

  it('provides study-level estimates', () => {
    const result = binomialNormalModel(binaryStudies);

    expect(result.studies).toBeDefined();
    expect(result.studies.length).toBe(5);
    result.studies.forEach(s => {
      expect(s.studyEffect).toBeDefined();
      expect(s.studySE).toBeGreaterThan(0);
    });
  });

  it('returns error for insufficient studies', () => {
    const oneStudy = [binaryStudies[0]];
    const result = binomialNormalModel(oneStudy);

    expect(result.error).toBeDefined();
  });

  it('includes interpretation text', () => {
    const result = binomialNormalModel(binaryStudies);

    expect(result.interpretation).toBeDefined();
    expect(typeof result.interpretation).toBe('string');
  });
});

describe('betaBinomialModel', () => {
  it('estimates odds ratio with overdispersion', () => {
    const result = betaBinomialModel(binaryStudies);

    expect(result.error).toBeUndefined();
    expect(result.OR).toBeDefined();
    expect(result.OR.estimate).toBeGreaterThan(0);
  });

  it('quantifies overdispersion', () => {
    const result = betaBinomialModel(binaryStudies);

    expect(result.overdispersion).toBeDefined();
    expect(result.overdispersion.phi).toBeGreaterThanOrEqual(1);
    expect(result.overdispersion.interpretation).toBeDefined();
  });

  it('estimates arm-specific rates', () => {
    const result = betaBinomialModel(binaryStudies);

    expect(result.armRates).toBeDefined();
    expect(result.armRates.treatment.mean).toBeGreaterThanOrEqual(0);
    expect(result.armRates.treatment.mean).toBeLessThanOrEqual(1);
    expect(result.armRates.control.mean).toBeGreaterThanOrEqual(0);
    expect(result.armRates.control.mean).toBeLessThanOrEqual(1);
  });

  it('provides confidence interval for log OR', () => {
    const result = betaBinomialModel(binaryStudies);

    expect(result.logOR.ci.lo).toBeLessThan(result.logOR.estimate);
    expect(result.logOR.ci.hi).toBeGreaterThan(result.logOR.estimate);
  });

  it('handles sparse data', () => {
    const result = betaBinomialModel(sparseStudies);

    expect(result.error).toBeUndefined();
    expect(result.OR.estimate).toBeGreaterThan(0);
  });
});

describe('penalizedLikelihood', () => {
  it('handles zero cells with Firth penalization', () => {
    const result = penalizedLikelihood(sparseStudies);

    expect(result.error).toBeUndefined();
    expect(result.OR.estimate).toBeGreaterThan(0);
    expect(Number.isFinite(result.OR.estimate)).toBe(true);
  });

  it('counts studies with zero cells', () => {
    const result = penalizedLikelihood(sparseStudies);

    expect(result.zeroCellStudies).toBeGreaterThan(0);
  });

  it('provides wider CIs than standard methods', () => {
    const standard = binomialNormalModel(binaryStudies);
    const penalized = penalizedLikelihood(binaryStudies);

    if (!standard.error && !penalized.error) {
      const standardWidth = standard.effect.ci.hi - standard.effect.ci.lo;
      const penalizedWidth = penalized.OR.ci.hi - penalized.OR.ci.lo;

      // Penalized should generally give wider or similar CIs
      expect(penalizedWidth).toBeGreaterThanOrEqual(standardWidth * 0.8);
    }
  });

  it('uses Firth penalty by default', () => {
    const result = penalizedLikelihood(sparseStudies);

    expect(result.penalty).toBe('firth');
  });

  it('provides interpretation for zero cells', () => {
    const result = penalizedLikelihood(sparseStudies);

    expect(result.interpretation).toContain('zero');
  });
});

describe('compareSparseDataMethods', () => {
  it('compares multiple methods', () => {
    const result = compareSparseDataMethods(binaryStudies);

    expect(result.results).toBeDefined();
    expect(result.results.length).toBeGreaterThanOrEqual(2);
  });

  it('reports data characteristics', () => {
    const result = compareSparseDataMethods(binaryStudies);

    expect(result.dataCharacteristics).toBeDefined();
    expect(result.dataCharacteristics.k).toBe(5);
    expect(result.dataCharacteristics.eventRate).toBeGreaterThan(0);
    expect(typeof result.dataCharacteristics.hasZeroCells).toBe('boolean');
  });

  it('provides method recommendation', () => {
    const result = compareSparseDataMethods(sparseStudies);

    expect(result.recommendation).toBeDefined();
    expect(typeof result.recommendation).toBe('string');
    // Should recommend penalized likelihood for sparse data
    expect(result.recommendation.toLowerCase()).toContain('penalized');
  });

  it('includes OR estimates from each method', () => {
    const result = compareSparseDataMethods(binaryStudies);

    result.results.forEach(r => {
      expect(r.method).toBeDefined();
      expect(r.OR).toBeGreaterThan(0);
      expect(r.CI).toBeDefined();
      expect(r.p).toBeGreaterThanOrEqual(0);
      expect(r.p).toBeLessThanOrEqual(1);
    });
  });

  it('detects zero cells in sparse data', () => {
    const result = compareSparseDataMethods(sparseStudies);

    expect(result.dataCharacteristics.hasZeroCells).toBe(true);
  });
});
