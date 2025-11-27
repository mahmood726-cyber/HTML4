/**
 * Unit tests for Multivariate Meta-Analysis module
 * Tests correlated outcomes analysis using Riley's method
 */

import { describe, it, expect } from 'vitest';
import {
  rileyMultivariateMA,
  fullMultivariateMA,
  compareUnivariateMultivariate
} from '../multivariateMA.js';

// Test dataset in wide format (required by the API)
const wideFormatStudies = [
  { id: 'S1', sens: 0.85, sens_se: 0.14, spec: 0.80, spec_se: 0.17, excluded: false },
  { id: 'S2', sens: 0.82, sens_se: 0.16, spec: 0.78, spec_se: 0.17, excluded: false },
  { id: 'S3', sens: 0.88, sens_se: 0.13, spec: 0.82, spec_se: 0.16, excluded: false },
  { id: 'S4', sens: 0.79, sens_se: 0.17, spec: 0.75, spec_se: 0.19, excluded: false },
  { id: 'S5', sens: 0.90, sens_se: 0.12, spec: 0.85, spec_se: 0.14, excluded: false }
];

const outcomeNames = ['sens', 'spec'];

describe('rileyMultivariateMA', () => {
  it('estimates pooled effects for each outcome', () => {
    const result = rileyMultivariateMA(wideFormatStudies, outcomeNames);

    expect(result.error).toBeUndefined();
    expect(result.outcomes).toBeDefined();
    expect(result.outcomes.length).toBe(2);
  });

  it('provides effect estimates with confidence intervals', () => {
    const result = rileyMultivariateMA(wideFormatStudies, outcomeNames);

    result.outcomes.forEach(outcome => {
      expect(outcome.es).toBeGreaterThan(0);
      expect(outcome.se).toBeGreaterThan(0);
      expect(outcome.ci.lo).toBeLessThan(outcome.es);
      expect(outcome.ci.hi).toBeGreaterThan(outcome.es);
    });
  });

  it('estimates between-outcome correlation', () => {
    const result = rileyMultivariateMA(wideFormatStudies, outcomeNames);

    expect(result.betweenStudyCorrelation).toBeDefined();
    expect(result.betweenStudyCorrelation).toBeGreaterThanOrEqual(-1);
    expect(result.betweenStudyCorrelation).toBeLessThanOrEqual(1);
  });

  it('calculates heterogeneity for each outcome', () => {
    const result = rileyMultivariateMA(wideFormatStudies, outcomeNames);

    result.outcomes.forEach(outcome => {
      expect(outcome.tau2).toBeGreaterThanOrEqual(0);
      expect(outcome.I2).toBeGreaterThanOrEqual(0);
      expect(outcome.I2).toBeLessThanOrEqual(100);
    });
  });

  it('reports borrowing of strength metric', () => {
    const result = rileyMultivariateMA(wideFormatStudies, outcomeNames);

    expect(result.borrowingOfStrength).toBeDefined();
    expect(result.borrowingOfStrength.metrics).toBeDefined();
  });

  it('returns error for insufficient outcomes', () => {
    const result = rileyMultivariateMA(wideFormatStudies, ['sens']);

    expect(result.error).toBeDefined();
  });

  it('returns error for insufficient data', () => {
    const tooFew = wideFormatStudies.slice(0, 1);
    const result = rileyMultivariateMA(tooFew, outcomeNames);

    expect(result.error).toBeDefined();
  });

  it('includes method name', () => {
    const result = rileyMultivariateMA(wideFormatStudies, outcomeNames);

    expect(result.method).toBe('Riley Multivariate Meta-Analysis');
  });

  it('provides joint test', () => {
    const result = rileyMultivariateMA(wideFormatStudies, outcomeNames);

    expect(result.jointTest).toBeDefined();
    expect(result.jointTest.chiSq).toBeGreaterThanOrEqual(0);
    expect(result.jointTest.p).toBeGreaterThanOrEqual(0);
    expect(result.jointTest.p).toBeLessThanOrEqual(1);
  });
});

describe('fullMultivariateMA', () => {
  it('estimates pooled effects when all outcomes present', () => {
    const result = fullMultivariateMA(wideFormatStudies, outcomeNames);

    expect(result.error).toBeUndefined();
    expect(result.outcomes).toBeDefined();
    expect(result.outcomes.length).toBe(2);
  });

  it('provides method name', () => {
    const result = fullMultivariateMA(wideFormatStudies, outcomeNames);

    expect(result.method).toBe('Full Multivariate Meta-Analysis');
  });

  it('provides between-study covariance matrix', () => {
    const result = fullMultivariateMA(wideFormatStudies, outcomeNames);

    if (!result.error) {
      expect(result.betweenStudyCovariance).toBeDefined();
      expect(result.betweenStudyCovariance.length).toBe(2);
    }
  });

  it('provides between-study correlation matrix', () => {
    const result = fullMultivariateMA(wideFormatStudies, outcomeNames);

    if (!result.error) {
      expect(result.betweenStudyCorrelation).toBeDefined();
      expect(result.betweenStudyCorrelation.length).toBe(2);
    }
  });
});

describe('compareUnivariateMultivariate', () => {
  it('compares univariate and multivariate analyses', () => {
    const result = compareUnivariateMultivariate(wideFormatStudies, outcomeNames);

    expect(result.univariate).toBeDefined();
    expect(result.multivariate).toBeDefined();
  });

  it('reports efficiency gains from multivariate approach', () => {
    const result = compareUnivariateMultivariate(wideFormatStudies, outcomeNames);

    expect(result.comparison).toBeDefined();
    result.comparison.forEach(comp => {
      expect(comp.univariate.se).toBeGreaterThan(0);
      expect(comp.multivariate.se).toBeGreaterThan(0);
    });
  });

  it('indicates when multivariate approach is beneficial', () => {
    const result = compareUnivariateMultivariate(wideFormatStudies, outcomeNames);

    expect(result.recommendation).toBeDefined();
    expect(typeof result.recommendation).toBe('string');
  });

  it('provides summary statistics', () => {
    const result = compareUnivariateMultivariate(wideFormatStudies, outcomeNames);

    expect(result.summary).toBeDefined();
    expect(result.summary.avgSEReduction).toBeDefined();
  });
});
