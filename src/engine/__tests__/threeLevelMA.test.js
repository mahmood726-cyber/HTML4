/**
 * Unit tests for Three-Level Meta-Analysis module
 * Tests hierarchical meta-analysis with nested effect sizes
 */

import { describe, it, expect } from 'vitest';
import {
  threeLevelMetaAnalysis,
  threeLevelMetaRegression,
  compareModels
} from '../threeLevelMA.js';

// Test dataset with hierarchical structure (effects nested within studies)
const hierarchicalStudies = [
  { id: 'S1_O1', studyId: 'S1', es: 0.30, vi: 0.04, excluded: false },
  { id: 'S1_O2', studyId: 'S1', es: 0.35, vi: 0.05, excluded: false },
  { id: 'S1_O3', studyId: 'S1', es: 0.28, vi: 0.04, excluded: false },
  { id: 'S2_O1', studyId: 'S2', es: 0.20, vi: 0.03, excluded: false },
  { id: 'S2_O2', studyId: 'S2', es: 0.25, vi: 0.04, excluded: false },
  { id: 'S3_O1', studyId: 'S3', es: 0.40, vi: 0.06, excluded: false },
  { id: 'S3_O2', studyId: 'S3', es: 0.45, vi: 0.05, excluded: false },
  { id: 'S4_O1', studyId: 'S4', es: 0.15, vi: 0.03, excluded: false },
  { id: 'S5_O1', studyId: 'S5', es: 0.50, vi: 0.07, excluded: false },
  { id: 'S5_O2', studyId: 'S5', es: 0.48, vi: 0.06, excluded: false }
];

describe('threeLevelMetaAnalysis', () => {
  it('calculates pooled effect estimate', () => {
    const result = threeLevelMetaAnalysis(hierarchicalStudies);

    expect(result.error).toBeUndefined();
    expect(result.mu.estimate).toBeGreaterThan(0);
    expect(result.mu.estimate).toBeLessThan(1);
    expect(Number.isFinite(result.mu.estimate)).toBe(true);
  });

  it('returns correct study and effect counts', () => {
    const result = threeLevelMetaAnalysis(hierarchicalStudies);

    expect(result.structure.nEffects).toBe(10);
    expect(result.structure.nClusters).toBe(5);
  });

  it('estimates both level-2 and level-3 variance', () => {
    const result = threeLevelMetaAnalysis(hierarchicalStudies);

    expect(result.varianceComponents.tau2_level2).toBeGreaterThanOrEqual(0);
    expect(result.varianceComponents.tau2_level3).toBeGreaterThanOrEqual(0);
    expect(result.varianceComponents.tau_level2).toBeGreaterThanOrEqual(0);
    expect(result.varianceComponents.tau_level3).toBeGreaterThanOrEqual(0);
  });

  it('calculates variance distribution (I² decomposition)', () => {
    const result = threeLevelMetaAnalysis(hierarchicalStudies);

    expect(result.heterogeneity.I2_total).toBeGreaterThanOrEqual(0);
    expect(result.heterogeneity.I2_level2).toBeGreaterThanOrEqual(0);
    expect(result.heterogeneity.I2_level3).toBeGreaterThanOrEqual(0);
  });

  it('provides confidence interval', () => {
    const result = threeLevelMetaAnalysis(hierarchicalStudies);

    expect(result.mu.ci.lo).toBeLessThan(result.mu.estimate);
    expect(result.mu.ci.hi).toBeGreaterThan(result.mu.estimate);
  });

  it('calculates test statistic and p-value', () => {
    const result = threeLevelMetaAnalysis(hierarchicalStudies);

    expect(Number.isFinite(result.mu.t)).toBe(true);
    expect(result.mu.p).toBeGreaterThanOrEqual(0);
    expect(result.mu.p).toBeLessThanOrEqual(1);
  });

  it('handles single effect per study correctly', () => {
    const singleEffectStudies = [
      { id: 'S1', studyId: 'S1', es: 0.30, vi: 0.04, excluded: false },
      { id: 'S2', studyId: 'S2', es: 0.20, vi: 0.03, excluded: false },
      { id: 'S3', studyId: 'S3', es: 0.40, vi: 0.06, excluded: false },
      { id: 'S4', studyId: 'S4', es: 0.15, vi: 0.03, excluded: false }
    ];

    const result = threeLevelMetaAnalysis(singleEffectStudies);

    expect(result.error).toBeUndefined();
  });

  it('returns error for insufficient studies', () => {
    const tooFew = hierarchicalStudies.slice(0, 2);
    const result = threeLevelMetaAnalysis(tooFew);

    expect(result.error).toBeDefined();
  });

  it('includes interpretation', () => {
    const result = threeLevelMetaAnalysis(hierarchicalStudies);

    expect(result.interpretation).toBeDefined();
    expect(typeof result.interpretation).toBe('string');
  });

  it('provides prediction interval', () => {
    const result = threeLevelMetaAnalysis(hierarchicalStudies);

    expect(result.predictionInterval).toBeDefined();
    expect(result.predictionInterval.lo).toBeLessThan(result.mu.estimate);
    expect(result.predictionInterval.hi).toBeGreaterThan(result.mu.estimate);
  });
});

describe('threeLevelMetaRegression', () => {
  const studiesWithCovariates = hierarchicalStudies.map((s, i) => ({
    ...s,
    year: 2010 + Math.floor(i / 2),
    intervention: i % 2 === 0 ? 1 : 0
  }));

  it('estimates regression coefficients', () => {
    const result = threeLevelMetaRegression(studiesWithCovariates, ['year']);

    expect(result.error).toBeUndefined();
    expect(result.coefficients.length).toBe(2);
    expect(result.coefficients[0].name).toBe('(Intercept)');
    expect(result.coefficients[1].name).toBe('year');
  });

  it('provides standard errors and p-values', () => {
    const result = threeLevelMetaRegression(studiesWithCovariates, ['year']);

    result.coefficients.forEach(coef => {
      expect(coef.se).toBeGreaterThan(0);
      expect(coef.p).toBeGreaterThanOrEqual(0);
      expect(coef.p).toBeLessThanOrEqual(1);
    });
  });

  it('handles multiple covariates', () => {
    const result = threeLevelMetaRegression(studiesWithCovariates, ['year', 'intervention']);

    expect(result.coefficients.length).toBe(3);
    expect(result.coefficients[2].name).toBe('intervention');
  });

  it('reports residual variance components', () => {
    const result = threeLevelMetaRegression(studiesWithCovariates, ['year']);

    expect(result.varianceComponents.tau2_level2).toBeGreaterThanOrEqual(0);
    expect(result.varianceComponents.tau2_level3).toBeGreaterThanOrEqual(0);
  });

  it('calculates R²', () => {
    const result = threeLevelMetaRegression(studiesWithCovariates, ['year']);

    expect(result.R2).toBeGreaterThanOrEqual(0);
    expect(result.R2).toBeLessThanOrEqual(100);
  });
});

describe('compareModels', () => {
  it('compares two-level and three-level models', () => {
    const result = compareModels(hierarchicalStudies);

    expect(result.twoLevel).toBeDefined();
    expect(result.threeLevel).toBeDefined();
  });

  it('provides AIC and BIC for model selection', () => {
    const result = compareModels(hierarchicalStudies);

    expect(Number.isFinite(result.twoLevel.AIC)).toBe(true);
    expect(Number.isFinite(result.threeLevel.AIC)).toBe(true);
    expect(Number.isFinite(result.twoLevel.BIC)).toBe(true);
    expect(Number.isFinite(result.threeLevel.BIC)).toBe(true);
  });

  it('calculates likelihood ratio test', () => {
    const result = compareModels(hierarchicalStudies);

    expect(result.comparison).toBeDefined();
    expect(result.comparison.LRT).toBeGreaterThanOrEqual(0);
    expect(result.comparison.pLRT).toBeGreaterThanOrEqual(0);
    expect(result.comparison.pLRT).toBeLessThanOrEqual(1);
  });

  it('provides model recommendation', () => {
    const result = compareModels(hierarchicalStudies);

    expect(result.recommendation).toBeDefined();
    expect(typeof result.recommendation).toBe('string');
  });

  it('includes variance components for both models', () => {
    const result = compareModels(hierarchicalStudies);

    expect(result.twoLevel.tau2).toBeGreaterThanOrEqual(0);
    expect(result.threeLevel.tau2_level2).toBeGreaterThanOrEqual(0);
    expect(result.threeLevel.tau2_level3).toBeGreaterThanOrEqual(0);
  });
});
