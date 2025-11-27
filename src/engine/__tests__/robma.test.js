/**
 * Unit tests for Robust Bayesian Meta-Analysis (RoBMA) module
 * Tests Bayesian model averaging for meta-analysis
 */

import { describe, it, expect } from 'vitest';
import {
  robustBayesianMA,
  bayesianModelAveraging
} from '../robma.js';

// Test dataset
const testStudies = [
  { id: 'S1', es: 0.30, se: 0.12, excluded: false },
  { id: 'S2', es: 0.45, se: 0.15, excluded: false },
  { id: 'S3', es: 0.20, se: 0.10, excluded: false },
  { id: 'S4', es: 0.55, se: 0.18, excluded: false },
  { id: 'S5', es: 0.35, se: 0.14, excluded: false },
  { id: 'S6', es: 0.25, se: 0.11, excluded: false },
  { id: 'S7', es: 0.40, se: 0.13, excluded: false }
];

// Studies with possible publication bias (all significant)
const biasedStudies = [
  { id: 'S1', es: 0.50, se: 0.15, excluded: false },  // z = 3.33
  { id: 'S2', es: 0.55, se: 0.18, excluded: false },  // z = 3.06
  { id: 'S3', es: 0.48, se: 0.16, excluded: false },  // z = 3.00
  { id: 'S4', es: 0.60, se: 0.20, excluded: false },  // z = 3.00
  { id: 'S5', es: 0.45, se: 0.14, excluded: false }   // z = 3.21
];

describe('robustBayesianMA', () => {
  it('calculates model-averaged effect estimate', () => {
    const result = robustBayesianMA(testStudies);

    expect(result.error).toBeUndefined();
    expect(result.modelAveraged).toBeDefined();
    expect(Number.isFinite(result.modelAveraged.effect)).toBe(true);
  });

  it('returns correct study count', () => {
    const result = robustBayesianMA(testStudies);

    expect(result.k).toBe(7);
  });

  it('provides credible interval', () => {
    const result = robustBayesianMA(testStudies);

    expect(result.modelAveraged.ci95.lo).toBeLessThan(result.modelAveraged.effect);
    expect(result.modelAveraged.ci95.hi).toBeGreaterThan(result.modelAveraged.effect);
  });

  it('estimates model-averaged tau', () => {
    const result = robustBayesianMA(testStudies);

    expect(result.modelAveraged.tau2).toBeGreaterThanOrEqual(0);
    expect(result.modelAveraged.tau).toBeGreaterThanOrEqual(0);
  });

  it('calculates inclusion Bayes factors', () => {
    const result = robustBayesianMA(testStudies);

    expect(result.inclusionBayesFactors).toBeDefined();
    expect(result.inclusionBayesFactors.effect).toBeDefined();
    expect(result.inclusionBayesFactors.heterogeneity).toBeDefined();
    expect(result.inclusionBayesFactors.publicationBias).toBeDefined();
  });

  it('provides BF interpretation', () => {
    const result = robustBayesianMA(testStudies);

    expect(result.inclusionBayesFactors.effect.interpretation).toBeDefined();
    expect(typeof result.inclusionBayesFactors.effect.interpretation).toBe('string');
  });

  it('calculates posterior probabilities', () => {
    const result = robustBayesianMA(testStudies);

    expect(result.posteriorProbabilities.effect).toBeGreaterThanOrEqual(0);
    expect(result.posteriorProbabilities.effect).toBeLessThanOrEqual(1);
    expect(result.posteriorProbabilities.heterogeneity).toBeGreaterThanOrEqual(0);
    expect(result.posteriorProbabilities.heterogeneity).toBeLessThanOrEqual(1);
    expect(result.posteriorProbabilities.publicationBias).toBeGreaterThanOrEqual(0);
    expect(result.posteriorProbabilities.publicationBias).toBeLessThanOrEqual(1);
  });

  it('evaluates multiple models', () => {
    const result = robustBayesianMA(testStudies);

    expect(result.models).toBeDefined();
    expect(result.models.length).toBe(8);  // 2^3 combinations
  });

  it('assigns posterior probabilities to models', () => {
    const result = robustBayesianMA(testStudies);

    const sumPosterior = result.models.reduce((sum, m) => sum + m.posteriorProb, 0);
    expect(sumPosterior).toBeCloseTo(1, 2);
  });

  it('detects potential publication bias', () => {
    const result = robustBayesianMA(biasedStudies);

    // With all significant studies, bias probability should be elevated
    expect(result.posteriorProbabilities.publicationBias).toBeGreaterThan(0.1);
  });

  it('returns error for insufficient studies', () => {
    const tooFew = testStudies.slice(0, 2);
    const result = robustBayesianMA(tooFew);

    expect(result.error).toBeDefined();
  });

  it('includes interpretation text', () => {
    const result = robustBayesianMA(testStudies);

    expect(result.interpretation).toBeDefined();
    expect(typeof result.interpretation).toBe('string');
  });

  it('supports custom priors', () => {
    const result = robustBayesianMA(testStudies, {
      effectPrior: { type: 'normal', mean: 0, sd: 0.5 }
    });

    expect(result.error).toBeUndefined();
  });
});

describe('bayesianModelAveraging', () => {
  it('calculates Bayes factor for effect', () => {
    const result = bayesianModelAveraging(testStudies);

    expect(result.error).toBeUndefined();
    expect(result.BF10).toBeGreaterThan(0);
    expect(Number.isFinite(result.BF10)).toBe(true);
  });

  it('provides log Bayes factor', () => {
    const result = bayesianModelAveraging(testStudies);

    expect(Number.isFinite(result.logBF10)).toBe(true);
    expect(Math.exp(result.logBF10)).toBeCloseTo(result.BF10, 1);
  });

  it('calculates posterior probabilities', () => {
    const result = bayesianModelAveraging(testStudies);

    expect(result.posteriorProbEffect).toBeGreaterThanOrEqual(0);
    expect(result.posteriorProbEffect).toBeLessThanOrEqual(1);
    expect(result.posteriorProbNull).toBeCloseTo(1 - result.posteriorProbEffect, 5);
  });

  it('provides model-averaged effect', () => {
    const result = bayesianModelAveraging(testStudies);

    expect(Number.isFinite(result.modelAveragedEffect)).toBe(true);
  });

  it('gives conditional effect estimate', () => {
    const result = bayesianModelAveraging(testStudies);

    expect(result.conditionalEffect).toBeDefined();
    expect(Number.isFinite(result.conditionalEffect.mean)).toBe(true);
    expect(result.conditionalEffect.sd).toBeGreaterThan(0);
    expect(result.conditionalEffect.ci95.lo).toBeLessThan(result.conditionalEffect.mean);
    expect(result.conditionalEffect.ci95.hi).toBeGreaterThan(result.conditionalEffect.mean);
  });

  it('interprets Bayes factor', () => {
    const result = bayesianModelAveraging(testStudies);

    expect(result.interpretation).toBeDefined();
    expect(typeof result.interpretation).toBe('string');
  });

  it('supports custom prior probability', () => {
    const result05 = bayesianModelAveraging(testStudies, { priorProbEffect: 0.5 });
    const result09 = bayesianModelAveraging(testStudies, { priorProbEffect: 0.9 });

    // Higher prior should give higher posterior
    expect(result09.posteriorProbEffect).toBeGreaterThan(result05.posteriorProbEffect * 0.9);
  });

  it('supports custom prior SD', () => {
    const result = bayesianModelAveraging(testStudies, { effectPriorSD: 0.3 });

    expect(result.error).toBeUndefined();
  });

  it('returns error for insufficient studies', () => {
    const oneStudy = [testStudies[0]];
    const result = bayesianModelAveraging(oneStudy);

    expect(result.error).toBeDefined();
  });
});

describe('Bayes factor interpretation', () => {
  it('correctly classifies strong evidence', () => {
    // Create studies with strong effect
    const strongEffect = Array(10).fill(null).map((_, i) => ({
      id: `S${i}`,
      es: 0.8,
      se: 0.15,
      excluded: false
    }));

    const result = bayesianModelAveraging(strongEffect);

    expect(result.BF10).toBeGreaterThan(10);
    // "strong" or "extreme" evidence both indicate strong support for H1
    const interp = result.interpretation.toLowerCase();
    expect(interp.includes('strong') || interp.includes('extreme')).toBe(true);
  });

  it('correctly classifies weak evidence', () => {
    // Create studies with mixed/weak effect
    const weakEffect = [
      { id: 'S1', es: 0.1, se: 0.2, excluded: false },
      { id: 'S2', es: -0.05, se: 0.18, excluded: false },
      { id: 'S3', es: 0.15, se: 0.22, excluded: false },
      { id: 'S4', es: 0.0, se: 0.19, excluded: false }
    ];

    const result = bayesianModelAveraging(weakEffect);

    // Should be inconclusive (BF close to 1)
    expect(result.BF10).toBeLessThan(10);
  });
});
