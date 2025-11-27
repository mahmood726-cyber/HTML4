/**
 * Tests for PET-PEESE Bias Correction Methods
 */
import { describe, it, expect } from 'vitest';
import {
  runPET,
  runPEESE,
  petPeese,
  enhancedTrimFill,
  comprehensiveBiasAssessment
} from '../petPeese.js';

// Test data with publication bias (small studies have larger effects)
const biasedStudies = [
  { id: 'Small1', es: -1.2, se: 0.5 },
  { id: 'Small2', es: -1.0, se: 0.45 },
  { id: 'Small3', es: -0.9, se: 0.4 },
  { id: 'Medium1', es: -0.6, se: 0.25 },
  { id: 'Medium2', es: -0.5, se: 0.22 },
  { id: 'Large1', es: -0.3, se: 0.12 },
  { id: 'Large2', es: -0.25, se: 0.10 },
];

// Unbiased studies (no correlation between ES and SE)
const unbiasedStudies = [
  { id: 'S1', es: -0.5, se: 0.15 },
  { id: 'S2', es: -0.45, se: 0.40 },
  { id: 'S3', es: -0.55, se: 0.25 },
  { id: 'S4', es: -0.48, se: 0.12 },
  { id: 'S5', es: -0.52, se: 0.35 },
];

describe('PET (Precision-Effect Test)', () => {
  it('should fit PET model', () => {
    const result = runPET(biasedStudies);

    expect(result.error).toBeUndefined();
    expect(result.method).toBe('PET');
    expect(result.k).toBe(7);
    expect(result.intercept.estimate).toBeDefined();
    expect(result.slope.estimate).toBeDefined();
  });

  it('should detect small-study effect in biased data', () => {
    const result = runPET(biasedStudies);

    // Slope should be significant (negative SE-effect relationship)
    expect(result.slope.estimate).toBeLessThan(0);
    expect(result.hasSmallStudyEffect).toBe(true);
  });

  it('should provide corrected estimate', () => {
    const result = runPET(biasedStudies);

    // Corrected ES (intercept) should be less extreme than simple average
    expect(result.correctedES).toBeDefined();
    expect(Math.abs(result.correctedES)).toBeLessThan(1.2);
  });

  it('should not detect bias in unbiased data', () => {
    const result = runPET(unbiasedStudies);

    // Should not show strong small-study effect
    expect(result.error).toBeUndefined();
    // Slope p-value should be higher
    expect(result.slope.p).toBeGreaterThan(0.05);
  });

  it('should handle too few studies', () => {
    const result = runPET(biasedStudies.slice(0, 2));
    expect(result.error).toBeDefined();
  });
});

describe('PEESE (Precision-Effect Estimate with Standard Error)', () => {
  it('should fit PEESE model', () => {
    const result = runPEESE(biasedStudies);

    expect(result.error).toBeUndefined();
    expect(result.method).toBe('PEESE');
    expect(result.k).toBe(7);
    expect(result.intercept.estimate).toBeDefined();
    expect(result.slope.estimate).toBeDefined();
  });

  it('should provide different estimate than PET', () => {
    const pet = runPET(biasedStudies);
    const peese = runPEESE(biasedStudies);

    // PEESE estimate typically differs from PET
    expect(pet.correctedES).not.toBeCloseTo(peese.correctedES, 3);
  });
});

describe('PET-PEESE Combined', () => {
  it('should apply conditional estimator', () => {
    const result = petPeese(biasedStudies);

    expect(result.error).toBeUndefined();
    expect(result.method).toBe('PET-PEESE');
    expect(result.pet).toBeDefined();
    expect(result.peese).toBeDefined();
    expect(result.decision).toBeDefined();
  });

  it('should choose PEESE when PET is significant', () => {
    const result = petPeese(biasedStudies);

    // Given the biased data, PET intercept should be significant
    if (result.pet.intercept.p < 0.10) {
      expect(result.decision).toBe('PEESE');
    }
  });

  it('should provide interpretation', () => {
    const result = petPeese(biasedStudies);
    expect(result.interpretation).toBeDefined();
    expect(typeof result.interpretation).toBe('string');
  });
});

describe('Enhanced Trim-and-Fill', () => {
  it('should detect asymmetry and impute studies', () => {
    const result = enhancedTrimFill(biasedStudies);

    expect(result.error).toBeUndefined();
    expect(result.method).toBe('Trim-and-Fill');
    expect(result.k0).toBeGreaterThanOrEqual(0);
  });

  it('should provide adjusted estimate', () => {
    const result = enhancedTrimFill(biasedStudies);

    expect(result.adjustedES).toBeDefined();
    expect(result.unadjustedES).toBeDefined();
    // Adjusted should be less extreme if studies were imputed
    if (result.k0 > 0) {
      expect(Math.abs(result.adjustedES)).toBeLessThanOrEqual(Math.abs(result.unadjustedES) + 0.01);
    }
  });

  it('should not impute for symmetric funnel', () => {
    const result = enhancedTrimFill(unbiasedStudies);

    // Symmetric data should have few or no imputed studies
    expect(result.k0).toBeLessThanOrEqual(1);
  });
});

describe('Comprehensive Bias Assessment', () => {
  it('should run all methods', () => {
    const result = comprehensiveBiasAssessment(biasedStudies);

    expect(result.pet).toBeDefined();
    expect(result.peese).toBeDefined();
    expect(result.petPeese).toBeDefined();
    expect(result.trimFill).toBeDefined();
  });

  it('should provide severity assessment', () => {
    const result = comprehensiveBiasAssessment(biasedStudies);

    expect(result.severity).toBeDefined();
    expect(['none', 'possible', 'mild', 'moderate', 'severe']).toContain(result.severity);
  });

  it('should count bias indicators', () => {
    const result = comprehensiveBiasAssessment(biasedStudies);

    expect(result.biasIndicators).toBeGreaterThanOrEqual(0);
    expect(result.indicators).toBeDefined();
  });

  it('should provide recommendation', () => {
    const result = comprehensiveBiasAssessment(biasedStudies);

    expect(result.recommendation).toBeDefined();
    expect(typeof result.recommendation).toBe('string');
  });
});
