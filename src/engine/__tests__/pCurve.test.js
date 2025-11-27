/**
 * Unit tests for P-Curve and P-Uniform module
 * Tests publication bias methods based on p-value distribution
 */

import { describe, it, expect } from 'vitest';
import {
  pCurveAnalysis,
  pUniform,
  pUniformStar,
  compareSelectionModels
} from '../pCurve.js';

// Test dataset with p-values
const studiesWithP = [
  { id: 'S1', es: 0.45, se: 0.15, p: 0.003, excluded: false },
  { id: 'S2', es: 0.38, se: 0.14, p: 0.007, excluded: false },
  { id: 'S3', es: 0.52, se: 0.18, p: 0.004, excluded: false },
  { id: 'S4', es: 0.30, se: 0.12, p: 0.012, excluded: false },
  { id: 'S5', es: 0.42, se: 0.16, p: 0.009, excluded: false },
  { id: 'S6', es: 0.35, se: 0.13, p: 0.007, excluded: false },
  { id: 'S7', es: 0.48, se: 0.17, p: 0.005, excluded: false }
];

// Studies without explicit p-values (will be calculated)
const studiesES = [
  { id: 'S1', es: 0.50, se: 0.15, excluded: false },  // z=3.33, p=0.0009
  { id: 'S2', es: 0.42, se: 0.14, excluded: false },  // z=3.00, p=0.0027
  { id: 'S3', es: 0.55, se: 0.18, excluded: false },  // z=3.06, p=0.0022
  { id: 'S4', es: 0.38, se: 0.12, excluded: false },  // z=3.17, p=0.0015
  { id: 'S5', es: 0.45, se: 0.16, excluded: false },  // z=2.81, p=0.0050
  { id: 'S6', es: 0.40, se: 0.13, excluded: false }   // z=3.08, p=0.0021
];

// Studies suggestive of p-hacking (many just-significant results)
const pHackedStudies = [
  { id: 'S1', es: 0.30, se: 0.15, p: 0.045, excluded: false },
  { id: 'S2', es: 0.28, se: 0.14, p: 0.046, excluded: false },
  { id: 'S3', es: 0.32, se: 0.16, p: 0.044, excluded: false },
  { id: 'S4', es: 0.29, se: 0.15, p: 0.048, excluded: false },
  { id: 'S5', es: 0.31, se: 0.16, p: 0.047, excluded: false }
];

describe('pCurveAnalysis', () => {
  it('analyzes p-value distribution', () => {
    const result = pCurveAnalysis(studiesWithP);

    expect(result.error).toBeUndefined();
    expect(result.method).toBe('P-Curve Analysis');
    expect(result.k).toBeGreaterThanOrEqual(3);
  });

  it('returns extracted p-values', () => {
    const result = pCurveAnalysis(studiesWithP);

    expect(result.pValues).toBeDefined();
    expect(result.pValues.length).toBe(7);
    result.pValues.forEach(p => {
      expect(p).toBeGreaterThan(0);
      expect(p).toBeLessThan(0.05);
    });
  });

  it('performs right-skew test', () => {
    const result = pCurveAnalysis(studiesWithP);

    expect(result.tests.rightSkew).toBeDefined();
    expect(result.tests.rightSkew.belowHalf).toBeGreaterThanOrEqual(0);
    expect(result.tests.rightSkew.aboveHalf).toBeGreaterThanOrEqual(0);
    expect(result.tests.rightSkew.binomialP).toBeGreaterThanOrEqual(0);
    expect(result.tests.rightSkew.binomialP).toBeLessThanOrEqual(1);
  });

  it('performs flatness test', () => {
    const result = pCurveAnalysis(studiesWithP);

    expect(result.tests.flatness).toBeDefined();
    expect(Number.isFinite(result.tests.flatness.stoufferZ)).toBe(true);
    expect(result.tests.flatness.stoufferP).toBeGreaterThanOrEqual(0);
    expect(result.tests.flatness.stoufferP).toBeLessThanOrEqual(1);
  });

  it('provides p-value histogram', () => {
    const result = pCurveAnalysis(studiesWithP);

    expect(result.histogram).toBeDefined();
    expect(result.histogram.length).toBe(5);  // 5 bins from 0-0.05

    result.histogram.forEach(bin => {
      expect(bin.range).toBeDefined();
      expect(bin.count).toBeGreaterThanOrEqual(0);
      expect(bin.expected).toBeGreaterThan(0);
    });
  });

  it('assesses evidential value', () => {
    const result = pCurveAnalysis(studiesWithP);

    expect(result.evidentialValue).toBeDefined();
    expect(['present', 'absent', 'inconclusive']).toContain(result.evidentialValue);
  });

  it('calculates p-values from effect sizes when not provided', () => {
    const result = pCurveAnalysis(studiesES);

    expect(result.error).toBeUndefined();
    expect(result.k).toBe(6);
    expect(result.pValues.length).toBe(6);
  });

  it('detects potential p-hacking', () => {
    const result = pCurveAnalysis(pHackedStudies);

    // P-hacked studies should show flat or left-skewed distribution
    // More values near 0.05 than near 0.01
    expect(result.tests.extremeRatio.near05).toBeGreaterThan(0);
  });

  it('includes only significant studies', () => {
    const mixedStudies = [
      ...studiesWithP,
      { id: 'NS1', es: 0.10, se: 0.15, p: 0.50, excluded: false },
      { id: 'NS2', es: 0.05, se: 0.20, p: 0.80, excluded: false }
    ];

    const result = pCurveAnalysis(mixedStudies);

    expect(result.k).toBe(7);  // Only the 7 significant ones
  });

  it('returns error for insufficient significant studies', () => {
    const tooFew = studiesWithP.slice(0, 2);
    const result = pCurveAnalysis(tooFew);

    expect(result.error).toBeDefined();
  });

  it('provides interpretation', () => {
    const result = pCurveAnalysis(studiesWithP);

    expect(result.interpretation).toBeDefined();
    expect(typeof result.interpretation).toBe('string');
  });
});

describe('pUniform', () => {
  it('estimates bias-corrected effect size', () => {
    const result = pUniform(studiesES);

    expect(result.error).toBeUndefined();
    expect(result.method).toBe('P-Uniform');
    expect(Number.isFinite(result.estimate.mu)).toBe(true);
  });

  it('provides standard error and CI', () => {
    const result = pUniform(studiesES);

    expect(result.estimate.se).toBeGreaterThan(0);
    expect(result.estimate.ci.lo).toBeLessThan(result.estimate.mu);
    expect(result.estimate.ci.hi).toBeGreaterThan(result.estimate.mu);
  });

  it('reports number of significant studies used', () => {
    const result = pUniform(studiesES);

    expect(result.k).toBeGreaterThanOrEqual(3);
    expect(result.kTotal).toBeGreaterThanOrEqual(result.k);
  });

  it('compares with naive estimate', () => {
    const result = pUniform(studiesES);

    expect(Number.isFinite(result.naiveEstimate)).toBe(true);
    expect(Number.isFinite(result.biasCorrectionAmount)).toBe(true);
  });

  it('detects publication bias', () => {
    const result = pUniform(studiesES);

    expect(typeof result.biasDetected).toBe('boolean');
  });

  it('provides test statistics', () => {
    const result = pUniform(studiesES);

    expect(Number.isFinite(result.estimate.z)).toBe(true);
    expect(result.estimate.p).toBeGreaterThanOrEqual(0);
    expect(result.estimate.p).toBeLessThanOrEqual(1);
  });

  it('returns error for insufficient studies', () => {
    const tooFew = studiesES.slice(0, 2);
    const result = pUniform(tooFew);

    expect(result.error).toBeDefined();
  });

  it('includes interpretation', () => {
    const result = pUniform(studiesES);

    expect(result.interpretation).toBeDefined();
  });
});

describe('pUniformStar', () => {
  it('estimates effect with heterogeneity', () => {
    const result = pUniformStar(studiesES);

    expect(result.error).toBeUndefined();
    expect(result.method).toBe('P-Uniform*');
    expect(Number.isFinite(result.estimate.mu)).toBe(true);
  });

  it('estimates tau² (heterogeneity)', () => {
    const result = pUniformStar(studiesES);

    expect(result.heterogeneity).toBeDefined();
    expect(result.heterogeneity.tau2).toBeGreaterThanOrEqual(0);
    expect(result.heterogeneity.tau).toBeGreaterThanOrEqual(0);
    expect(result.heterogeneity.I2).toBeGreaterThanOrEqual(0);
    expect(result.heterogeneity.I2).toBeLessThanOrEqual(100);
  });

  it('provides confidence interval', () => {
    const result = pUniformStar(studiesES);

    expect(result.estimate.ci.lo).toBeLessThan(result.estimate.mu);
    expect(result.estimate.ci.hi).toBeGreaterThan(result.estimate.mu);
  });

  it('compares with p-uniform estimate', () => {
    const result = pUniformStar(studiesES);

    expect(Number.isFinite(result.pUniformEstimate)).toBe(true);
    expect(Number.isFinite(result.naiveEstimate)).toBe(true);
  });

  it('requires more studies than p-uniform', () => {
    const threeStudies = studiesES.slice(0, 3);
    const result = pUniformStar(threeStudies);

    // p-uniform* requires at least 4 significant studies
    expect(result.error).toBeDefined();
  });

  it('includes interpretation', () => {
    const result = pUniformStar(studiesES);

    expect(result.interpretation).toBeDefined();
    expect(result.interpretation).toContain('τ²');
  });
});

describe('compareSelectionModels', () => {
  it('compares all selection model methods', () => {
    const result = compareSelectionModels(studiesES);

    expect(result.estimates).toBeDefined();
    expect(result.estimates.naive).toBeDefined();
  });

  it('includes p-curve results when valid', () => {
    const result = compareSelectionModels(studiesES);

    if (result.estimates.pCurve !== null) {
      expect(result.estimates.pCurve.evidentialValue).toBeDefined();
    }
  });

  it('includes p-uniform estimates when valid', () => {
    const result = compareSelectionModels(studiesES);

    if (result.estimates.pUniform !== null) {
      expect(Number.isFinite(result.estimates.pUniform)).toBe(true);
    }
  });

  it('assesses estimate concordance', () => {
    const result = compareSelectionModels(studiesES);

    expect(result.concordance).toBeDefined();
    expect(result.concordance.range).toBeGreaterThanOrEqual(0);
    expect(typeof result.concordance.concordant).toBe('boolean');
    expect(result.concordance.interpretation).toBeDefined();
  });

  it('provides recommendation', () => {
    const result = compareSelectionModels(studiesES);

    expect(result.recommendation).toBeDefined();
    expect(typeof result.recommendation).toBe('string');
  });
});
