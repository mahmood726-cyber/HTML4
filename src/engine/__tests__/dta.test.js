/**
 * Tests for Diagnostic Test Accuracy (DTA) Meta-Analysis
 */
import { describe, it, expect } from 'vitest';
import {
  calculateAccuracyMeasures,
  bivariateModel,
  generateSROCCurve,
  hsrocModel,
  dtaMetaRegression
} from '../dta.js';

// Test DTA studies with 2x2 data
const dtaStudies = [
  { id: 'S1', tp: 80, fp: 10, fn: 20, tn: 90 },
  { id: 'S2', tp: 70, fp: 15, fn: 30, tn: 85 },
  { id: 'S3', tp: 85, fp: 8, fn: 15, tn: 92 },
  { id: 'S4', tp: 75, fp: 12, fn: 25, tn: 88 },
  { id: 'S5', tp: 78, fp: 11, fn: 22, tn: 89 },
];

// Studies with covariate
const dtaWithCovariate = dtaStudies.map((s, i) => ({
  ...s,
  year: 2000 + i * 4
}));

describe('Calculate Accuracy Measures', () => {
  it('should calculate sensitivity and specificity', () => {
    const result = calculateAccuracyMeasures({ tp: 80, fp: 10, fn: 20, tn: 90 });

    expect(result.error).toBeUndefined();
    expect(result.sensitivity.estimate).toBeCloseTo(0.8, 2);  // 80/100
    expect(result.specificity.estimate).toBeCloseTo(0.9, 2);  // 90/100
  });

  it('should calculate confidence intervals', () => {
    const result = calculateAccuracyMeasures({ tp: 80, fp: 10, fn: 20, tn: 90 });

    expect(result.sensitivity.ci.lo).toBeLessThan(0.8);
    expect(result.sensitivity.ci.hi).toBeGreaterThan(0.8);
    expect(result.specificity.ci.lo).toBeLessThan(0.9);
    expect(result.specificity.ci.hi).toBeGreaterThan(0.9);
  });

  it('should calculate likelihood ratios', () => {
    const result = calculateAccuracyMeasures({ tp: 80, fp: 10, fn: 20, tn: 90 });

    // LR+ = sens / (1 - spec) = 0.8 / 0.1 = 8
    expect(result.lrPlus.estimate).toBeCloseTo(8, 1);
    // LR- = (1 - sens) / spec = 0.2 / 0.9 ≈ 0.22
    expect(result.lrMinus.estimate).toBeCloseTo(0.22, 1);
  });

  it('should calculate diagnostic odds ratio', () => {
    const result = calculateAccuracyMeasures({ tp: 80, fp: 10, fn: 20, tn: 90 });

    // DOR = (TP * TN) / (FP * FN) = (80*90)/(10*20) = 36
    expect(result.dor.estimate).toBeCloseTo(36, 0);
  });

  it('should calculate predictive values', () => {
    const result = calculateAccuracyMeasures({ tp: 80, fp: 10, fn: 20, tn: 90 });

    // PPV = TP / (TP + FP) = 80/90 ≈ 0.89
    expect(result.ppv.estimate).toBeCloseTo(0.889, 2);
    // NPV = TN / (TN + FN) = 90/110 ≈ 0.82
    expect(result.npv.estimate).toBeCloseTo(0.818, 2);
  });

  it('should handle zero cells gracefully', () => {
    const result = calculateAccuracyMeasures({ tp: 50, fp: 0, fn: 10, tn: 40 });

    expect(result.sensitivity.estimate).toBeCloseTo(0.833, 2);
    expect(result.specificity.estimate).toBe(1);  // No FP
  });
});

describe('Bivariate Model', () => {
  it('should fit bivariate random effects model', () => {
    const result = bivariateModel(dtaStudies);

    expect(result.error).toBeUndefined();
    expect(result.method).toBe('Bivariate Random Effects');
    expect(result.k).toBe(5);
  });

  it('should estimate summary sensitivity and specificity', () => {
    const result = bivariateModel(dtaStudies);

    expect(result.sensitivity.estimate).toBeGreaterThan(0.7);
    expect(result.sensitivity.estimate).toBeLessThan(0.9);
    expect(result.specificity.estimate).toBeGreaterThan(0.8);
    expect(result.specificity.estimate).toBeLessThan(1);
  });

  it('should provide confidence intervals', () => {
    const result = bivariateModel(dtaStudies);

    expect(result.sensitivity.ci.lo).toBeLessThan(result.sensitivity.estimate);
    expect(result.sensitivity.ci.hi).toBeGreaterThan(result.sensitivity.estimate);
  });

  it('should calculate heterogeneity', () => {
    const result = bivariateModel(dtaStudies);

    expect(result.heterogeneity.tau2Sens).toBeGreaterThanOrEqual(0);
    expect(result.heterogeneity.tau2Spec).toBeGreaterThanOrEqual(0);
    expect(result.heterogeneity.correlation).toBeDefined();
  });

  it('should calculate summary DOR', () => {
    const result = bivariateModel(dtaStudies);

    expect(result.summaryDOR).toBeGreaterThan(1);  // Test better than chance
  });

  it('should require minimum studies', () => {
    const result = bivariateModel(dtaStudies.slice(0, 2));
    expect(result.error).toBeDefined();
  });
});

describe('SROC Curve', () => {
  it('should generate SROC curve data', () => {
    const bivar = bivariateModel(dtaStudies);
    const sroc = generateSROCCurve(bivar);

    expect(sroc.error).toBeUndefined();
    expect(sroc.curve.length).toBeGreaterThan(0);
    expect(sroc.summaryPoint).toBeDefined();
  });

  it('should include study points', () => {
    const bivar = bivariateModel(dtaStudies);
    const sroc = generateSROCCurve(bivar);

    expect(sroc.studies.length).toBe(5);
    sroc.studies.forEach(s => {
      expect(s.fpr).toBeGreaterThanOrEqual(0);
      expect(s.fpr).toBeLessThanOrEqual(1);
      expect(s.sens).toBeGreaterThanOrEqual(0);
      expect(s.sens).toBeLessThanOrEqual(1);
    });
  });

  it('should calculate AUC', () => {
    const bivar = bivariateModel(dtaStudies);
    const sroc = generateSROCCurve(bivar);

    expect(sroc.auc).toBeGreaterThan(0.5);  // Better than chance
    expect(sroc.auc).toBeLessThanOrEqual(1);
  });

  it('should include confidence region', () => {
    const bivar = bivariateModel(dtaStudies);
    const sroc = generateSROCCurve(bivar);

    expect(sroc.confRegion.length).toBeGreaterThan(0);
    sroc.confRegion.forEach(point => {
      expect(point.sensLo).toBeLessThanOrEqual(point.sensHi);
    });
  });
});

describe('HSROC Model', () => {
  it('should fit HSROC model', () => {
    const result = hsrocModel(dtaStudies);

    expect(result.error).toBeUndefined();
    expect(result.method).toBe('HSROC');
    expect(result.k).toBe(5);
  });

  it('should estimate HSROC parameters', () => {
    const result = hsrocModel(dtaStudies);

    expect(result.parameters.Lambda).toBeDefined();  // Accuracy
    expect(result.parameters.ThetaStar).toBeDefined();  // Threshold
    expect(result.parameters.beta).toBeDefined();  // Asymmetry
  });

  it('should provide summary estimates', () => {
    const result = hsrocModel(dtaStudies);

    expect(result.sensitivity.estimate).toBeGreaterThan(0);
    expect(result.sensitivity.estimate).toBeLessThanOrEqual(1);
    expect(result.specificity.estimate).toBeGreaterThan(0);
    expect(result.specificity.estimate).toBeLessThanOrEqual(1);
  });

  it('should interpret curve symmetry', () => {
    const result = hsrocModel(dtaStudies);

    expect(result.interpretation.symmetry).toBeDefined();
    expect(result.interpretation.accuracy).toBeDefined();
  });
});

describe('DTA Meta-Regression', () => {
  it('should perform covariate analysis', () => {
    const result = dtaMetaRegression(dtaWithCovariate, 'year');

    expect(result.error).toBeUndefined();
    expect(result.method).toBe('DTA Meta-Regression');
    expect(result.covariate).toBe('year');
  });

  it('should estimate covariate effects on sens and spec', () => {
    const result = dtaMetaRegression(dtaWithCovariate, 'year');

    expect(result.sensitivity.slope).toBeDefined();
    expect(result.specificity.slope).toBeDefined();
    expect(result.sensitivity.R2).toBeDefined();
    expect(result.specificity.R2).toBeDefined();
  });

  it('should provide interpretation', () => {
    const result = dtaMetaRegression(dtaWithCovariate, 'year');

    expect(result.interpretation).toBeDefined();
    expect(typeof result.interpretation).toBe('string');
  });

  it('should require minimum studies', () => {
    const result = dtaMetaRegression(dtaWithCovariate.slice(0, 3), 'year');
    expect(result.error).toBeDefined();
  });
});
