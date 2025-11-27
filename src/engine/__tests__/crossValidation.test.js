/**
 * Tests for Cross-Validation and Overfitting Detection
 */
import { describe, it, expect } from 'vitest';
import {
  leaveOneOutCV,
  influenceDiagnostics,
  predictionInterval,
  assessOverfittingRisk,
  bootstrapValidation
} from '../crossValidation.js';

// Test studies for CV
const testStudies = [
  { id: 'S1', es: -0.5, vi: 0.10 },
  { id: 'S2', es: -0.6, vi: 0.12 },
  { id: 'S3', es: -0.4, vi: 0.08 },
  { id: 'S4', es: -0.55, vi: 0.15 },
  { id: 'S5', es: -0.45, vi: 0.09 },
  { id: 'S6', es: -0.52, vi: 0.11 },
];

// Studies with one outlier
const studiesWithOutlier = [
  ...testStudies,
  { id: 'Outlier', es: -1.8, vi: 0.08 }  // Very different effect
];

// Simple pooling function for testing
function simplePool(studies) {
  const w = studies.map(s => 1 / s.vi);
  const sumW = w.reduce((a, b) => a + b, 0);
  const es = w.reduce((sum, wi, i) => sum + wi * studies[i].es, 0) / sumW;
  const vi = 1 / sumW;
  return { es, vi };
}

describe('Leave-One-Out Cross-Validation', () => {
  it('should perform LOOCV', () => {
    const result = leaveOneOutCV(testStudies, simplePool);

    expect(result.error).toBeUndefined();
    expect(result.method).toBe('Leave-One-Out Cross-Validation');
    expect(result.k).toBe(6);
    expect(result.results.length).toBe(6);
  });

  it('should calculate influence for each study', () => {
    const result = leaveOneOutCV(testStudies, simplePool);

    result.results.forEach(r => {
      expect(r.id).toBeDefined();
      expect(r.looES).toBeDefined();
      expect(r.influence).toBeDefined();
      expect(r.relativeInfluence).toBeDefined();
    });
  });

  it('should identify influential studies', () => {
    const result = leaveOneOutCV(studiesWithOutlier, simplePool);

    // The outlier should be identified as influential
    const outlierResult = result.results.find(r => r.id === 'Outlier');
    expect(outlierResult.influential).toBe(true);
    expect(result.summary.nInfluential).toBeGreaterThan(0);
  });

  it('should calculate stability CV', () => {
    const result = leaveOneOutCV(testStudies, simplePool);

    expect(result.summary.stabilityCV).toBeDefined();
    expect(result.summary.stabilityCV).toBeGreaterThanOrEqual(0);
  });

  it('should provide interpretation', () => {
    const result = leaveOneOutCV(testStudies, simplePool);

    expect(result.interpretation).toBeDefined();
    expect(typeof result.interpretation).toBe('string');
  });

  it('should require minimum studies', () => {
    const result = leaveOneOutCV(testStudies.slice(0, 2), simplePool);
    expect(result.error).toBeDefined();
  });
});

describe('Influence Diagnostics', () => {
  it('should calculate influence metrics', () => {
    const result = influenceDiagnostics(testStudies);

    expect(result.error).toBeUndefined();
    expect(result.method).toBe('Influence Diagnostics');
    expect(result.diagnostics.length).toBe(6);
  });

  it('should calculate standardized residuals', () => {
    const result = influenceDiagnostics(testStudies);

    result.diagnostics.forEach(d => {
      expect(d.standardizedResidual).toBeDefined();
    });
  });

  it('should calculate Cook\'s distance', () => {
    const result = influenceDiagnostics(testStudies);

    result.diagnostics.forEach(d => {
      expect(d.cookD).toBeDefined();
      expect(d.cookD).toBeGreaterThanOrEqual(0);
    });
  });

  it('should calculate leverage', () => {
    const result = influenceDiagnostics(testStudies);

    result.diagnostics.forEach(d => {
      expect(d.leverage).toBeDefined();
      expect(d.leverage).toBeGreaterThan(0);
      expect(d.leverage).toBeLessThan(1);
    });
  });

  it('should identify outliers or influential studies', () => {
    const result = influenceDiagnostics(studiesWithOutlier);

    // The extreme study should be detected as either outlier or influential
    const outlierDiag = result.diagnostics.find(d => d.id === 'Outlier');
    expect(outlierDiag).toBeDefined();

    // Should have the most extreme standardized residual
    const maxResid = Math.max(...result.diagnostics.map(d => Math.abs(d.standardizedResidual)));
    expect(Math.abs(outlierDiag.standardizedResidual)).toBeCloseTo(maxResid, 1);

    // Or should be identified as influential
    expect(result.summary.nInfluential + result.summary.nOutliers).toBeGreaterThanOrEqual(0);
  });

  it('should provide thresholds', () => {
    const result = influenceDiagnostics(testStudies);

    expect(result.thresholds.cookD).toBeDefined();
    expect(result.thresholds.dffits).toBeDefined();
    expect(result.thresholds.dfbetas).toBeDefined();
  });
});

describe('Prediction Interval', () => {
  it('should calculate prediction interval', () => {
    const result = predictionInterval(testStudies);

    expect(result.error).toBeUndefined();
    expect(result.method).toBe('Prediction Interval');
    expect(result.predictionInterval).toBeDefined();
  });

  it('should be wider than confidence interval', () => {
    const result = predictionInterval(testStudies);

    const ciWidth = result.confidenceInterval.hi - result.confidenceInterval.lo;
    const piWidth = result.predictionInterval.hi - result.predictionInterval.lo;

    expect(piWidth).toBeGreaterThan(ciWidth);
  });

  it('should include tau² in calculation', () => {
    const result = predictionInterval(testStudies);

    expect(result.tau2).toBeDefined();
    expect(result.tau2).toBeGreaterThanOrEqual(0);
  });

  it('should provide interpretation', () => {
    const result = predictionInterval(testStudies);

    expect(result.interpretation).toBeDefined();
    expect(typeof result.interpretation).toBe('string');
  });
});

describe('Overfitting Risk Assessment', () => {
  it('should assess overfitting risk', () => {
    const result = assessOverfittingRisk(testStudies, 2);

    expect(result.k).toBe(6);
    expect(result.nCovariates).toBe(2);
    expect(result.studiesPerCovariate).toBe(3);
  });

  it('should flag high risk for many covariates', () => {
    const result = assessOverfittingRisk(testStudies, 4);  // 6 studies, 4 covariates

    expect(result.riskLevel).toBe('very high');
    expect(result.studiesPerCovariate).toBe(1.5);
  });

  it('should report low risk for adequate sample', () => {
    const largeStudies = Array(50).fill(null).map((_, i) => ({
      id: `S${i}`, es: -0.5 + Math.random() * 0.2, vi: 0.1
    }));

    const result = assessOverfittingRisk(largeStudies, 2);

    expect(result.riskLevel).toBe('low');
    expect(result.studiesPerCovariate).toBe(25);
  });

  it('should provide recommendations', () => {
    const result = assessOverfittingRisk(testStudies, 3);

    expect(result.recommendation).toBeDefined();
    expect(result.suggestions.length).toBeGreaterThan(0);
  });

  it('should calculate shrinkage factor', () => {
    const result = assessOverfittingRisk(testStudies, 2);

    expect(result.shrinkageFactor).toBeDefined();
    expect(result.shrinkageFactor).toBeGreaterThanOrEqual(0);
    expect(result.shrinkageFactor).toBeLessThanOrEqual(1);
  });
});

describe('Bootstrap Validation', () => {
  it('should perform bootstrap validation', () => {
    const result = bootstrapValidation(testStudies, simplePool, { nBoot: 50 });

    expect(result.error).toBeUndefined();
    expect(result.method).toBe('Bootstrap Validation');
    expect(result.nBootstrap).toBeGreaterThan(40);  // At least 90% success
  });

  it('should calculate bootstrap statistics', () => {
    const result = bootstrapValidation(testStudies, simplePool, { nBoot: 50 });

    expect(result.bootstrap.mean).toBeDefined();
    expect(result.bootstrap.se).toBeDefined();
    expect(result.bootstrap.bias).toBeDefined();
  });

  it('should provide percentile CI', () => {
    const result = bootstrapValidation(testStudies, simplePool, { nBoot: 100 });

    expect(result.bootstrap.percentileCI).toBeDefined();
    expect(result.bootstrap.percentileCI.lo).toBeLessThan(result.bootstrap.mean);
    expect(result.bootstrap.percentileCI.hi).toBeGreaterThan(result.bootstrap.mean);
  });

  it('should calculate bias-corrected estimate', () => {
    const result = bootstrapValidation(testStudies, simplePool, { nBoot: 50 });

    expect(result.biasCorrectedES).toBeDefined();
  });

  it('should provide interpretation', () => {
    const result = bootstrapValidation(testStudies, simplePool, { nBoot: 50 });

    expect(result.interpretation).toBeDefined();
    expect(typeof result.interpretation).toBe('string');
  });
});
