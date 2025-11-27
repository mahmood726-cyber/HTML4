/**
 * Unit tests for Robust Variance Estimation (RVE) module
 * Tests cluster-robust variance estimation for dependent effect sizes
 */

import { describe, it, expect } from 'vitest';
import {
  constructVMatrix,
  rveMetaAnalysis,
  rveMetaRegression,
  rveSensitivityAnalysis
} from '../rve.js';

// Test dataset with clustered effects (multiple outcomes per study)
const clusteredStudies = [
  { id: 'Study1_a', clusterId: 'Study1', es: 0.3, vi: 0.04, se: 0.2, excluded: false },
  { id: 'Study1_b', clusterId: 'Study1', es: 0.35, vi: 0.05, se: 0.224, excluded: false },
  { id: 'Study2_a', clusterId: 'Study2', es: 0.2, vi: 0.03, se: 0.173, excluded: false },
  { id: 'Study2_b', clusterId: 'Study2', es: 0.25, vi: 0.04, se: 0.2, excluded: false },
  { id: 'Study2_c', clusterId: 'Study2', es: 0.22, vi: 0.035, se: 0.187, excluded: false },
  { id: 'Study3_a', clusterId: 'Study3', es: 0.4, vi: 0.06, se: 0.245, excluded: false },
  { id: 'Study4_a', clusterId: 'Study4', es: 0.15, vi: 0.03, se: 0.173, excluded: false },
  { id: 'Study4_b', clusterId: 'Study4', es: 0.18, vi: 0.032, se: 0.179, excluded: false }
];

describe('constructVMatrix', () => {
  it('creates variance-covariance matrix with correct dimensions', () => {
    const result = constructVMatrix(clusteredStudies);

    expect(result.V.length).toBe(8);
    expect(result.V[0].length).toBe(8);
    expect(result.n).toBe(8);
    expect(result.m).toBe(4);
  });

  it('identifies correct number of clusters', () => {
    const result = constructVMatrix(clusteredStudies);

    expect(result.clusters.length).toBe(4);
    expect(result.clusters[0].id).toBe('Study1');
    expect(result.clusters[0].size).toBe(2);
    expect(result.clusters[1].size).toBe(3);
  });

  it('uses diagonal variance for same-effect entries', () => {
    const result = constructVMatrix(clusteredStudies);

    // First effect variance should be on diagonal
    expect(result.V[0][0]).toBeCloseTo(0.04, 4);
  });

  it('uses correlation for within-cluster off-diagonal entries', () => {
    const result = constructVMatrix(clusteredStudies, 0.8);

    // Off-diagonal within Study1 cluster
    const expectedCov = 0.8 * Math.sqrt(0.04 * 0.05);
    expect(result.V[0][1]).toBeCloseTo(expectedCov, 4);
    expect(result.V[1][0]).toBeCloseTo(expectedCov, 4);
  });

  it('uses zero for between-cluster entries', () => {
    const result = constructVMatrix(clusteredStudies);

    // Between Study1 (indices 0-1) and Study2 (indices 2-4)
    expect(result.V[0][2]).toBe(0);
    expect(result.V[1][3]).toBe(0);
  });

  it('respects different rho values', () => {
    const result05 = constructVMatrix(clusteredStudies, 0.5);
    const result08 = constructVMatrix(clusteredStudies, 0.8);

    // Higher rho should give higher off-diagonal values
    expect(result08.V[0][1]).toBeGreaterThan(result05.V[0][1]);
  });
});

describe('rveMetaAnalysis', () => {
  it('calculates pooled effect estimate', () => {
    const result = rveMetaAnalysis(clusteredStudies);

    expect(result.error).toBeUndefined();
    expect(result.es).toBeGreaterThan(0);
    expect(result.es).toBeLessThan(1);
    expect(Number.isFinite(result.es)).toBe(true);
  });

  it('returns correct effect and cluster counts', () => {
    const result = rveMetaAnalysis(clusteredStudies);

    expect(result.nEffects).toBe(8);
    expect(result.nClusters).toBe(4);
  });

  it('calculates robust standard error', () => {
    const result = rveMetaAnalysis(clusteredStudies);

    expect(result.se).toBeGreaterThan(0);
    expect(Number.isFinite(result.se)).toBe(true);
  });

  it('provides confidence interval', () => {
    const result = rveMetaAnalysis(clusteredStudies);

    expect(result.ci.lo).toBeLessThan(result.es);
    expect(result.ci.hi).toBeGreaterThan(result.es);
  });

  it('calculates Satterthwaite degrees of freedom', () => {
    const result = rveMetaAnalysis(clusteredStudies);

    expect(result.df).toBeGreaterThan(0);
    expect(result.df).toBeLessThanOrEqual(result.nClusters);
  });

  it('supports different small-sample corrections', () => {
    const cr0 = rveMetaAnalysis(clusteredStudies, { smallSampleCorrection: 'CR0' });
    const cr1 = rveMetaAnalysis(clusteredStudies, { smallSampleCorrection: 'CR1' });
    const cr2 = rveMetaAnalysis(clusteredStudies, { smallSampleCorrection: 'CR2' });

    expect(cr0.correction).toBe('CR0');
    expect(cr1.correction).toBe('CR1');
    expect(cr2.correction).toBe('CR2');

    // CR2 typically gives larger SE than CR0
    expect(cr2.se).toBeGreaterThanOrEqual(cr0.se * 0.9);
  });

  it('returns error for insufficient clusters', () => {
    const singleCluster = clusteredStudies.filter(s => s.clusterId === 'Study1');
    const result = rveMetaAnalysis(singleCluster);

    expect(result.error).toBeDefined();
  });

  it('provides heterogeneity approximation', () => {
    const result = rveMetaAnalysis(clusteredStudies);

    expect(result.heterogeneity.tau2Approx).toBeGreaterThanOrEqual(0);
    expect(result.heterogeneity.I2Approx).toBeGreaterThanOrEqual(0);
    expect(result.heterogeneity.I2Approx).toBeLessThanOrEqual(100);
  });

  it('includes interpretation text', () => {
    const result = rveMetaAnalysis(clusteredStudies);

    expect(result.interpretation).toBeDefined();
    expect(typeof result.interpretation).toBe('string');
    expect(result.interpretation.length).toBeGreaterThan(0);
  });
});

describe('rveMetaRegression', () => {
  const studiesWithCovariates = clusteredStudies.map((s, i) => ({
    ...s,
    year: 2010 + i,
    quality: Math.random() * 10
  }));

  it('estimates regression coefficients', () => {
    const result = rveMetaRegression(studiesWithCovariates, ['year']);

    expect(result.error).toBeUndefined();
    expect(result.coefficients.length).toBe(2); // Intercept + year
    expect(result.coefficients[0].name).toBe('(Intercept)');
    expect(result.coefficients[1].name).toBe('year');
  });

  it('provides standard errors for each coefficient', () => {
    const result = rveMetaRegression(studiesWithCovariates, ['year']);

    result.coefficients.forEach(coef => {
      expect(coef.se).toBeGreaterThan(0);
      expect(Number.isFinite(coef.se)).toBe(true);
    });
  });

  it('calculates confidence intervals', () => {
    const result = rveMetaRegression(studiesWithCovariates, ['year']);

    result.coefficients.forEach(coef => {
      expect(coef.ci.lo).toBeLessThan(coef.estimate);
      expect(coef.ci.hi).toBeGreaterThan(coef.estimate);
    });
  });

  it('calculates R² for model fit', () => {
    const result = rveMetaRegression(studiesWithCovariates, ['year']);

    expect(result.R2).toBeGreaterThanOrEqual(0);
    expect(result.R2).toBeLessThanOrEqual(100);
  });

  it('returns error for insufficient studies', () => {
    const fewStudies = studiesWithCovariates.slice(0, 2);
    const result = rveMetaRegression(fewStudies, ['year', 'quality']);

    expect(result.error).toBeDefined();
  });
});

describe('rveSensitivityAnalysis', () => {
  it('analyzes sensitivity to rho values', () => {
    const result = rveSensitivityAnalysis(clusteredStudies);

    expect(result.results.length).toBeGreaterThan(0);
    expect(result.results[0].rho).toBeDefined();
    expect(result.results[0].es).toBeDefined();
  });

  it('tests multiple rho values', () => {
    const result = rveSensitivityAnalysis(clusteredStudies, {
      rhoValues: [0, 0.5, 1.0]
    });

    expect(result.results.length).toBe(3);
    expect(result.results[0].rho).toBe(0);
    expect(result.results[2].rho).toBe(1.0);
  });

  it('calculates estimate range', () => {
    const result = rveSensitivityAnalysis(clusteredStudies);

    expect(result.sensitivity.estimateRange).toBeGreaterThanOrEqual(0);
    expect(result.sensitivity.minES).toBeLessThanOrEqual(result.sensitivity.maxES);
  });

  it('assesses conclusion robustness', () => {
    const result = rveSensitivityAnalysis(clusteredStudies);

    expect(typeof result.sensitivity.conclusionRobust).toBe('boolean');
    expect(result.interpretation).toBeDefined();
  });
});
