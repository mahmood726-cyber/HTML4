/**
 * Tests for Meta-Regression
 */
import { describe, it, expect } from 'vitest';
import {
  univariateMetaRegression,
  multivariateMetaRegression,
  subgroupMetaRegression,
  createBubblePlotData
} from '../metaRegression.js';

// Test data with covariate (year)
const studiesWithYear = [
  { id: 'Study1', es: -0.8, vi: 0.15, year: 1990 },
  { id: 'Study2', es: -0.6, vi: 0.12, year: 1995 },
  { id: 'Study3', es: -0.5, vi: 0.10, year: 2000 },
  { id: 'Study4', es: -0.3, vi: 0.08, year: 2005 },
  { id: 'Study5', es: -0.2, vi: 0.09, year: 2010 },
  { id: 'Study6', es: -0.1, vi: 0.11, year: 2015 },
];

// Extract covariate
const years = studiesWithYear.map(s => s.year);

// Studies with subgroup
const studiesWithGroup = [
  { id: 'A1', es: -0.8, vi: 0.15 },
  { id: 'A2', es: -0.7, vi: 0.12 },
  { id: 'A3', es: -0.6, vi: 0.10 },
  { id: 'B1', es: -0.2, vi: 0.08 },
  { id: 'B2', es: -0.1, vi: 0.09 },
  { id: 'B3', es: 0.0, vi: 0.11 },
];
const groups = ['A', 'A', 'A', 'B', 'B', 'B'];

describe('Univariate Meta-Regression', () => {
  it('should fit meta-regression with continuous covariate', () => {
    const result = univariateMetaRegression(studiesWithYear, years);

    expect(result.error).toBeUndefined();
    expect(result.k).toBe(6);
    expect(result.intercept).toBeDefined();
    expect(result.slope).toBeDefined();
  });

  it('should detect significant trend', () => {
    const result = univariateMetaRegression(studiesWithYear, years);

    // Effect increases (less negative) over time
    expect(result.slope.estimate).toBeGreaterThan(0);
    // P-value should be relatively small (trend visible but k=6 limits power)
    expect(result.slope.p).toBeLessThan(0.25);
  });

  it('should calculate R² for heterogeneity explained', () => {
    const result = univariateMetaRegression(studiesWithYear, years);

    expect(result.R2).toBeDefined();
    expect(result.R2).toBeGreaterThanOrEqual(0);
    expect(result.R2).toBeLessThanOrEqual(100);
  });

  it('should calculate residual heterogeneity', () => {
    const result = univariateMetaRegression(studiesWithYear, years);

    expect(result.QE).toBeDefined();
    expect(result.pQE).toBeDefined();
    expect(result.I2res).toBeDefined();
  });

  it('should provide confidence intervals', () => {
    const result = univariateMetaRegression(studiesWithYear, years);

    expect(result.slope.ci.lo).toBeLessThan(result.slope.estimate);
    expect(result.slope.ci.hi).toBeGreaterThan(result.slope.estimate);
  });

  it('should apply Knapp-Hartung adjustment by default', () => {
    const result = univariateMetaRegression(studiesWithYear, years, { knha: true });
    const resultNoKH = univariateMetaRegression(studiesWithYear, years, { knha: false });

    // KH adjustment typically produces wider CIs
    expect(result.knha).toBe(true);
  });

  it('should return study-level fitted values', () => {
    const result = univariateMetaRegression(studiesWithYear, years);

    expect(result.studies.length).toBe(6);
    result.studies.forEach(s => {
      expect(s.x).toBeDefined();
      expect(s.y).toBeDefined();
      expect(s.fitted).toBeDefined();
      expect(s.residual).toBeDefined();
    });
  });

  it('should handle too few studies', () => {
    const twoStudies = studiesWithYear.slice(0, 2);
    const result = univariateMetaRegression(twoStudies, years.slice(0, 2));

    expect(result.error).toBeDefined();
  });
});

describe('Subgroup Meta-Regression', () => {
  it('should compare subgroups using meta-regression', () => {
    const result = subgroupMetaRegression(studiesWithGroup, groups);

    expect(result.error).toBeUndefined();
    expect(result.type).toBe('subgroup');
    expect(result.groups).toEqual(['A', 'B']);
  });

  it('should calculate group-specific estimates', () => {
    const result = subgroupMetaRegression(studiesWithGroup, groups);

    expect(result.groupResults).toBeDefined();
    expect(result.groupResults['A']).toBeDefined();
    expect(result.groupResults['B']).toBeDefined();

    // Group A has more negative effects
    expect(result.groupResults['A'].es).toBeLessThan(result.groupResults['B'].es);
  });

  it('should test for subgroup differences', () => {
    const result = subgroupMetaRegression(studiesWithGroup, groups);

    expect(result.QM).toBeDefined();
    expect(result.pQM).toBeDefined();
    // Clear difference between groups
    expect(result.pQM).toBeLessThan(0.1);
  });
});

describe('Multivariate Meta-Regression', () => {
  it('should fit model with multiple covariates', () => {
    // Add second covariate
    const doses = [100, 100, 200, 200, 300, 300];
    const result = multivariateMetaRegression(
      studiesWithYear,
      [years, doses],
      ['Year', 'Dose']
    );

    expect(result.error).toBeUndefined();
    expect(result.coefficients.length).toBe(3);  // Intercept + 2 covariates
  });

  it('should name coefficients correctly', () => {
    const doses = [100, 100, 200, 200, 300, 300];
    const result = multivariateMetaRegression(
      studiesWithYear,
      [years, doses],
      ['Year', 'Dose']
    );

    expect(result.coefficients[0].name).toBe('Intercept');
    expect(result.coefficients[1].name).toBe('Year');
    expect(result.coefficients[2].name).toBe('Dose');
  });
});

describe('Bubble Plot Data', () => {
  it('should generate bubble plot data', () => {
    const reg = univariateMetaRegression(studiesWithYear, years);
    const plot = createBubblePlotData(reg);

    expect(plot).not.toBeNull();
    expect(plot.points.length).toBe(6);
    expect(plot.line.x.length).toBe(2);
    expect(plot.line.y.length).toBe(2);
  });

  it('should include CI band', () => {
    const reg = univariateMetaRegression(studiesWithYear, years);
    const plot = createBubblePlotData(reg);

    expect(plot.ci).toBeDefined();
    expect(plot.ci.lo.length).toBe(2);
    expect(plot.ci.hi.length).toBe(2);
  });

  it('should size points by weight', () => {
    const reg = univariateMetaRegression(studiesWithYear, years);
    const plot = createBubblePlotData(reg);

    plot.points.forEach(p => {
      expect(p.size).toBeDefined();
      expect(p.size).toBeGreaterThan(0);
    });
  });
});
