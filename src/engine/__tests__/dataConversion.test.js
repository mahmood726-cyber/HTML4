/**
 * Tests for Data Conversion Utilities
 */
import { describe, it, expect } from 'vitest';
import {
  medianIQRToMeanSD,
  medianRangeToMeanSD,
  combinedToMeanSD,
  oddsRatioToRiskRatio,
  riskRatioToOddsRatio,
  cohensD_to_HedgesG,
  correlationToD,
  dToCorrelation,
  fisherZ,
  fisherZWithSE,
  inverseFisherZ,
  tToEffectSize,
  fToEffectSize,
  standardizeStudyData
} from '../dataConversion.js';

describe('Median/IQR to Mean/SD Conversion', () => {
  it('should convert median and IQR to mean and SD', () => {
    // Normal distribution: median ≈ mean
    const result = medianIQRToMeanSD(50, 40, 60, 100);

    expect(result.error).toBeUndefined();
    expect(result.method).toBe('Wan et al. (2014)');
    expect(result.mean).toBeCloseTo(50, 0);  // Mean ≈ median for symmetric
    expect(result.sd).toBeGreaterThan(0);
  });

  it('should estimate SD from IQR correctly', () => {
    // For normal: IQR ≈ 1.35 * SD, so SD ≈ IQR / 1.35
    const result = medianIQRToMeanSD(100, 86, 114, 200);  // IQR = 28

    expect(result.sd).toBeCloseTo(28 / 1.35, 1);  // ≈ 20.7
  });

  it('should handle small samples', () => {
    const result = medianIQRToMeanSD(50, 45, 55, 10);

    expect(result.error).toBeUndefined();
    expect(result.mean).toBeDefined();
    expect(result.sd).toBeGreaterThan(0);
  });
});

describe('Median/Range to Mean/SD Conversion', () => {
  it('should convert median and range to mean and SD', () => {
    const result = medianRangeToMeanSD(50, 20, 80, 50);

    expect(result.error).toBeUndefined();
    expect(result.mean).toBeCloseTo(50, 0);
    expect(result.sd).toBeGreaterThan(0);
  });

  it('should give reasonable SD estimates', () => {
    // For large n, range ≈ 4*SD (roughly)
    const result = medianRangeToMeanSD(100, 60, 140, 100);  // Range = 80

    expect(result.sd).toBeGreaterThan(10);
    expect(result.sd).toBeLessThan(30);
  });
});

describe('Combined Quantile Conversion', () => {
  it('should use both IQR and range', () => {
    const result = combinedToMeanSD({
      median: 50,
      q1: 40,
      q3: 60,
      min: 20,
      max: 80,
      n: 100
    });

    expect(result.error).toBeUndefined();
    expect(result.method).toBe('Luo et al. (2018)');
    expect(result.mean).toBeCloseTo(50, 0);
  });
});

describe('Odds Ratio / Risk Ratio Conversion', () => {
  it('should convert OR to RR', () => {
    // With baseline risk 0.1, OR 2.0 → RR ≈ 1.82
    const result = oddsRatioToRiskRatio(2.0, 0.1);

    expect(result.error).toBeUndefined();
    expect(result.rr).toBeLessThan(result.or);  // RR < OR when risk > 0
    expect(result.rr).toBeGreaterThan(1);
  });

  it('should convert RR to OR', () => {
    const result = riskRatioToOddsRatio(1.5, 0.2);

    expect(result.error).toBeUndefined();
    expect(result.or).toBeGreaterThan(result.rr);
  });

  it('should be consistent in both directions', () => {
    const p0 = 0.15;
    const originalOR = 2.5;

    const rr = oddsRatioToRiskRatio(originalOR, p0).rr;
    const backToOR = riskRatioToOddsRatio(rr, p0).or;

    expect(backToOR).toBeCloseTo(originalOR, 2);
  });
});

describe('Cohen\'s d to Hedges\' g', () => {
  it('should apply small sample correction', () => {
    const result = cohensD_to_HedgesG(0.5, 20, 20);

    expect(result.g).toBeLessThan(result.d);  // g is smaller due to correction
    expect(result.J).toBeLessThan(1);  // Correction factor < 1
  });

  it('should be close to d for large samples', () => {
    const result = cohensD_to_HedgesG(0.5, 500, 500);

    expect(result.g).toBeCloseTo(result.d, 2);  // Very close for large n
    expect(result.J).toBeCloseTo(1, 2);
  });

  it('should calculate SE of g', () => {
    const result = cohensD_to_HedgesG(0.5, 30, 30);

    expect(result.seG).toBeGreaterThan(0);
  });
});

describe('Correlation to d Conversions', () => {
  it('should convert r to d', () => {
    const result = correlationToD(0.3, 0.5);

    expect(result.d).toBeGreaterThan(0);
    // For r = 0.3: d ≈ 0.63
    expect(result.d).toBeCloseTo(0.63, 1);
  });

  it('should convert d to r', () => {
    const result = dToCorrelation(0.5, 50, 50);

    expect(result.r).toBeGreaterThan(0);
    expect(result.r).toBeLessThan(1);
  });

  it('should be approximately reversible', () => {
    const r = 0.4;
    const d = correlationToD(r, 0.5).d;
    const rBack = dToCorrelation(d, 100, 100).r;

    expect(rBack).toBeCloseTo(r, 1);
  });
});

describe('Fisher\'s z Transformation', () => {
  it('should transform r to z', () => {
    const result = fisherZ(0.5);

    expect(result.z).toBeGreaterThan(result.r);  // z > r for positive r
    expect(result.z).toBeCloseTo(0.549, 2);  // Known value
  });

  it('should provide SE with sample size', () => {
    const result = fisherZWithSE(0.5, 100);

    expect(result.se).toBeCloseTo(1 / Math.sqrt(97), 3);
    expect(result.ci).toBeDefined();
  });

  it('should inverse correctly', () => {
    const r = 0.6;
    const z = fisherZ(r).z;
    const rBack = inverseFisherZ(z).r;

    expect(rBack).toBeCloseTo(r, 5);
  });

  it('should handle negative correlations', () => {
    const result = fisherZ(-0.4);

    expect(result.z).toBeLessThan(0);
  });
});

describe('t-statistic to Effect Size', () => {
  it('should convert t to d', () => {
    const result = tToEffectSize(2.5, 30, 30);

    expect(result.d).toBeGreaterThan(0);
    expect(result.g).toBeLessThan(result.d);  // Hedges correction
    expect(result.seD).toBeGreaterThan(0);
  });

  it('should give larger d for larger t', () => {
    const small = tToEffectSize(1.5, 30, 30);
    const large = tToEffectSize(3.0, 30, 30);

    expect(large.d).toBeGreaterThan(small.d);
  });
});

describe('F-statistic to Effect Size', () => {
  it('should convert F to d via t', () => {
    const result = fToEffectSize(6.25, 30, 30);  // F = 6.25 = 2.5²

    // Should match t = 2.5
    const tResult = tToEffectSize(2.5, 30, 30);
    expect(result.d).toBeCloseTo(tResult.d, 5);
  });
});

describe('Standardize Study Data', () => {
  it('should convert mean/SD to SMD', () => {
    const studies = [
      {
        id: 'S1',
        mean1: 105, sd1: 15, n1: 50,
        mean2: 100, sd2: 15, n2: 50
      }
    ];

    const result = standardizeStudyData(studies);

    expect(result[0].es).toBeCloseTo(0.333, 2);  // d = (105-100)/15
    expect(result[0].vi).toBeDefined();
    expect(result[0].conversionMethod).toBe('mean/SD to SMD');
  });

  it('should convert t-statistic to SMD', () => {
    const studies = [
      { id: 'S1', t: 2.0, n1: 25, n2: 25 }
    ];

    const result = standardizeStudyData(studies);

    expect(result[0].es).toBeDefined();
    expect(result[0].vi).toBeDefined();
    expect(result[0].conversionMethod).toBe('t-statistic to SMD');
  });

  it('should convert correlation to Fisher z', () => {
    const studies = [
      { id: 'S1', r: 0.3, n: 100 }
    ];

    const result = standardizeStudyData(studies);

    expect(result[0].es).toBeCloseTo(fisherZ(0.3).z, 3);
    expect(result[0].conversionMethod).toBe('r to Fisher z');
  });

  it('should pass through already standardized data', () => {
    const studies = [
      { id: 'S1', es: 0.5, vi: 0.04 }
    ];

    const result = standardizeStudyData(studies);

    expect(result[0].es).toBe(0.5);
    expect(result[0].vi).toBe(0.04);
    expect(result[0].conversionMethod).toBeUndefined();
  });
});
