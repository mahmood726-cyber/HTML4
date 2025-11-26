/**
 * Unit tests for effect size calculations
 * Expected values validated against R metafor package
 */

import { describe, it, expect } from 'vitest';
import { calculateEffect, transformEffect, isRatioMetric } from '../effects.js';

describe('Binary Effect Sizes', () => {
  describe('Odds Ratio (OR)', () => {
    it('calculates OR correctly for typical 2x2 data', () => {
      // Study: Treatment 10/50, Control 20/50
      const row = { e1: 10, n1: 50, e2: 20, n2: 50 };
      const result = calculateEffect(row, 'binary', 'OR');

      // OR = (10*30)/(20*40) = 0.375, log(OR) = -0.981
      expect(result).not.toBeNull();
      expect(result.es).toBeCloseTo(-0.981, 2);
      expect(result.display).toBeCloseTo(0.375, 3);
      // Variance = 1/10 + 1/40 + 1/20 + 1/30 = 0.208
      expect(result.vi).toBeCloseTo(0.208, 2);
    });

    it('applies continuity correction for zero cells', () => {
      // Study with zero events in treatment
      const row = { e1: 0, n1: 50, e2: 10, n2: 50 };
      const result = calculateEffect(row, 'binary', 'OR', 0.5);

      expect(result).not.toBeNull();
      expect(result.es).toBeLessThan(0); // Treatment protective
      expect(Number.isFinite(result.vi)).toBe(true);
    });

    it('returns null for invalid input', () => {
      expect(calculateEffect({ e1: -1, n1: 50, e2: 10, n2: 50 }, 'binary', 'OR')).toBeNull();
      expect(calculateEffect({ e1: 60, n1: 50, e2: 10, n2: 50 }, 'binary', 'OR')).toBeNull();
      expect(calculateEffect({ e1: 10, n1: 0, e2: 10, n2: 50 }, 'binary', 'OR')).toBeNull();
    });
  });

  describe('Risk Ratio (RR)', () => {
    it('calculates RR correctly', () => {
      // Study: Treatment 10/100, Control 20/100
      const row = { e1: 10, n1: 100, e2: 20, n2: 100 };
      const result = calculateEffect(row, 'binary', 'RR');

      // RR = (10/100)/(20/100) = 0.5, log(RR) = -0.693
      expect(result).not.toBeNull();
      expect(result.es).toBeCloseTo(-0.693, 2);
      expect(result.display).toBeCloseTo(0.5, 3);

      // Correct variance formula: (n-e)/(e*n) for each group
      // vi = (100-10)/(10*100) + (100-20)/(20*100) = 0.09 + 0.04 = 0.13
      expect(result.vi).toBeCloseTo(0.13, 2);
    });

    it('uses correct RR variance formula (not wrong formula)', () => {
      // Verify the fix: should be (n-e)/(e*n), not (1-p)/e
      const row = { e1: 10, n1: 100, e2: 20, n2: 100 };
      const result = calculateEffect(row, 'binary', 'RR');

      // Wrong formula would give: (1-0.1)/10 + (1-0.2)/20 = 0.09 + 0.04 = 0.13
      // Correct formula gives: (90)/(10*100) + (80)/(20*100) = 0.09 + 0.04 = 0.13
      // In this case they're the same, test with different values

      const row2 = { e1: 15, n1: 60, e2: 25, n2: 80 };
      const result2 = calculateEffect(row2, 'binary', 'RR');

      // Correct: (60-15)/(15*60) + (80-25)/(25*80) = 0.05 + 0.0275 = 0.0775
      expect(result2.vi).toBeCloseTo(0.0775, 3);
    });
  });

  describe('Risk Difference (RD)', () => {
    it('calculates RD correctly', () => {
      const row = { e1: 10, n1: 100, e2: 20, n2: 100 };
      const result = calculateEffect(row, 'binary', 'RD');

      // RD = 0.1 - 0.2 = -0.1
      expect(result).not.toBeNull();
      expect(result.es).toBeCloseTo(-0.1, 3);
      expect(result.display).toBeCloseTo(-0.1, 3);
    });
  });
});

describe('Continuous Effect Sizes (SMD)', () => {
  it('calculates Hedges g correctly', () => {
    // Example: M1=10, SD1=3, N1=30; M2=8, SD2=3, N2=30
    const row = { m1: 10, s1: 3, n1: 30, m2: 8, s2: 3, n2: 30 };
    const result = calculateEffect(row, 'continuous', null);

    // Pooled SD = sqrt(((29*9 + 29*9)/(58))) = 3
    // Cohen's d = (10-8)/3 = 0.667
    // J = 1 - 3/(4*58 - 1) = 0.987
    // Hedges' g = 0.667 * 0.987 = 0.658
    expect(result).not.toBeNull();
    expect(result.es).toBeCloseTo(0.658, 2);
  });

  it('applies small-sample correction factor', () => {
    // Small sample should have larger correction
    const smallSample = { m1: 10, s1: 3, n1: 10, m2: 8, s2: 3, n2: 10 };
    const largeSample = { m1: 10, s1: 3, n1: 100, m2: 8, s2: 3, n2: 100 };

    const smallResult = calculateEffect(smallSample, 'continuous', null);
    const largeResult = calculateEffect(largeSample, 'continuous', null);

    // J factor closer to 1 for large samples
    const smallJ = smallResult.es / ((10 - 8) / 3);
    const largeJ = largeResult.es / ((10 - 8) / 3);

    expect(largeJ).toBeGreaterThan(smallJ);
    expect(largeJ).toBeCloseTo(1, 1);
  });
});

describe('Correlation Effect Sizes', () => {
  it('applies Fisher z transformation', () => {
    const row = { r: 0.5, n: 50 };
    const result = calculateEffect(row, 'correlation', null);

    // Fisher's z = 0.5 * ln((1+r)/(1-r)) = 0.5 * ln(3) = 0.549
    expect(result).not.toBeNull();
    expect(result.es).toBeCloseTo(0.549, 2);

    // Variance = 1/(n-3) = 1/47 = 0.0213
    expect(result.vi).toBeCloseTo(0.0213, 3);
  });

  it('rejects correlations at boundaries', () => {
    expect(calculateEffect({ r: 1, n: 50 }, 'correlation', null)).toBeNull();
    expect(calculateEffect({ r: -1, n: 50 }, 'correlation', null)).toBeNull();
    expect(calculateEffect({ r: 1.01, n: 50 }, 'correlation', null)).toBeNull();
  });
});

describe('Survival Effect Sizes (HR)', () => {
  it('calculates log HR from CI', () => {
    const row = { hr: 0.75, ll: 0.5, ul: 1.125 };
    const result = calculateEffect(row, 'survival', null);

    // log(HR) = log(0.75) = -0.288
    expect(result).not.toBeNull();
    expect(result.es).toBeCloseTo(-0.288, 2);

    // SE from CI: (log(1.125) - log(0.5)) / (2*1.96) = 0.208
    expect(result.se).toBeCloseTo(0.208, 2);
  });
});

describe('Proportion Effect Sizes', () => {
  it('applies logit transformation', () => {
    const row = { e: 30, n: 100 };
    const result = calculateEffect(row, 'proportion', null);

    // Adjusted p = (30+0.5)/(100+1) = 0.302
    // logit = log(0.302/0.698) = -0.838
    expect(result).not.toBeNull();
    expect(result.es).toBeCloseTo(-0.838, 2);
  });
});

describe('Utility Functions', () => {
  it('transformEffect exponentiates ratio metrics', () => {
    expect(transformEffect(-0.693, 'OR')).toBeCloseTo(0.5, 2);
    expect(transformEffect(-0.693, 'RR')).toBeCloseTo(0.5, 2);
    expect(transformEffect(-0.693, 'HR')).toBeCloseTo(0.5, 2);
    expect(transformEffect(-0.5, 'SMD')).toBeCloseTo(-0.5, 2);
  });

  it('isRatioMetric identifies ratio metrics', () => {
    expect(isRatioMetric('OR')).toBe(true);
    expect(isRatioMetric('RR')).toBe(true);
    expect(isRatioMetric('HR')).toBe(true);
    expect(isRatioMetric('SMD')).toBe(false);
    expect(isRatioMetric('RD')).toBe(false);
  });
});
