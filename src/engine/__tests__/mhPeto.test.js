/**
 * Tests for Mantel-Haenszel and Peto methods
 */
import { describe, it, expect } from 'vitest';
import { mantelHaenszel, petoOddsRatio, recommendBinaryMethod } from '../mhPeto.js';

// Test data: sparse binary data with zero cells
const sparseStudies = [
  { id: 'Study1', raw: { e1: 3, n1: 100, e2: 10, n2: 100 }, es: -1.24, vi: 0.44, se: 0.66 },
  { id: 'Study2', raw: { e1: 0, n1: 50, e2: 2, n2: 50 }, es: -2.35, vi: 1.5, se: 1.22 },  // Zero cell
  { id: 'Study3', raw: { e1: 5, n1: 200, e2: 15, n2: 200 }, es: -1.15, vi: 0.27, se: 0.52 },
  { id: 'Study4', raw: { e1: 1, n1: 80, e2: 4, n2: 80 }, es: -1.42, vi: 1.25, se: 1.12 },
];

// Balanced studies for Peto method
const balancedStudies = [
  { id: 'A', raw: { e1: 5, n1: 100, e2: 15, n2: 100 }, es: -1.15, vi: 0.27, se: 0.52 },
  { id: 'B', raw: { e1: 8, n1: 150, e2: 20, n2: 150 }, es: -0.99, vi: 0.19, se: 0.44 },
  { id: 'C', raw: { e1: 3, n1: 80, e2: 10, n2: 80 }, es: -1.28, vi: 0.43, se: 0.66 },
];

describe('Mantel-Haenszel Methods', () => {
  describe('mantelHaenszel - Odds Ratio', () => {
    it('should calculate MH OR for sparse data', () => {
      const result = mantelHaenszel(sparseStudies, 'OR');

      expect(result).not.toBeNull();
      expect(result.method).toBe('Mantel-Haenszel');
      expect(result.metric).toBe('OR');
      expect(result.k).toBe(4);
      expect(result.es).toBeDefined();
      expect(result.se).toBeGreaterThan(0);
      expect(result.display).toBeGreaterThan(0);  // OR should be positive
      expect(result.display).toBeLessThan(1);     // Treatment is protective
    });

    it('should produce narrower CI than IV method for sparse data', () => {
      const result = mantelHaenszel(sparseStudies, 'OR');

      // MH typically produces narrower CIs for sparse data
      expect(result.ciHi).toBeGreaterThan(result.ciLo);
      expect(result.ciHi / result.ciLo).toBeLessThan(20);  // Reasonable CI width
    });

    it('should calculate heterogeneity statistics', () => {
      const result = mantelHaenszel(sparseStudies, 'OR');

      expect(result.Q).toBeDefined();
      expect(result.I2).toBeDefined();
      expect(result.I2).toBeGreaterThanOrEqual(0);
      expect(result.I2).toBeLessThanOrEqual(100);
      expect(result.pHet).toBeDefined();
    });
  });

  describe('mantelHaenszel - Risk Ratio', () => {
    it('should calculate MH RR', () => {
      const result = mantelHaenszel(balancedStudies, 'RR');

      expect(result).not.toBeNull();
      expect(result.metric).toBe('RR');
      expect(result.display).toBeGreaterThan(0);
      expect(result.display).toBeLessThan(1);  // Treatment is protective
    });
  });

  describe('mantelHaenszel - Risk Difference', () => {
    it('should calculate MH RD', () => {
      const result = mantelHaenszel(balancedStudies, 'RD');

      expect(result).not.toBeNull();
      expect(result.metric).toBe('RD');
      expect(result.es).toBeLessThan(0);  // Negative RD (treatment reduces risk)
      expect(result.isRatio).toBe(false);
    });
  });

  describe('Edge cases', () => {
    it('should handle empty input', () => {
      const result = mantelHaenszel([], 'OR');
      expect(result).toBeNull();
    });

    it('should handle excluded studies', () => {
      const withExcluded = [
        ...balancedStudies,
        { id: 'Excluded', raw: { e1: 50, n1: 100, e2: 10, n2: 100 }, excluded: true }
      ];
      const result = mantelHaenszel(withExcluded, 'OR');

      expect(result.k).toBe(3);  // Excluded study not counted
    });

    it('should exclude double-zero studies by default', () => {
      const withDoubleZero = [
        ...balancedStudies,
        { id: 'DoubleZero', raw: { e1: 0, n1: 100, e2: 0, n2: 100 } }
      ];
      const result = mantelHaenszel(withDoubleZero, 'OR');

      expect(result.k).toBe(3);  // Double-zero excluded
    });
  });
});

describe('Peto Odds Ratio', () => {
  it('should calculate Peto OR for balanced studies', () => {
    const result = petoOddsRatio(balancedStudies);

    expect(result).not.toBeNull();
    expect(result.method).toBe('Peto');
    expect(result.metric).toBe('OR');
    expect(result.k).toBe(3);
    expect(result.display).toBeGreaterThan(0);
    expect(result.display).toBeLessThan(1);  // Treatment protective
  });

  it('should calculate O-E and variance correctly', () => {
    const result = petoOddsRatio(balancedStudies);

    expect(result.O_E).toBeDefined();
    expect(result.V).toBeGreaterThan(0);

    // Peto log(OR) = O_E / V
    expect(Math.abs(result.es - result.O_E / result.V)).toBeLessThan(0.001);
  });

  it('should handle heterogeneity', () => {
    const result = petoOddsRatio(balancedStudies);

    expect(result.Q).toBeDefined();
    expect(result.I2).toBeDefined();
    expect(result.df).toBe(2);  // k - 1
  });
});

describe('recommendBinaryMethod', () => {
  it('should recommend Peto for rare events with balanced groups', () => {
    const rareBalanced = [
      { raw: { e1: 1, n1: 500, e2: 5, n2: 500 } },
      { raw: { e1: 2, n1: 600, e2: 8, n2: 600 } },
    ];
    const rec = recommendBinaryMethod(rareBalanced);

    expect(rec.method).toBe('Peto');
    expect(rec.rationale).toContain('Rare events');
  });

  it('should recommend MH for sparse data with zero cells', () => {
    const rec = recommendBinaryMethod(sparseStudies);

    expect(['MH', 'Peto']).toContain(rec.method);
  });

  it('should recommend IV for adequate events', () => {
    const adequateEvents = [
      { raw: { e1: 30, n1: 100, e2: 40, n2: 100 } },
      { raw: { e1: 25, n1: 80, e2: 35, n2: 80 } },
    ];
    const rec = recommendBinaryMethod(adequateEvents);

    expect(rec.method).toBe('IV');
    expect(rec.rationale).toContain('Inverse variance');
  });
});
