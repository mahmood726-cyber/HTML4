/**
 * Tests for Network Meta-Analysis
 */
import { describe, it, expect } from 'vitest';
import {
  createNetwork,
  isNetworkConnected,
  bucherIndirect,
  frequentistNMA,
  nodeSplitting,
  globalInconsistencyTest,
  generateLeagueTable,
  getNetworkDiagramData
} from '../nma.js';

// Test network: A vs B vs C (star network)
const starNetwork = [
  { id: 'S1', t1: 'A', t2: 'B', es: -0.5, vi: 0.1, se: 0.316 },
  { id: 'S2', t1: 'A', t2: 'B', es: -0.6, vi: 0.12, se: 0.346 },
  { id: 'S3', t1: 'A', t2: 'C', es: -0.3, vi: 0.08, se: 0.283 },
  { id: 'S4', t1: 'A', t2: 'C', es: -0.4, vi: 0.09, se: 0.300 },
];

// Triangle network (A-B-C with all direct comparisons)
const triangleNetwork = [
  { id: 'S1', t1: 'A', t2: 'B', es: -0.5, vi: 0.1, se: 0.316 },
  { id: 'S2', t1: 'A', t2: 'C', es: -0.3, vi: 0.08, se: 0.283 },
  { id: 'S3', t1: 'B', t2: 'C', es: 0.2, vi: 0.09, se: 0.300 },
];

// Disconnected network
const disconnectedNetwork = [
  { id: 'S1', t1: 'A', t2: 'B', es: -0.5, vi: 0.1, se: 0.316 },
  { id: 'S2', t1: 'C', t2: 'D', es: -0.3, vi: 0.08, se: 0.283 },
];

describe('Network Creation', () => {
  it('should create network structure from studies', () => {
    const network = createNetwork(starNetwork);

    expect(network.treatments).toEqual(['A', 'B', 'C']);
    expect(network.nTreatments).toBe(3);
    expect(network.nStudies).toBe(4);
    expect(network.nComparisons).toBe(2);  // A-B and A-C
  });

  it('should build adjacency matrix', () => {
    const network = createNetwork(starNetwork);

    expect(network.adjacency).toBeDefined();
    expect(network.adjacency.length).toBe(3);

    // A connected to B and C
    const aIdx = network.treatments.indexOf('A');
    const bIdx = network.treatments.indexOf('B');
    const cIdx = network.treatments.indexOf('C');

    expect(network.adjacency[aIdx][bIdx]).toBe(2);  // 2 studies
    expect(network.adjacency[aIdx][cIdx]).toBe(2);  // 2 studies
    expect(network.adjacency[bIdx][cIdx]).toBe(0);  // No direct
  });

  it('should detect connected network', () => {
    const network = createNetwork(starNetwork);
    expect(isNetworkConnected(network)).toBe(true);
  });

  it('should detect disconnected network', () => {
    const network = createNetwork(disconnectedNetwork);
    expect(isNetworkConnected(network)).toBe(false);
  });
});

describe('Bucher Indirect Comparison', () => {
  it('should calculate indirect comparison', () => {
    const directAB = { es: -0.5, vi: 0.1 };  // A vs B
    const directCB = { es: 0.2, vi: 0.09 };  // C vs B

    // A vs C = (A vs B) - (C vs B)
    const indirect = bucherIndirect(directAB, directCB);

    expect(indirect.method).toBe('Bucher');
    expect(indirect.es).toBeCloseTo(-0.7, 2);  // -0.5 - 0.2
    expect(indirect.vi).toBeCloseTo(0.19, 2);  // 0.1 + 0.09
    expect(indirect.se).toBeCloseTo(Math.sqrt(0.19), 2);
  });

  it('should calculate confidence intervals', () => {
    const directAB = { es: -0.5, vi: 0.1 };
    const directCB = { es: 0.2, vi: 0.09 };
    const indirect = bucherIndirect(directAB, directCB, { confLevel: 0.95 });

    expect(indirect.ciLo).toBeDefined();
    expect(indirect.ciHi).toBeDefined();
    expect(indirect.ciLo).toBeLessThan(indirect.es);
    expect(indirect.ciHi).toBeGreaterThan(indirect.es);
  });
});

describe('Frequentist NMA', () => {
  it('should run NMA on star network', () => {
    const network = createNetwork(starNetwork);
    const result = frequentistNMA(network, { reference: 'A' });

    expect(result.error).toBeUndefined();
    expect(result.method).toBe('Frequentist NMA');
    expect(result.reference).toBe('A');
    expect(result.nTreatments).toBe(3);
  });

  it('should calculate treatment effects relative to reference', () => {
    const network = createNetwork(starNetwork);
    const result = frequentistNMA(network, { reference: 'A' });

    expect(result.treatmentEffects.length).toBe(3);

    // Reference should have effect = 0
    const refEffect = result.treatmentEffects.find(t => t.treatment === 'A');
    expect(refEffect.es).toBe(0);
    expect(refEffect.isReference).toBe(true);

    // Other treatments should have non-zero effects
    const bEffect = result.treatmentEffects.find(t => t.treatment === 'B');
    expect(bEffect.es).not.toBe(0);
  });

  it('should calculate all pairwise comparisons', () => {
    const network = createNetwork(triangleNetwork);
    const result = frequentistNMA(network);

    // 3 treatments = 3 pairwise comparisons
    expect(result.pairwise.length).toBe(3);

    // Check A vs B
    const ab = result.pairwise.find(p =>
      (p.t1 === 'A' && p.t2 === 'B') || (p.t1 === 'B' && p.t2 === 'A')
    );
    expect(ab).toBeDefined();
    expect(ab.hasDirect).toBe(true);
  });

  it('should calculate P-scores for ranking', () => {
    const network = createNetwork(starNetwork);
    const result = frequentistNMA(network);

    expect(result.pscores.length).toBe(3);
    result.pscores.forEach(p => {
      expect(p.pscore).toBeGreaterThanOrEqual(0);
      expect(p.pscore).toBeLessThanOrEqual(100);
      expect(p.rank).toBeGreaterThanOrEqual(1);
      expect(p.rank).toBeLessThanOrEqual(3);
    });
  });

  it('should fail on disconnected network', () => {
    const network = createNetwork(disconnectedNetwork);
    const result = frequentistNMA(network);

    expect(result.error).toBeDefined();
    expect(result.error).toContain('disconnected');
  });

  it('should calculate heterogeneity statistics', () => {
    const network = createNetwork(starNetwork);
    const result = frequentistNMA(network);

    expect(result.Q).toBeDefined();
    expect(result.I2).toBeDefined();
    expect(result.tau2).toBeDefined();
  });
});

describe('Inconsistency Assessment', () => {
  it('should perform node-splitting analysis', () => {
    const network = createNetwork(triangleNetwork);
    const nma = frequentistNMA(network);
    const splits = nodeSplitting(nma);

    // Should have results for each comparison with direct evidence
    expect(splits.length).toBeGreaterThan(0);

    splits.forEach(split => {
      expect(split.direct.es).toBeDefined();
      expect(split.indirect.es).toBeDefined();
      expect(split.difference.pVal).toBeDefined();
    });
  });

  it('should perform global inconsistency test', () => {
    const network = createNetwork(triangleNetwork);
    const nma = frequentistNMA(network);
    const global = globalInconsistencyTest(nma);

    expect(global.QTotal).toBeDefined();
    expect(global.QInconsistency).toBeDefined();
    expect(global.pInconsistency).toBeDefined();
    expect(global.hasInconsistency).toBeDefined();
  });
});

describe('League Table', () => {
  it('should generate league table matrix', () => {
    const network = createNetwork(triangleNetwork);
    const nma = frequentistNMA(network);
    const table = generateLeagueTable(nma);

    expect(table.treatments.length).toBe(3);
    expect(table.matrix.length).toBe(3);
    expect(table.matrix[0].length).toBe(3);

    // Diagonal should contain treatment info
    table.treatments.forEach((t, i) => {
      expect(table.matrix[i][i].treatment).toBe(t);
    });
  });
});

describe('Network Diagram Data', () => {
  it('should generate visualization data', () => {
    const network = createNetwork(starNetwork);
    const viz = getNetworkDiagramData(network);

    expect(viz.nodes.length).toBe(3);
    expect(viz.edges.length).toBe(2);

    viz.nodes.forEach(node => {
      expect(node.id).toBeDefined();
      expect(node.x).toBeDefined();
      expect(node.y).toBeDefined();
    });

    viz.edges.forEach(edge => {
      expect(edge.source).toBeDefined();
      expect(edge.target).toBeDefined();
      expect(edge.weight).toBeGreaterThan(0);
    });
  });
});
