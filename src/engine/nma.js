/**
 * Network Meta-Analysis (NMA) module
 * Implements frequentist NMA methods for indirect treatment comparisons
 * References:
 * - Bucher et al. (1997) - Simple indirect comparison
 * - Lumley (2002) - Network meta-analysis
 * - Rücker (2012) - Graph-theoretical approach
 * - Salanti et al. (2008) - Consistency assessment
 */

import { tCritical, pFromZ, chiSqP } from './distributions.js';
import { tauDL, tauREML } from './tau.js';

// ============ DATA STRUCTURES ============

/**
 * Create network structure from study data
 * @param {Array} studies - Array of study objects with treatment comparisons
 * @returns {Object} Network structure
 */
export function createNetwork(studies) {
  const treatments = new Set();
  const comparisons = [];
  const edges = new Map();  // Map of "A:B" -> array of study data

  studies.forEach((study, i) => {
    if (!study.t1 || !study.t2) return;

    const t1 = study.t1;
    const t2 = study.t2;
    treatments.add(t1);
    treatments.add(t2);

    // Normalize edge key (alphabetical order)
    const [tA, tB] = [t1, t2].sort();
    const key = `${tA}:${tB}`;
    const reversed = tA !== t1;

    // Effect size (reverse if needed)
    const es = reversed ? -study.es : study.es;
    const comparison = {
      studyId: study.id || `Study${i + 1}`,
      t1: tA,
      t2: tB,
      es,
      vi: study.vi,
      se: study.se || Math.sqrt(study.vi),
      reversed
    };

    comparisons.push(comparison);

    if (!edges.has(key)) {
      edges.set(key, []);
    }
    edges.get(key).push(comparison);
  });

  const treatmentList = [...treatments].sort();

  // Build adjacency matrix
  const n = treatmentList.length;
  const adjacency = Array(n).fill(null).map(() => Array(n).fill(0));

  edges.forEach((studies, key) => {
    const [t1, t2] = key.split(':');
    const i = treatmentList.indexOf(t1);
    const j = treatmentList.indexOf(t2);
    adjacency[i][j] = studies.length;
    adjacency[j][i] = studies.length;
  });

  return {
    treatments: treatmentList,
    comparisons,
    edges,
    adjacency,
    nTreatments: treatmentList.length,
    nStudies: comparisons.length,
    nComparisons: edges.size
  };
}

/**
 * Check if network is connected (all treatments reachable)
 * @param {Object} network - Network structure
 * @returns {boolean} True if connected
 */
export function isNetworkConnected(network) {
  const { adjacency, nTreatments } = network;
  if (nTreatments <= 1) return true;

  const visited = new Set([0]);
  const queue = [0];

  while (queue.length > 0) {
    const current = queue.shift();
    for (let j = 0; j < nTreatments; j++) {
      if (adjacency[current][j] > 0 && !visited.has(j)) {
        visited.add(j);
        queue.push(j);
      }
    }
  }

  return visited.size === nTreatments;
}

// ============ BUCHER METHOD ============

/**
 * Bucher indirect comparison (simple adjusted indirect comparison)
 * For comparing A vs C through common comparator B: A vs C = (A vs B) - (C vs B)
 * Reference: Bucher et al. (1997)
 *
 * @param {Object} directAB - Direct comparison A vs B { es, vi }
 * @param {Object} directCB - Direct comparison C vs B { es, vi }
 * @param {Object} options - Analysis options
 * @returns {Object} Indirect comparison A vs C
 */
export function bucherIndirect(directAB, directCB, options = {}) {
  const { confLevel = 0.95 } = options;

  // Indirect estimate: d_AC = d_AB - d_CB
  const es = directAB.es - directCB.es;

  // Variance: Var(d_AC) = Var(d_AB) + Var(d_CB)
  const vi = directAB.vi + directCB.vi;
  const se = Math.sqrt(vi);

  // Test statistic
  const z = es / se;
  const pVal = pFromZ(z);

  // Confidence interval
  const alpha = 1 - confLevel;
  const zCrit = Math.abs(normalQuantile(alpha / 2));
  const ciLo = es - zCrit * se;
  const ciHi = es + zCrit * se;

  return {
    method: 'Bucher',
    es,
    se,
    vi,
    z,
    pVal,
    ciLo,
    ciHi,
    confLevel,
    direct: {
      AB: directAB,
      CB: directCB
    }
  };
}

/**
 * Perform all possible Bucher indirect comparisons in a network
 * @param {Object} network - Network structure
 * @param {string} reference - Reference treatment (common comparator)
 * @returns {Array} Array of indirect comparisons
 */
export function allBucherComparisons(network, reference = null) {
  const { treatments, edges } = network;

  // If no reference specified, use the most connected treatment
  if (!reference) {
    const connections = treatments.map(t => {
      let count = 0;
      edges.forEach((_, key) => {
        if (key.includes(t)) count++;
      });
      return { t, count };
    });
    connections.sort((a, b) => b.count - a.count);
    reference = connections[0].t;
  }

  // Find all treatments with direct comparison to reference
  const directToRef = [];
  edges.forEach((studies, key) => {
    const [t1, t2] = key.split(':');
    if (t1 === reference) {
      // Pool direct evidence
      const pooled = poolDirectEvidence(studies);
      directToRef.push({ treatment: t2, ...pooled });
    } else if (t2 === reference) {
      const pooled = poolDirectEvidence(studies, true);  // Reverse
      directToRef.push({ treatment: t1, ...pooled });
    }
  });

  // Generate all pairwise indirect comparisons
  const indirects = [];
  for (let i = 0; i < directToRef.length; i++) {
    for (let j = i + 1; j < directToRef.length; j++) {
      const A = directToRef[i];
      const B = directToRef[j];

      const indirect = bucherIndirect(
        { es: A.es, vi: A.vi },  // A vs Reference
        { es: B.es, vi: B.vi }   // B vs Reference
      );

      indirects.push({
        comparison: `${A.treatment} vs ${B.treatment}`,
        t1: A.treatment,
        t2: B.treatment,
        reference,
        ...indirect
      });
    }
  }

  return indirects;
}

// ============ FREQUENTIST NMA ============

/**
 * Frequentist network meta-analysis using weighted least squares
 * Based on graph-theoretical approach (Rücker 2012)
 *
 * @param {Object} network - Network structure from createNetwork()
 * @param {Object} options - Analysis options
 * @returns {Object} NMA results
 */
export function frequentistNMA(network, options = {}) {
  const {
    reference = null,  // Reference treatment
    tauMethod = 'DL',  // Tau estimator
    confLevel = 0.95,
    small = true       // Small sample correction
  } = options;

  const { treatments, comparisons, edges, nTreatments, nStudies } = network;

  if (nTreatments < 2) {
    return { error: 'Need at least 2 treatments' };
  }

  if (!isNetworkConnected(network)) {
    return { error: 'Network is disconnected - cannot perform NMA' };
  }

  // Set reference treatment (default: first alphabetically)
  const refTreat = reference || treatments[0];
  const refIndex = treatments.indexOf(refTreat);

  // Step 1: Pool direct evidence for each comparison
  const directEstimates = new Map();
  edges.forEach((studies, key) => {
    const pooled = poolDirectEvidence(studies);
    directEstimates.set(key, pooled);
  });

  // Step 2: Estimate between-study variance τ²
  const allEffects = [];
  directEstimates.forEach((est) => {
    if (est.k > 1) {
      allEffects.push(...est.studies.map(s => ({ es: s.es, vi: s.vi })));
    }
  });

  let tau2 = 0;
  if (allEffects.length > 1) {
    tau2 = tauMethod === 'REML' ? tauREML(allEffects) : tauDL(allEffects);
  }

  // Step 3: Build design matrix and weight matrix
  // Design matrix B maps comparisons to treatment effects
  const m = edges.size;  // Number of direct comparisons
  const p = nTreatments - 1;  // Parameters (excluding reference)

  const y = [];  // Observed effect sizes
  const w = [];  // Weights (1/variance)
  const B = [];  // Design matrix (m x p)

  const comparisonList = [];
  edges.forEach((studies, key) => {
    const [t1, t2] = key.split(':');
    const est = directEstimates.get(key);

    y.push(est.es);
    w.push(1 / (est.vi + tau2 / Math.max(1, est.k)));

    // Design matrix row
    const row = Array(p).fill(0);
    const i1 = treatments.indexOf(t1);
    const i2 = treatments.indexOf(t2);

    // Effect is t2 relative to t1: d_t2 - d_t1
    // All effects are relative to reference
    if (i1 !== refIndex && i1 < refIndex) {
      row[i1] = -1;
    } else if (i1 !== refIndex) {
      row[i1 - 1] = -1;
    }

    if (i2 !== refIndex && i2 < refIndex) {
      row[i2] = 1;
    } else if (i2 !== refIndex) {
      row[i2 - 1] = 1;
    }

    B.push(row);
    comparisonList.push({ t1, t2, direct: est });
  });

  // Step 4: Weighted least squares
  // β = (B'WB)^(-1) B'Wy

  // B'WB (p x p)
  const BtWB = Array(p).fill(null).map(() => Array(p).fill(0));
  for (let i = 0; i < p; i++) {
    for (let j = 0; j < p; j++) {
      for (let k = 0; k < m; k++) {
        BtWB[i][j] += B[k][i] * w[k] * B[k][j];
      }
    }
  }

  // B'Wy (p x 1)
  const BtWy = Array(p).fill(0);
  for (let i = 0; i < p; i++) {
    for (let k = 0; k < m; k++) {
      BtWy[i] += B[k][i] * w[k] * y[k];
    }
  }

  // Invert B'WB
  const BtWBinv = invertMatrix(BtWB);
  if (!BtWBinv) {
    return { error: 'Singular design matrix - check network structure' };
  }

  // Calculate treatment effects (relative to reference)
  const beta = Array(p).fill(0);
  for (let i = 0; i < p; i++) {
    for (let j = 0; j < p; j++) {
      beta[i] += BtWBinv[i][j] * BtWy[j];
    }
  }

  // Step 5: Calculate all pairwise comparisons
  const alpha = 1 - confLevel;
  const zCrit = Math.abs(normalQuantile(alpha / 2));

  const treatmentEffects = [];
  treatments.forEach((t, i) => {
    let es = 0;
    let se = 0;

    if (i === refIndex) {
      // Reference treatment
      es = 0;
      se = 0;
    } else {
      const idx = i < refIndex ? i : i - 1;
      es = beta[idx];
      se = Math.sqrt(BtWBinv[idx][idx]);
    }

    treatmentEffects.push({
      treatment: t,
      es,
      se,
      ciLo: es - zCrit * se,
      ciHi: es + zCrit * se,
      isReference: i === refIndex
    });
  });

  // All pairwise comparisons (league table)
  const pairwise = [];
  for (let i = 0; i < nTreatments; i++) {
    for (let j = i + 1; j < nTreatments; j++) {
      const t1 = treatments[i];
      const t2 = treatments[j];

      const idx1 = i < refIndex ? i : (i === refIndex ? -1 : i - 1);
      const idx2 = j < refIndex ? j : (j === refIndex ? -1 : j - 1);

      // d_ij = d_i - d_j (both relative to reference)
      let es, vi;

      if (idx1 === -1) {
        // t1 is reference
        es = -beta[idx2];
        vi = BtWBinv[idx2][idx2];
      } else if (idx2 === -1) {
        // t2 is reference
        es = beta[idx1];
        vi = BtWBinv[idx1][idx1];
      } else {
        es = beta[idx1] - beta[idx2];
        vi = BtWBinv[idx1][idx1] + BtWBinv[idx2][idx2] - 2 * BtWBinv[idx1][idx2];
      }

      const se = Math.sqrt(Math.max(0, vi));
      const z = se > 0 ? es / se : 0;
      const pVal = pFromZ(z);

      // Check for direct evidence
      const directKey = [t1, t2].sort().join(':');
      const hasDirect = edges.has(directKey);
      const directEst = hasDirect ? directEstimates.get(directKey) : null;

      pairwise.push({
        t1,
        t2,
        comparison: `${t1} vs ${t2}`,
        es,
        se,
        ciLo: es - zCrit * se,
        ciHi: es + zCrit * se,
        z,
        pVal,
        hasDirect,
        directEs: directEst?.es,
        directSe: directEst?.se
      });
    }
  }

  // Step 6: Calculate heterogeneity statistics
  // Q statistic for NMA
  const fitted = B.map((row, k) => {
    return row.reduce((sum, bij, j) => sum + bij * beta[j], 0);
  });

  const Q = y.reduce((sum, yi, k) => sum + w[k] * (yi - fitted[k]) ** 2, 0);
  const dfQ = m - p;
  const pQ = dfQ > 0 ? chiSqP(Q, dfQ) : 1;
  const I2 = dfQ > 0 ? Math.max(0, (Q - dfQ) / Q * 100) : 0;

  // Step 7: Treatment rankings (P-scores)
  const pscores = calculatePScores(treatmentEffects, BtWBinv, refIndex);

  return {
    method: 'Frequentist NMA',
    reference: refTreat,
    nTreatments,
    nStudies,
    nComparisons: m,
    tau2,
    tau: Math.sqrt(tau2),
    treatmentEffects,
    pairwise,
    pscores,
    Q,
    dfQ,
    pQ,
    I2,
    confLevel,
    network,
    designMatrix: B,
    weights: w,
    directEstimates: Object.fromEntries(directEstimates)
  };
}

/**
 * Calculate P-scores for treatment ranking
 * P-score = average probability that treatment is better than others
 */
function calculatePScores(effects, covMatrix, refIndex) {
  const n = effects.length;
  const pscores = [];

  effects.forEach((effect, i) => {
    let sumP = 0;
    let count = 0;

    effects.forEach((other, j) => {
      if (i === j) return;

      // Probability that treatment i is better than treatment j
      // P(d_i > d_j) = Φ((d_i - d_j) / SE(d_i - d_j))

      const diff = effect.es - other.es;
      let varDiff;

      // Calculate variance of difference
      const idx_i = i < refIndex ? i : (i === refIndex ? -1 : i - 1);
      const idx_j = j < refIndex ? j : (j === refIndex ? -1 : j - 1);

      if (idx_i === -1 && idx_j === -1) {
        varDiff = 0;
      } else if (idx_i === -1) {
        varDiff = covMatrix[idx_j][idx_j];
      } else if (idx_j === -1) {
        varDiff = covMatrix[idx_i][idx_i];
      } else {
        varDiff = covMatrix[idx_i][idx_i] + covMatrix[idx_j][idx_j] - 2 * (covMatrix[idx_i][idx_j] || 0);
      }

      const seDiff = Math.sqrt(Math.max(0, varDiff));
      const z = seDiff > 0 ? diff / seDiff : (diff > 0 ? Infinity : -Infinity);

      // Assuming lower effect is better (like OR < 1)
      // For higher-is-better, use normalCDF(-z)
      sumP += normalCDF(z);
      count++;
    });

    const pscore = count > 0 ? (sumP / count) * 100 : 50;
    pscores.push({
      treatment: effect.treatment,
      pscore,
      rank: 0  // Will be filled in
    });
  });

  // Assign ranks
  pscores.sort((a, b) => b.pscore - a.pscore);
  pscores.forEach((p, i) => p.rank = i + 1);

  // Sort back by treatment name
  pscores.sort((a, b) => a.treatment.localeCompare(b.treatment));

  return pscores;
}

// ============ INCONSISTENCY ASSESSMENT ============

/**
 * Node-splitting analysis for inconsistency
 * Separates direct and indirect evidence for each comparison
 * Reference: Dias et al. (2010)
 *
 * @param {Object} nmaResult - Result from frequentistNMA()
 * @returns {Array} Node-splitting results for each comparison
 */
export function nodeSplitting(nmaResult) {
  const { network, pairwise, tau2, confLevel } = nmaResult;
  const { edges, treatments } = network;

  const results = [];
  const alpha = 1 - confLevel;
  const zCrit = Math.abs(normalQuantile(alpha / 2));

  // For each comparison with direct evidence
  edges.forEach((studies, key) => {
    const [t1, t2] = key.split(':');

    // Direct estimate
    const direct = poolDirectEvidence(studies);

    // Find the NMA estimate for this comparison
    const nmaEst = pairwise.find(p =>
      (p.t1 === t1 && p.t2 === t2) || (p.t1 === t2 && p.t2 === t1)
    );

    if (!nmaEst) return;

    // Indirect estimate = NMA estimate without direct evidence
    // Using Bucher-type formula: indirect = NMA - (w_direct/w_total) * (direct - NMA)
    // Simplified: indirect = (NMA * w_total - direct * w_direct) / w_indirect

    const directVar = direct.vi + tau2 / Math.max(1, direct.k);
    const nmaVar = nmaEst.se ** 2;

    // Assuming consistent model, derive indirect
    // w_NMA = w_direct + w_indirect
    // When network is simple: indirect ≈ NMA (weighted combination)

    // For proper indirect, we'd need to refit NMA excluding this comparison
    // Approximation: use difference in precision
    const indirectVar = Math.max(nmaVar, directVar * 2);  // Conservative
    const indirectEs = (nmaEst.es * (1/nmaVar) - direct.es * (1/directVar)) /
                       (1/nmaVar - 1/directVar + 1e-10);
    const indirectSe = Math.sqrt(indirectVar);

    // Difference (inconsistency)
    const diff = direct.es - indirectEs;
    const diffVar = directVar + indirectVar;
    const diffSe = Math.sqrt(diffVar);
    const z = diff / diffSe;
    const pVal = pFromZ(z);

    results.push({
      comparison: `${t1} vs ${t2}`,
      t1,
      t2,
      direct: {
        es: direct.es,
        se: Math.sqrt(directVar),
        k: direct.k
      },
      indirect: {
        es: indirectEs,
        se: indirectSe
      },
      nma: {
        es: nmaEst.es,
        se: nmaEst.se
      },
      difference: {
        es: diff,
        se: diffSe,
        ciLo: diff - zCrit * diffSe,
        ciHi: diff + zCrit * diffSe,
        z,
        pVal
      },
      inconsistent: pVal < 0.1  // Common threshold
    });
  });

  return results;
}

/**
 * Global inconsistency test (Q statistic for inconsistency)
 * @param {Object} nmaResult - NMA result
 * @returns {Object} Inconsistency test result
 */
export function globalInconsistencyTest(nmaResult) {
  const { Q, dfQ, pQ, I2, network } = nmaResult;

  // Design-by-treatment interaction test
  // Q_inconsistency = Q_total - Q_heterogeneity

  // For simple networks, Q_heterogeneity ≈ sum of within-comparison Q
  let QHet = 0;
  const { edges } = network;

  edges.forEach((studies) => {
    if (studies.length > 1) {
      const pooled = poolDirectEvidence(studies);
      QHet += pooled.Q || 0;
    }
  });

  const QIncon = Math.max(0, Q - QHet);
  const dfHet = nmaResult.nStudies - edges.size;
  const dfIncon = dfQ - dfHet;
  const pIncon = dfIncon > 0 ? chiSqP(QIncon, dfIncon) : 1;

  return {
    QTotal: Q,
    QHeterogeneity: QHet,
    QInconsistency: QIncon,
    dfTotal: dfQ,
    dfHeterogeneity: dfHet,
    dfInconsistency: dfIncon,
    pHeterogeneity: dfHet > 0 ? chiSqP(QHet, dfHet) : 1,
    pInconsistency: pIncon,
    I2: I2,
    hasInconsistency: pIncon < 0.1
  };
}

/**
 * Generate league table (all pairwise comparisons matrix)
 * @param {Object} nmaResult - NMA result
 * @returns {Object} League table data
 */
export function generateLeagueTable(nmaResult) {
  const { network, pairwise, treatmentEffects, pscores } = nmaResult;
  const treatments = network?.treatments || [];
  const n = treatments.length;

  // Create n x n matrix
  const table = Array(n).fill(null).map(() => Array(n).fill(null));

  // Fill diagonal with treatment names/P-scores
  treatments.forEach((t, i) => {
    const pscore = pscores.find(p => p.treatment === t);
    table[i][i] = {
      treatment: t,
      pscore: pscore?.pscore,
      rank: pscore?.rank
    };
  });

  // Fill upper triangle with comparisons (row vs column)
  pairwise.forEach(comp => {
    const i = treatments.indexOf(comp.t1);
    const j = treatments.indexOf(comp.t2);

    // Upper triangle: row beats column
    table[i][j] = {
      es: comp.es,
      se: comp.se,
      ciLo: comp.ciLo,
      ciHi: comp.ciHi,
      hasDirect: comp.hasDirect
    };

    // Lower triangle: reversed (column beats row)
    table[j][i] = {
      es: -comp.es,
      se: comp.se,
      ciLo: -comp.ciHi,
      ciHi: -comp.ciLo,
      hasDirect: comp.hasDirect
    };
  });

  return {
    treatments,
    matrix: table,
    pscores
  };
}

// ============ NETWORK VISUALIZATION DATA ============

/**
 * Generate data for network diagram visualization
 * @param {Object} network - Network structure
 * @returns {Object} Visualization data
 */
export function getNetworkDiagramData(network) {
  const { treatments, edges, adjacency, nTreatments } = network;

  // Calculate node positions (circular layout)
  const nodes = treatments.map((t, i) => {
    const angle = (2 * Math.PI * i) / nTreatments - Math.PI / 2;
    const radius = 100;

    return {
      id: t,
      label: t,
      x: radius * Math.cos(angle),
      y: radius * Math.sin(angle),
      // Size based on number of connections
      size: [...edges.keys()].filter(k => k.includes(t)).reduce((sum, k) => sum + edges.get(k).length, 0)
    };
  });

  // Create edge data
  const edgeData = [];
  edges.forEach((studies, key) => {
    const [t1, t2] = key.split(':');
    const pooled = poolDirectEvidence(studies);

    edgeData.push({
      source: t1,
      target: t2,
      weight: studies.length,  // Number of studies
      es: pooled.es,
      se: Math.sqrt(pooled.vi),
      k: studies.length
    });
  });

  return { nodes, edges: edgeData };
}

// ============ HELPER FUNCTIONS ============

/**
 * Pool direct evidence for a comparison using random effects
 */
function poolDirectEvidence(studies, reverse = false) {
  if (studies.length === 0) return null;

  const effects = studies.map(s => ({
    es: reverse ? -s.es : s.es,
    vi: s.vi
  }));

  if (effects.length === 1) {
    return {
      es: effects[0].es,
      vi: effects[0].vi,
      se: Math.sqrt(effects[0].vi),
      k: 1,
      Q: 0,
      studies: effects
    };
  }

  // DL random effects pooling
  const tau2 = tauDL(effects);
  const w = effects.map(e => 1 / (e.vi + tau2));
  const sumW = w.reduce((a, b) => a + b, 0);

  const es = effects.reduce((acc, e, i) => acc + w[i] * e.es, 0) / sumW;
  const vi = 1 / sumW;

  // Q statistic
  const Q = effects.reduce((acc, e, i) => acc + w[i] * (e.es - es) ** 2, 0);

  return {
    es,
    vi,
    se: Math.sqrt(vi),
    k: effects.length,
    Q,
    tau2,
    studies: effects
  };
}

/**
 * Matrix inversion using Gauss-Jordan elimination
 */
function invertMatrix(matrix) {
  const n = matrix.length;
  if (n === 0) return null;

  const aug = matrix.map((row, i) =>
    [...row, ...Array(n).fill(0).map((_, j) => i === j ? 1 : 0)]
  );

  for (let i = 0; i < n; i++) {
    let maxRow = i;
    for (let k = i + 1; k < n; k++) {
      if (Math.abs(aug[k][i]) > Math.abs(aug[maxRow][i])) maxRow = k;
    }
    [aug[i], aug[maxRow]] = [aug[maxRow], aug[i]];

    if (Math.abs(aug[i][i]) < 1e-10) return null;

    const pivot = aug[i][i];
    for (let j = 0; j < 2 * n; j++) aug[i][j] /= pivot;

    for (let k = 0; k < n; k++) {
      if (k !== i) {
        const factor = aug[k][i];
        for (let j = 0; j < 2 * n; j++) aug[k][j] -= factor * aug[i][j];
      }
    }
  }

  return aug.map(row => row.slice(n));
}

/**
 * Normal CDF
 */
function normalCDF(x) {
  const a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741;
  const a4 = -1.453152027, a5 = 1.061405429, p = 0.3275911;

  const sign = x < 0 ? -1 : 1;
  x = Math.abs(x) / Math.sqrt(2);
  const t = 1 / (1 + p * x);
  const y = 1 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-x * x);

  return 0.5 * (1 + sign * y);
}

/**
 * Normal quantile (inverse CDF)
 */
function normalQuantile(p) {
  if (p <= 0) return -Infinity;
  if (p >= 1) return Infinity;
  if (p === 0.5) return 0;

  const a = [
    -3.969683028665376e+01, 2.209460984245205e+02,
    -2.759285104469687e+02, 1.383577518672690e+02,
    -3.066479806614716e+01, 2.506628277459239e+00
  ];
  const b = [
    -5.447609879822406e+01, 1.615858368580409e+02,
    -1.556989798598866e+02, 6.680131188771972e+01,
    -1.328068155288572e+01
  ];
  const c = [
    -7.784894002430293e-03, -3.223964580411365e-01,
    -2.400758277161838e+00, -2.549732539343734e+00,
    4.374664141464968e+00, 2.938163982698783e+00
  ];
  const d = [
    7.784695709041462e-03, 3.224671290700398e-01,
    2.445134137142996e+00, 3.754408661907416e+00
  ];

  const pLow = 0.02425;
  const pHigh = 1 - pLow;
  let q, r;

  if (p < pLow) {
    q = Math.sqrt(-2 * Math.log(p));
    return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
           ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
  } else if (p <= pHigh) {
    q = p - 0.5;
    r = q * q;
    return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q /
           (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
  } else {
    q = Math.sqrt(-2 * Math.log(1 - p));
    return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) /
            ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
  }
}
