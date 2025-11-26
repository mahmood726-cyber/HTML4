/**
 * Component Network Meta-Analysis (CNMA)
 * For analyzing complex interventions composed of multiple components
 * References:
 * - Welton et al. (2009) - Component NMA concept
 * - Rücker et al. (2020) - Additive CNMA
 * - Tsokani et al. (2024) - Cochrane webinar on CNMA
 */

import { chiSqP, pFromZ } from './distributions.js';

/**
 * Additive Component NMA
 * Assumes treatment effects are sum of component effects: d_AB = d_A + d_B
 *
 * @param {Array} studies - Array of study comparisons
 *        Each: { id, t1: ['A'], t2: ['A', 'B'], es, vi } where treatments are arrays of components
 * @param {Object} options - Analysis options
 * @returns {Object} CNMA results
 */
export function additiveCNMA(studies, options = {}) {
  const {
    reference = null,  // Reference (usually placebo/control with no components)
    confLevel = 0.95,
    includeInteractions = false
  } = options;

  // Extract all unique components
  const allComponents = new Set();
  studies.forEach(s => {
    (s.t1 || []).forEach(c => allComponents.add(c));
    (s.t2 || []).forEach(c => allComponents.add(c));
  });

  const components = [...allComponents].sort();
  const nComponents = components.length;

  if (nComponents === 0) {
    return { error: 'No components found in study data' };
  }

  // Build design matrix
  // Each row encodes the contrast: (components in t2) - (components in t1)
  const validStudies = studies.filter(s => s.es !== null && s.vi > 0);
  const m = validStudies.length;

  if (m < nComponents) {
    return { error: `Need at least ${nComponents} studies for ${nComponents} components` };
  }

  const y = [];  // Effect sizes
  const w = [];  // Weights
  const B = [];  // Design matrix

  validStudies.forEach(s => {
    y.push(s.es);
    w.push(1 / s.vi);

    // Create design row
    const row = Array(nComponents).fill(0);
    const t1Components = s.t1 || [];
    const t2Components = s.t2 || [];

    components.forEach((c, i) => {
      const inT1 = t1Components.includes(c) ? 1 : 0;
      const inT2 = t2Components.includes(c) ? 1 : 0;
      row[i] = inT2 - inT1;  // +1 if added in t2, -1 if removed, 0 if same
    });

    B.push(row);
  });

  // Check for identifiability
  // Design matrix must have full column rank
  const rank = estimateRank(B);
  if (rank < nComponents) {
    return {
      error: 'Network is not identifiable - some components cannot be estimated',
      identifiableComponents: rank,
      totalComponents: nComponents
    };
  }

  // Weighted least squares: β = (B'WB)⁻¹ B'Wy
  const BtWB = Array(nComponents).fill(null).map(() => Array(nComponents).fill(0));
  const BtWy = Array(nComponents).fill(0);

  for (let i = 0; i < nComponents; i++) {
    for (let j = 0; j < nComponents; j++) {
      for (let k = 0; k < m; k++) {
        BtWB[i][j] += B[k][i] * w[k] * B[k][j];
      }
    }
    for (let k = 0; k < m; k++) {
      BtWy[i] += B[k][i] * w[k] * y[k];
    }
  }

  const BtWBinv = invertMatrix(BtWB);
  if (!BtWBinv) {
    return { error: 'Singular design matrix - check component definitions' };
  }

  // Component effects
  const beta = Array(nComponents).fill(0);
  for (let i = 0; i < nComponents; i++) {
    for (let j = 0; j < nComponents; j++) {
      beta[i] += BtWBinv[i][j] * BtWy[j];
    }
  }

  // Standard errors and CIs
  const alpha = 1 - confLevel;
  const zCrit = Math.abs(normalQuantile(alpha / 2));

  const componentEffects = components.map((c, i) => {
    const se = Math.sqrt(BtWBinv[i][i]);
    const z = beta[i] / se;
    const pVal = pFromZ(z);

    return {
      component: c,
      es: beta[i],
      se,
      ciLo: beta[i] - zCrit * se,
      ciHi: beta[i] + zCrit * se,
      z,
      pVal
    };
  });

  // Fitted values and residuals
  const fitted = B.map((row) => row.reduce((sum, bij, j) => sum + bij * beta[j], 0));
  const residuals = y.map((yi, i) => yi - fitted[i]);

  // Q statistic for model fit
  const Q = residuals.reduce((sum, r, i) => sum + w[i] * r * r, 0);
  const dfQ = m - nComponents;
  const pQ = dfQ > 0 ? chiSqP(Q, dfQ) : 1;
  const I2 = dfQ > 0 ? Math.max(0, (Q - dfQ) / Q * 100) : 0;

  // Generate all possible treatment combinations
  const treatments = generateTreatmentCombinations(components, validStudies);

  // Calculate effects for each treatment combination
  const treatmentEffects = treatments.map(t => {
    let es = 0;
    let varES = 0;

    t.components.forEach((c, i) => {
      const idx = components.indexOf(c);
      if (idx >= 0) {
        es += beta[idx];
        // Add variance contribution
        varES += BtWBinv[idx][idx];

        // Add covariances
        t.components.forEach((c2, j) => {
          if (j > i) {
            const idx2 = components.indexOf(c2);
            if (idx2 >= 0) {
              varES += 2 * BtWBinv[idx][idx2];
            }
          }
        });
      }
    });

    const se = Math.sqrt(Math.max(0, varES));

    return {
      treatment: t.label,
      components: t.components,
      es,
      se,
      ciLo: es - zCrit * se,
      ciHi: es + zCrit * se
    };
  });

  // Component rankings
  const rankings = componentEffects
    .map(c => ({ ...c }))
    .sort((a, b) => a.es - b.es)  // Lower is better for OR
    .map((c, i) => ({ ...c, rank: i + 1 }));

  return {
    method: 'Additive Component NMA',
    nStudies: m,
    nComponents,
    components,
    componentEffects,
    rankings,
    treatmentEffects,
    modelFit: {
      Q,
      dfQ,
      pQ,
      I2
    },
    confLevel,
    designMatrix: B,
    fitted,
    residuals
  };
}

/**
 * Interaction Component NMA
 * Allows for synergistic or antagonistic effects between components
 *
 * @param {Array} studies - Study data
 * @param {Array<Array>} interactions - Pairs of components to model interactions for
 * @param {Object} options - Analysis options
 */
export function interactionCNMA(studies, interactions = [], options = {}) {
  const { confLevel = 0.95 } = options;

  // Get base additive model
  const baseResult = additiveCNMA(studies, { ...options, includeInteractions: false });
  if (baseResult.error) return baseResult;

  const { components, componentEffects } = baseResult;
  const validStudies = studies.filter(s => s.es !== null && s.vi > 0);
  const m = validStudies.length;

  // If no interactions specified, detect potential interactions
  if (interactions.length === 0) {
    interactions = detectPotentialInteractions(validStudies, components);
  }

  const nInteractions = interactions.length;
  if (nInteractions === 0) {
    return { ...baseResult, interactions: [], message: 'No interactions modeled' };
  }

  const nParams = components.length + nInteractions;
  if (m < nParams) {
    return {
      ...baseResult,
      interactions: [],
      warning: `Insufficient studies (${m}) to estimate ${nInteractions} interactions`
    };
  }

  // Build extended design matrix
  const y = validStudies.map(s => s.es);
  const w = validStudies.map(s => 1 / s.vi);
  const B = [];

  validStudies.forEach(s => {
    const t1 = s.t1 || [];
    const t2 = s.t2 || [];

    // Main effects
    const row = components.map(c => {
      const in1 = t1.includes(c) ? 1 : 0;
      const in2 = t2.includes(c) ? 1 : 0;
      return in2 - in1;
    });

    // Interaction effects
    interactions.forEach(([c1, c2]) => {
      const hasInt1 = t1.includes(c1) && t1.includes(c2) ? 1 : 0;
      const hasInt2 = t2.includes(c1) && t2.includes(c2) ? 1 : 0;
      row.push(hasInt2 - hasInt1);
    });

    B.push(row);
  });

  // WLS estimation
  const BtWB = Array(nParams).fill(null).map(() => Array(nParams).fill(0));
  const BtWy = Array(nParams).fill(0);

  for (let i = 0; i < nParams; i++) {
    for (let j = 0; j < nParams; j++) {
      for (let k = 0; k < m; k++) {
        BtWB[i][j] += B[k][i] * w[k] * B[k][j];
      }
    }
    for (let k = 0; k < m; k++) {
      BtWy[i] += B[k][i] * w[k] * y[k];
    }
  }

  const BtWBinv = invertMatrix(BtWB);
  if (!BtWBinv) {
    return { ...baseResult, interactions: [], warning: 'Could not estimate interactions - singular matrix' };
  }

  const beta = Array(nParams).fill(0);
  for (let i = 0; i < nParams; i++) {
    for (let j = 0; j < nParams; j++) {
      beta[i] += BtWBinv[i][j] * BtWy[j];
    }
  }

  const alpha = 1 - confLevel;
  const zCrit = Math.abs(normalQuantile(alpha / 2));

  // Component effects (updated)
  const updatedComponentEffects = components.map((c, i) => {
    const se = Math.sqrt(BtWBinv[i][i]);
    return {
      component: c,
      es: beta[i],
      se,
      ciLo: beta[i] - zCrit * se,
      ciHi: beta[i] + zCrit * se,
      z: beta[i] / se,
      pVal: pFromZ(beta[i] / se)
    };
  });

  // Interaction effects
  const interactionEffects = interactions.map(([c1, c2], i) => {
    const idx = components.length + i;
    const se = Math.sqrt(BtWBinv[idx][idx]);
    const effect = beta[idx];

    return {
      interaction: `${c1} × ${c2}`,
      components: [c1, c2],
      es: effect,
      se,
      ciLo: effect - zCrit * se,
      ciHi: effect + zCrit * se,
      z: effect / se,
      pVal: pFromZ(effect / se),
      type: effect > 0 ? 'antagonistic' : (effect < 0 ? 'synergistic' : 'none')
    };
  });

  // Model fit
  const fitted = B.map((row) => row.reduce((sum, bij, j) => sum + bij * beta[j], 0));
  const Q = validStudies.reduce((sum, s, i) => sum + w[i] * (y[i] - fitted[i]) ** 2, 0);
  const dfQ = m - nParams;

  // Compare to additive model
  const QDiff = baseResult.modelFit.Q - Q;
  const dfDiff = nInteractions;
  const pInteraction = chiSqP(QDiff, dfDiff);

  return {
    method: 'Interaction Component NMA',
    nStudies: m,
    nComponents: components.length,
    nInteractions,
    components,
    componentEffects: updatedComponentEffects,
    interactionEffects,
    modelFit: {
      Q,
      dfQ,
      pQ: dfQ > 0 ? chiSqP(Q, dfQ) : 1,
      I2: dfQ > 0 ? Math.max(0, (Q - dfQ) / Q * 100) : 0
    },
    interactionTest: {
      QDiff,
      dfDiff,
      pValue: pInteraction,
      significant: pInteraction < 0.05
    },
    additiveModel: baseResult,
    confLevel
  };
}

/**
 * Detect potential interactions from data patterns
 */
function detectPotentialInteractions(studies, components) {
  const interactions = [];

  // Find component pairs that appear together
  const pairCounts = {};

  studies.forEach(s => {
    const allComponents = [...(s.t1 || []), ...(s.t2 || [])];
    const unique = [...new Set(allComponents)];

    for (let i = 0; i < unique.length; i++) {
      for (let j = i + 1; j < unique.length; j++) {
        const pair = [unique[i], unique[j]].sort().join('|');
        pairCounts[pair] = (pairCounts[pair] || 0) + 1;
      }
    }
  });

  // Include pairs that appear in at least 2 studies
  Object.entries(pairCounts).forEach(([pair, count]) => {
    if (count >= 2) {
      const [c1, c2] = pair.split('|');
      if (components.includes(c1) && components.includes(c2)) {
        interactions.push([c1, c2]);
      }
    }
  });

  return interactions;
}

/**
 * Generate all treatment combinations seen in the data
 */
function generateTreatmentCombinations(components, studies) {
  const seen = new Set();
  const treatments = [];

  // Add reference (no components)
  treatments.push({ label: 'Reference', components: [] });
  seen.add('');

  // Add individual components
  components.forEach(c => {
    const key = c;
    if (!seen.has(key)) {
      treatments.push({ label: c, components: [c] });
      seen.add(key);
    }
  });

  // Add combinations from studies
  studies.forEach(s => {
    [s.t1, s.t2].forEach(t => {
      if (t && t.length > 0) {
        const sorted = [...t].sort();
        const key = sorted.join('+');
        if (!seen.has(key)) {
          treatments.push({ label: key, components: sorted });
          seen.add(key);
        }
      }
    });
  });

  return treatments.sort((a, b) => {
    // Sort by number of components, then alphabetically
    if (a.components.length !== b.components.length) {
      return a.components.length - b.components.length;
    }
    return a.label.localeCompare(b.label);
  });
}

/**
 * CNMA-UpSet plot data preparation
 * For visualizing component combinations
 */
export function prepareCNMAUpsetData(studies, components) {
  const combinations = new Map();

  studies.forEach(s => {
    [{ arm: 't1', components: s.t1 }, { arm: 't2', components: s.t2 }].forEach(({ arm, components: comps }) => {
      if (!comps || comps.length === 0) return;

      const key = [...comps].sort().join('|');
      if (!combinations.has(key)) {
        combinations.set(key, {
          components: comps,
          studies: [],
          count: 0
        });
      }
      combinations.get(key).studies.push(s.id);
      combinations.get(key).count++;
    });
  });

  // Convert to array and sort by frequency
  const data = [...combinations.values()]
    .map(c => ({
      ...c,
      set: c.components.map(comp => components.indexOf(comp)),
      label: c.components.join(' + ')
    }))
    .sort((a, b) => b.count - a.count);

  return {
    components,
    combinations: data,
    matrix: data.map(d => components.map(c => d.components.includes(c) ? 1 : 0))
  };
}

// ============ HELPER FUNCTIONS ============

function estimateRank(matrix) {
  const m = matrix.length;
  if (m === 0) return 0;
  const n = matrix[0].length;

  // Copy matrix
  const A = matrix.map(row => [...row]);

  let rank = 0;
  const rowUsed = Array(m).fill(false);

  for (let j = 0; j < n; j++) {
    // Find pivot
    let pivotRow = -1;
    for (let i = 0; i < m; i++) {
      if (!rowUsed[i] && Math.abs(A[i][j]) > 1e-10) {
        pivotRow = i;
        break;
      }
    }

    if (pivotRow === -1) continue;

    rowUsed[pivotRow] = true;
    rank++;

    // Eliminate
    for (let i = 0; i < m; i++) {
      if (i !== pivotRow && Math.abs(A[i][j]) > 1e-10) {
        const factor = A[i][j] / A[pivotRow][j];
        for (let k = j; k < n; k++) {
          A[i][k] -= factor * A[pivotRow][k];
        }
      }
    }
  }

  return rank;
}

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

function normalQuantile(p) {
  if (p <= 0) return -Infinity;
  if (p >= 1) return Infinity;
  if (p === 0.5) return 0;

  const a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
             1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00];
  const b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
             6.680131188771972e+01, -1.328068155288572e+01];
  const c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
             -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00];
  const d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
             3.754408661907416e+00];

  const pLow = 0.02425, pHigh = 1 - pLow;
  let q, r;

  if (p < pLow) {
    q = Math.sqrt(-2 * Math.log(p));
    return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
  } else if (p <= pHigh) {
    q = p - 0.5;
    r = q * q;
    return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1);
  } else {
    q = Math.sqrt(-2 * Math.log(1 - p));
    return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1);
  }
}
