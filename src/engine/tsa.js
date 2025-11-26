/**
 * Trial Sequential Analysis (TSA) calculations
 */

import { tauDL } from './tau.js';
import { normQuantile } from './distributions.js';

/**
 * Calculate Trial Sequential Analysis
 * @param {Array} effects - Array of effect objects
 * @param {Object} params - TSA parameters
 * @returns {Object|null} TSA results
 */
export function calculateTSA(effects, params = {}) {
  const {
    alpha = 0.05,
    beta = 0.20,
    rrr = 0.20, // Relative risk reduction
    cer = 0.10  // Control event rate
  } = params;

  const active = effects.filter(e => !e.excluded);
  if (active.length < 2) return null;

  // Calculate diversity (D²) - adjusted heterogeneity
  const tau2 = tauDL(active);
  const avgVi = active.reduce((a, e) => a + e.vi, 0) / active.length;
  const D2 = tau2 / (tau2 + avgVi) * 100;

  // Z-values for boundaries
  const zAlpha = normQuantile(1 - alpha / 2); // Two-sided
  const zBeta = normQuantile(1 - beta);

  // Calculate Required Information Size (RIS)
  const RIS = calculateRIS(cer, rrr, zAlpha, zBeta, D2);

  // Accrued information (total participants)
  const accrued = active.reduce((acc, s) => {
    if (s.raw) {
      return acc + (s.raw.n1 || 0) + (s.raw.n2 || 0);
    }
    return acc;
  }, 0);

  // Calculate cumulative statistics using RE weights
  const cumulative = calculateCumulativeZ(active, tau2);

  // Calculate monitoring and futility boundaries
  const boundaries = calculateBoundaries(cumulative, RIS, zAlpha, zBeta);

  // Determine TSA conclusion
  const conclusion = determineTSAConclusion(cumulative, boundaries, RIS);

  return {
    RIS,
    accrued,
    pctRIS: ((accrued / RIS) * 100).toFixed(1),
    D2: D2.toFixed(1),
    tau2,
    cumZ: cumulative.z,
    cumN: cumulative.n,
    infoFracs: cumulative.infoFrac,
    monitoringBoundary: boundaries.monitoring,
    futilityBoundary: boundaries.futility,
    conclusion: conclusion.verdict,
    boundary: conclusion.description,
    crossedBenefit: conclusion.crossedBenefit,
    crossedHarm: conclusion.crossedHarm,
    crossedFutility: conclusion.crossedFutility
  };
}

/**
 * Calculate Required Information Size
 * @param {number} cer - Control event rate
 * @param {number} rrr - Relative risk reduction
 * @param {number} zAlpha - Z for alpha
 * @param {number} zBeta - Z for beta
 * @param {number} D2 - Diversity percentage
 * @returns {number} Required information size
 */
function calculateRIS(cer, rrr, zAlpha, zBeta, D2) {
  const pC = cer;
  const pE = cer * (1 - rrr);
  const pooledP = (pC + pE) / 2;

  // Pooled variance for binary outcome
  const variance = 2 * pooledP * (1 - pooledP);

  // Effect size squared
  const delta2 = (pC - pE) ** 2;

  // Sample size per group for single trial (equal allocation)
  const nPerGroup = variance * (zAlpha + zBeta) ** 2 / delta2;

  // Total sample size adjusted for diversity
  const diversityFactor = 1 - D2 / 100;
  const RIS = Math.ceil((2 * nPerGroup) / Math.max(0.01, diversityFactor));

  return Math.max(RIS, 100); // Minimum RIS
}

/**
 * Calculate cumulative Z-scores using random-effects weights
 * @param {Array} active - Active effect objects
 * @param {number} tau2 - Between-study variance estimate
 * @returns {Object} Cumulative statistics
 */
function calculateCumulativeZ(active, tau2 = 0) {
  const z = [];
  const n = [];
  const es = [];
  const infoFrac = [];

  let cumES = 0;
  let cumW = 0;

  active.forEach((s, i) => {
    // Use random-effects weights (vi + tau²)
    const w = 1 / (s.vi + tau2);
    cumW += w;
    cumES = (cumES * (cumW - w) + w * s.es) / cumW;
    const cumSE = Math.sqrt(1 / cumW);

    z.push(cumES / cumSE);
    es.push(cumES);

    const cumulativeN = active.slice(0, i + 1).reduce((acc, x) => {
      if (x.raw) {
        return acc + (x.raw.n1 || 0) + (x.raw.n2 || 0);
      }
      return acc + 100; // Default if no raw data
    }, 0);

    n.push(cumulativeN);
  });

  return { z, n, es, infoFrac };
}

/**
 * Calculate O'Brien-Fleming monitoring boundaries
 * @param {Object} cumulative - Cumulative statistics
 * @param {number} RIS - Required information size
 * @param {number} zAlpha - Z for alpha
 * @param {number} zBeta - Z for beta
 * @returns {Object} Monitoring and futility boundaries
 */
function calculateBoundaries(cumulative, RIS, zAlpha, zBeta) {
  // Information fractions
  const infoFracs = cumulative.n.map(n => n / RIS);

  // O'Brien-Fleming spending function for monitoring boundary
  const monitoring = infoFracs.map(t => {
    if (t <= 0) return 8; // Very large boundary at start
    return zAlpha / Math.sqrt(t);
  });

  // Futility boundary (beta-spending)
  const futility = infoFracs.map(t => {
    if (t <= 0) return 0;
    // Simplified inner wedge
    return -0.5 + zBeta * Math.sqrt(t);
  });

  return { monitoring, futility, infoFracs };
}

/**
 * Determine TSA conclusion based on boundaries
 * @param {Object} cumulative - Cumulative statistics
 * @param {Object} boundaries - Calculated boundaries
 * @param {number} RIS - Required information size
 * @returns {Object} Conclusion
 */
function determineTSAConclusion(cumulative, boundaries, RIS) {
  const lastIdx = cumulative.z.length - 1;
  const lastZ = cumulative.z[lastIdx];
  const lastN = cumulative.n[lastIdx];
  const lastMonitor = boundaries.monitoring[lastIdx];
  const lastFutility = boundaries.futility[lastIdx];

  let verdict = 'Inconclusive';
  let description = 'Continue accrual';
  let crossedBenefit = false;
  let crossedHarm = false;
  let crossedFutility = false;

  // Check if crossed monitoring boundary
  if (lastZ >= lastMonitor) {
    verdict = 'Benefit';
    description = 'Crossed monitoring boundary (benefit)';
    crossedBenefit = true;
  } else if (lastZ <= -lastMonitor) {
    verdict = 'Harm';
    description = 'Crossed monitoring boundary (harm)';
    crossedHarm = true;
  } else if (Math.abs(lastZ) <= Math.abs(lastFutility) && lastN >= RIS * 0.5) {
    verdict = 'Futility';
    description = 'Crossed futility boundary';
    crossedFutility = true;
  } else if (lastN >= RIS) {
    verdict = 'RIS Reached';
    description = 'Required information size reached';
    if (lastZ > 1.96) {
      verdict = 'Benefit (RIS)';
      crossedBenefit = true;
    } else if (lastZ < -1.96) {
      verdict = 'Harm (RIS)';
      crossedHarm = true;
    } else {
      verdict = 'No effect (RIS)';
      crossedFutility = true;
    }
  } else if (lastN >= RIS * 0.9) {
    verdict = 'Near RIS';
    description = 'Near RIS without crossing boundaries';
  }

  return { verdict, description, crossedBenefit, crossedHarm, crossedFutility };
}

/**
 * Generate TSA plot data for visualization
 * @param {Object} tsaResult - TSA calculation result
 * @returns {Object} Plot data
 */
export function getTSAPlotData(tsaResult) {
  if (!tsaResult) return null;

  const maxN = Math.max(tsaResult.RIS, Math.max(...tsaResult.cumN)) * 1.1;

  return {
    cumN: tsaResult.cumN,
    cumZ: tsaResult.cumZ,
    monitoringUpper: tsaResult.monitoringBoundary,
    monitoringLower: tsaResult.monitoringBoundary.map(b => -b),
    futilityUpper: tsaResult.futilityBoundary,
    futilityLower: tsaResult.futilityBoundary.map(b => -b),
    RIS: tsaResult.RIS,
    maxN,
    traditionalZ: 1.96
  };
}
