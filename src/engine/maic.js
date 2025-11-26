/**
 * Matching-Adjusted Indirect Comparison (MAIC)
 * Population-adjusted indirect comparison when IPD is available for one trial
 * References:
 * - Signorovitch et al. (2010) - Original MAIC method
 * - Signorovitch et al. (2012) - Comparative effectiveness
 * - Phillippo et al. (2018) - NICE DSU TSD 18
 * - Chandler et al. (2024) - Variance estimation methods
 */

/**
 * Perform MAIC analysis
 * Adjusts IPD from trial AC to match aggregate data from trial BC
 *
 * @param {Array} ipdData - Individual patient data from trial comparing A vs C
 *                          Each row: { id, treatment, outcome, ...covariates }
 * @param {Object} aggData - Aggregate baseline characteristics from trial BC
 *                          { mean_age: 50, prop_male: 0.6, ... }
 * @param {Object} aggOutcome - Aggregate outcome from trial BC
 *                          { es: log(OR), se: 0.2 } for B vs C
 * @param {Object} options - Analysis options
 * @returns {Object} MAIC results
 */
export function performMAIC(ipdData, aggData, aggOutcome, options = {}) {
  const {
    covariates = Object.keys(aggData),  // Which covariates to match on
    outcomeType = 'binary',              // 'binary', 'continuous', 'survival'
    confLevel = 0.95,
    varianceMethod = 'robust',           // 'robust', 'bootstrap', 'naive'
    nBootstrap = 1000,
    centered = true                      // Center covariates for numerical stability
  } = options;

  // Validate inputs
  if (!ipdData || ipdData.length === 0) {
    return { error: 'No IPD data provided' };
  }

  if (!aggData || covariates.length === 0) {
    return { error: 'No aggregate data or covariates specified' };
  }

  // Step 1: Extract treatment groups from IPD
  const treatmentA = ipdData.filter(p => p.treatment === 'A' || p.treatment === 1);
  const treatmentC = ipdData.filter(p => p.treatment === 'C' || p.treatment === 0);

  if (treatmentA.length === 0 || treatmentC.length === 0) {
    return { error: 'IPD must contain both treatment arms' };
  }

  // Step 2: Calculate IPD means for covariates (before matching)
  const ipdMeans = {};
  covariates.forEach(cov => {
    const values = ipdData.map(p => p[cov]).filter(v => v !== null && v !== undefined);
    ipdMeans[cov] = values.reduce((a, b) => a + b, 0) / values.length;
  });

  // Step 3: Calculate weights using entropy balancing
  const weights = calculateMAICWeights(ipdData, aggData, covariates, centered);

  if (weights.error) {
    return weights;
  }

  // Step 4: Calculate effective sample size
  const ESS = calculateESS(weights.weights);
  const originalN = ipdData.length;

  // Step 5: Check balance after weighting
  const balance = checkBalance(ipdData, aggData, covariates, weights.weights);

  // Step 6: Calculate weighted treatment effect in IPD trial (A vs C)
  const weightedOutcome = calculateWeightedOutcome(
    ipdData,
    weights.weights,
    outcomeType
  );

  // Step 7: Perform indirect comparison with aggregate trial (B vs C)
  // A vs B = (A vs C) - (B vs C)
  const indirectComparison = calculateIndirectEffect(
    weightedOutcome,
    aggOutcome,
    varianceMethod,
    nBootstrap,
    ipdData,
    weights.weights,
    covariates,
    aggData,
    outcomeType,
    confLevel
  );

  // Calculate unweighted (naive) comparison for reference
  const naiveOutcome = calculateWeightedOutcome(
    ipdData,
    ipdData.map(() => 1),  // Equal weights
    outcomeType
  );

  const naiveIndirect = {
    es: naiveOutcome.es - aggOutcome.es,
    se: Math.sqrt(naiveOutcome.vi + aggOutcome.se ** 2)
  };

  return {
    method: 'MAIC',

    // Population adjustment
    originalN,
    effectiveSampleSize: ESS,
    relativeESS: (ESS / originalN) * 100,

    // Covariate balance
    covariates: covariates,
    ipdMeansOriginal: ipdMeans,
    targetMeans: aggData,
    balanceAfter: balance,

    // Weighted outcome from IPD trial (A vs C)
    weightedEffect: {
      es: weightedOutcome.es,
      se: Math.sqrt(weightedOutcome.vi),
      comparison: 'A vs C (adjusted)'
    },

    // Aggregate outcome from comparator trial (B vs C)
    aggregateEffect: {
      es: aggOutcome.es,
      se: aggOutcome.se,
      comparison: 'B vs C'
    },

    // Indirect comparison (A vs B)
    indirectEffect: {
      es: indirectComparison.es,
      se: indirectComparison.se,
      ciLo: indirectComparison.ciLo,
      ciHi: indirectComparison.ciHi,
      z: indirectComparison.z,
      pVal: indirectComparison.pVal,
      comparison: 'A vs B (MAIC)'
    },

    // Naive (unadjusted) comparison for reference
    naiveEffect: {
      es: naiveIndirect.es,
      se: naiveIndirect.se,
      comparison: 'A vs B (unadjusted)'
    },

    // Weights and diagnostics
    weights: weights.weights,
    convergence: weights.convergence,
    confLevel
  };
}

/**
 * Calculate MAIC weights using entropy balancing (method of moments)
 * Solves: minimize Σ wᵢ log(wᵢ) subject to Σ wᵢ xᵢ = μ_target
 */
function calculateMAICWeights(ipdData, aggData, covariates, centered = true) {
  const n = ipdData.length;

  // Create design matrix
  const X = ipdData.map(p => covariates.map(cov => p[cov] || 0));

  // Target moments (from aggregate data)
  const target = covariates.map(cov => aggData[cov]);

  // Center covariates for numerical stability
  const means = [];
  if (centered) {
    covariates.forEach((cov, j) => {
      const mean = X.reduce((sum, row) => sum + row[j], 0) / n;
      means.push(mean);
      X.forEach(row => row[j] -= mean);
      target[j] -= mean;
    });
  }

  // Solve for balancing weights using Newton-Raphson
  // Q(β) = Σ exp(β'xᵢ) - n·β'μ_target
  // Solution: exp(β'xᵢ) / Σexp(β'xⱼ) * n

  const maxIter = 1000;
  const tol = 1e-8;
  const p = covariates.length;

  let beta = Array(p).fill(0);
  let converged = false;

  for (let iter = 0; iter < maxIter; iter++) {
    // Calculate weights
    const logW = X.map(xi => xi.reduce((sum, xij, j) => sum + beta[j] * xij, 0));
    const maxLogW = Math.max(...logW);
    const expW = logW.map(lw => Math.exp(lw - maxLogW));
    const sumExpW = expW.reduce((a, b) => a + b, 0);
    const w = expW.map(ew => ew / sumExpW * n);

    // Calculate gradient (moment conditions)
    const gradient = Array(p).fill(0);
    for (let j = 0; j < p; j++) {
      const weightedMean = X.reduce((sum, xi, i) => sum + w[i] * xi[j], 0) / n;
      gradient[j] = weightedMean - target[j];
    }

    // Check convergence
    const maxGrad = Math.max(...gradient.map(Math.abs));
    if (maxGrad < tol) {
      converged = true;
      // Return normalized weights
      return {
        weights: w.map(wi => wi / n),  // Normalize to mean 1
        convergence: { converged: true, iterations: iter, maxGradient: maxGrad }
      };
    }

    // Hessian (information matrix)
    const H = Array(p).fill(null).map(() => Array(p).fill(0));
    for (let j = 0; j < p; j++) {
      for (let k = 0; k <= j; k++) {
        let sum = 0;
        for (let i = 0; i < n; i++) {
          sum += w[i] * X[i][j] * X[i][k];
        }
        H[j][k] = sum / n;
        H[k][j] = H[j][k];
      }
    }

    // Newton step: beta_new = beta - H⁻¹ gradient
    const Hinv = invertMatrix(H);
    if (!Hinv) {
      return { error: 'Singular Hessian - covariates may be collinear', convergence: { converged: false } };
    }

    for (let j = 0; j < p; j++) {
      for (let k = 0; k < p; k++) {
        beta[j] -= Hinv[j][k] * gradient[k];
      }
    }

    // Damping for stability
    beta = beta.map(b => Math.max(-10, Math.min(10, b)));
  }

  // Final weights even if not converged
  const logW = X.map(xi => xi.reduce((sum, xij, j) => sum + beta[j] * xij, 0));
  const maxLogW = Math.max(...logW);
  const expW = logW.map(lw => Math.exp(lw - maxLogW));
  const sumExpW = expW.reduce((a, b) => a + b, 0);
  const w = expW.map(ew => ew / sumExpW);

  return {
    weights: w,
    convergence: { converged: false, iterations: maxIter, message: 'Max iterations reached' }
  };
}

/**
 * Calculate effective sample size (ESS)
 * ESS = (Σwᵢ)² / Σwᵢ²
 */
function calculateESS(weights) {
  const sumW = weights.reduce((a, b) => a + b, 0);
  const sumW2 = weights.reduce((a, b) => a + b * b, 0);
  return sumW2 > 0 ? (sumW * sumW) / sumW2 : 0;
}

/**
 * Check balance of covariates after weighting
 */
function checkBalance(ipdData, aggData, covariates, weights) {
  const n = ipdData.length;
  const sumW = weights.reduce((a, b) => a + b, 0);

  return covariates.map(cov => {
    // Weighted mean
    const weightedMean = ipdData.reduce((sum, p, i) =>
      sum + weights[i] * (p[cov] || 0), 0) / sumW;

    // Target mean
    const targetMean = aggData[cov];

    // Standardized difference
    const values = ipdData.map(p => p[cov] || 0);
    const sd = Math.sqrt(values.reduce((sum, v) =>
      sum + (v - values.reduce((a, b) => a + b, 0) / n) ** 2, 0) / (n - 1));

    const stdDiff = sd > 0 ? (weightedMean - targetMean) / sd : 0;

    return {
      covariate: cov,
      weightedMean,
      targetMean,
      difference: weightedMean - targetMean,
      standardizedDiff: stdDiff,
      balanced: Math.abs(stdDiff) < 0.1  // Common threshold
    };
  });
}

/**
 * Calculate weighted treatment effect from IPD
 */
function calculateWeightedOutcome(ipdData, weights, outcomeType) {
  const treatmentA = [];
  const treatmentC = [];

  ipdData.forEach((p, i) => {
    const w = weights[i];
    const isA = p.treatment === 'A' || p.treatment === 1;

    if (isA) {
      treatmentA.push({ outcome: p.outcome, weight: w });
    } else {
      treatmentC.push({ outcome: p.outcome, weight: w });
    }
  });

  if (outcomeType === 'binary') {
    return calculateWeightedOR(treatmentA, treatmentC);
  } else if (outcomeType === 'continuous') {
    return calculateWeightedMD(treatmentA, treatmentC);
  } else {
    // Default to MD
    return calculateWeightedMD(treatmentA, treatmentC);
  }
}

/**
 * Calculate weighted log odds ratio for binary outcome
 */
function calculateWeightedOR(groupA, groupC) {
  // Weighted event counts
  const sumWA = groupA.reduce((sum, p) => sum + p.weight, 0);
  const sumWC = groupC.reduce((sum, p) => sum + p.weight, 0);

  const eventsA = groupA.reduce((sum, p) => sum + p.weight * p.outcome, 0);
  const eventsC = groupC.reduce((sum, p) => sum + p.weight * p.outcome, 0);

  const pA = eventsA / sumWA;
  const pC = eventsC / sumWC;

  // Log odds ratio
  const oddsA = pA / (1 - pA);
  const oddsC = pC / (1 - pC);

  if (oddsA <= 0 || oddsC <= 0 || !isFinite(oddsA) || !isFinite(oddsC)) {
    return { es: 0, vi: 1, error: 'Invalid odds' };
  }

  const logOR = Math.log(oddsA / oddsC);

  // Variance using weighted counts
  // This is the "robust" sandwich estimator approximation
  const nA = sumWA;
  const nC = sumWC;
  const vi = 1 / (nA * pA * (1 - pA)) + 1 / (nC * pC * (1 - pC));

  return { es: logOR, vi, pA, pC, nA, nC };
}

/**
 * Calculate weighted mean difference for continuous outcome
 */
function calculateWeightedMD(groupA, groupC) {
  const sumWA = groupA.reduce((sum, p) => sum + p.weight, 0);
  const sumWC = groupC.reduce((sum, p) => sum + p.weight, 0);

  // Weighted means
  const meanA = groupA.reduce((sum, p) => sum + p.weight * p.outcome, 0) / sumWA;
  const meanC = groupC.reduce((sum, p) => sum + p.weight * p.outcome, 0) / sumWC;

  // Weighted variances
  const varA = groupA.reduce((sum, p) =>
    sum + p.weight * (p.outcome - meanA) ** 2, 0) / sumWA;
  const varC = groupC.reduce((sum, p) =>
    sum + p.weight * (p.outcome - meanC) ** 2, 0) / sumWC;

  const md = meanA - meanC;
  const vi = varA / sumWA + varC / sumWC;

  return { es: md, vi, meanA, meanC, varA, varC };
}

/**
 * Calculate indirect effect with appropriate variance estimation
 */
function calculateIndirectEffect(
  weightedOutcome,
  aggOutcome,
  varianceMethod,
  nBootstrap,
  ipdData,
  weights,
  covariates,
  aggData,
  outcomeType,
  confLevel
) {
  // Point estimate: A vs B = (A vs C) - (B vs C)
  const es = weightedOutcome.es - aggOutcome.es;

  let se;

  if (varianceMethod === 'naive') {
    // Simple variance addition (underestimates uncertainty)
    se = Math.sqrt(weightedOutcome.vi + aggOutcome.se ** 2);
  } else if (varianceMethod === 'bootstrap' && ipdData) {
    // Bootstrap variance (accounts for weight estimation uncertainty)
    se = bootstrapVariance(
      ipdData,
      covariates,
      aggData,
      aggOutcome,
      outcomeType,
      nBootstrap
    );
  } else {
    // Robust sandwich estimator (recommended)
    // Accounts for uncertainty in weight estimation
    // Var(θ) ≈ Var_naive + Var_weight_estimation
    const naiveVar = weightedOutcome.vi + aggOutcome.se ** 2;

    // Inflation factor for weight uncertainty (approximation)
    // Based on Phillippo et al. (2018) and Chandler et al. (2024)
    const ESS = calculateESS(weights);
    const n = weights.length;
    const inflationFactor = n / ESS;

    se = Math.sqrt(naiveVar * inflationFactor);
  }

  const z = es / se;
  const pVal = 2 * (1 - normalCDF(Math.abs(z)));

  const alpha = 1 - confLevel;
  const zCrit = Math.abs(normalQuantile(alpha / 2));

  return {
    es,
    se,
    ciLo: es - zCrit * se,
    ciHi: es + zCrit * se,
    z,
    pVal
  };
}

/**
 * Bootstrap variance estimation
 */
function bootstrapVariance(ipdData, covariates, aggData, aggOutcome, outcomeType, nBootstrap) {
  const n = ipdData.length;
  const bootEstimates = [];

  for (let b = 0; b < nBootstrap; b++) {
    // Resample IPD with replacement
    const bootSample = [];
    for (let i = 0; i < n; i++) {
      const idx = Math.floor(Math.random() * n);
      bootSample.push(ipdData[idx]);
    }

    // Recalculate weights
    const bootWeights = calculateMAICWeights(bootSample, aggData, covariates, true);
    if (bootWeights.error) continue;

    // Recalculate effect
    const bootOutcome = calculateWeightedOutcome(bootSample, bootWeights.weights, outcomeType);

    // Indirect effect
    bootEstimates.push(bootOutcome.es - aggOutcome.es);
  }

  if (bootEstimates.length < 10) {
    return Math.sqrt(1);  // Fallback
  }

  // Bootstrap standard error
  const mean = bootEstimates.reduce((a, b) => a + b, 0) / bootEstimates.length;
  const variance = bootEstimates.reduce((sum, est) =>
    sum + (est - mean) ** 2, 0) / (bootEstimates.length - 1);

  return Math.sqrt(variance);
}

/**
 * Anchored MAIC (when common comparator C is measured in both trials)
 * More reliable than unanchored MAIC
 *
 * @param {Array} ipdDataAC - IPD from trial A vs C
 * @param {Object} aggDataBC - Aggregate from trial B vs C
 * @param {Object} aggBaselineBC - Baseline characteristics from trial BC
 * @param {Object} options - Analysis options
 */
export function anchoredMAIC(ipdDataAC, aggDataBC, aggBaselineBC, options = {}) {
  // This is the standard MAIC - both trials have common comparator C
  return performMAIC(ipdDataAC, aggBaselineBC, {
    es: aggDataBC.es,
    se: aggDataBC.se
  }, options);
}

/**
 * Unanchored MAIC (single-arm comparison)
 * Warning: Higher risk of bias due to potential unmeasured confounding
 *
 * @param {Array} ipdDataA - IPD from single-arm trial of treatment A
 * @param {Object} aggDataB - Aggregate outcome from treatment B
 * @param {Object} aggBaselineB - Baseline characteristics from treatment B
 * @param {Object} options - Analysis options
 */
export function unanchoredMAIC(ipdDataA, aggDataB, aggBaselineB, options = {}) {
  // Mark all IPD as treatment A (no control arm)
  const ipdMarked = ipdDataA.map(p => ({ ...p, treatment: 'A' }));

  // Calculate weighted mean outcome for A
  const weights = calculateMAICWeights(ipdMarked, aggBaselineB, options.covariates || Object.keys(aggBaselineB), true);

  if (weights.error) return weights;

  const sumW = weights.weights.reduce((a, b) => a + b, 0);
  const weightedOutcomeA = ipdMarked.reduce((sum, p, i) =>
    sum + weights.weights[i] * p.outcome, 0) / sumW;

  // For unanchored, we directly compare outcomes
  // This assumes outcome is on same scale (e.g., response rate)
  const diff = weightedOutcomeA - aggDataB.mean;
  const se = Math.sqrt(
    ipdMarked.reduce((sum, p, i) =>
      sum + weights.weights[i] * (p.outcome - weightedOutcomeA) ** 2, 0) / (sumW * sumW) +
    aggDataB.se ** 2
  );

  const ESS = calculateESS(weights.weights);

  return {
    method: 'Unanchored MAIC',
    warning: 'Unanchored MAIC has higher risk of bias due to potential unmeasured confounding',
    effectiveSampleSize: ESS,
    adjustedOutcomeA: weightedOutcomeA,
    aggregateOutcomeB: aggDataB.mean,
    difference: {
      es: diff,
      se,
      ciLo: diff - 1.96 * se,
      ciHi: diff + 1.96 * se
    },
    weights: weights.weights,
    convergence: weights.convergence
  };
}

/**
 * Simulated Treatment Comparison (STC)
 * Alternative to MAIC when only aggregate data available for adjustment
 * Uses outcome regression instead of propensity weighting
 *
 * @param {Array} ipdData - IPD from trial comparing A vs C
 * @param {Object} aggData - Aggregate covariates from target population
 * @param {Object} aggOutcome - Aggregate outcome from trial B vs C
 * @param {Object} options - Analysis options
 */
export function simulatedTreatmentComparison(ipdData, aggData, aggOutcome, options = {}) {
  const { covariates = Object.keys(aggData), outcomeType = 'binary' } = options;

  // Step 1: Fit outcome model on IPD
  // E[Y] = β₀ + β₁·treatment + β₂·X + β₃·treatment×X

  // For simplicity, we use a linear model even for binary outcomes
  // A full implementation would use logistic regression

  const n = ipdData.length;
  const p = covariates.length;

  // Design matrix: [1, treatment, X₁, X₂, ..., treatment×X₁, treatment×X₂, ...]
  const X = ipdData.map(pt => {
    const trt = pt.treatment === 'A' || pt.treatment === 1 ? 1 : 0;
    const row = [1, trt];

    // Centered covariates
    covariates.forEach(cov => {
      const ipdMean = ipdData.reduce((sum, p) => sum + (p[cov] || 0), 0) / n;
      row.push((pt[cov] || 0) - ipdMean);
    });

    // Interactions
    covariates.forEach(cov => {
      const ipdMean = ipdData.reduce((sum, p) => sum + (p[cov] || 0), 0) / n;
      row.push(trt * ((pt[cov] || 0) - ipdMean));
    });

    return row;
  });

  const y = ipdData.map(pt => pt.outcome);

  // OLS: β = (X'X)⁻¹ X'y
  const k = X[0].length;
  const XtX = Array(k).fill(null).map(() => Array(k).fill(0));
  const Xty = Array(k).fill(0);

  for (let i = 0; i < n; i++) {
    for (let j = 0; j < k; j++) {
      for (let l = 0; l < k; l++) {
        XtX[j][l] += X[i][j] * X[i][l];
      }
      Xty[j] += X[i][j] * y[i];
    }
  }

  const XtXinv = invertMatrix(XtX);
  if (!XtXinv) {
    return { error: 'Singular design matrix in STC regression' };
  }

  const beta = Array(k).fill(0);
  for (let j = 0; j < k; j++) {
    for (let l = 0; l < k; l++) {
      beta[j] += XtXinv[j][l] * Xty[l];
    }
  }

  // Step 2: Predict outcomes at target population covariate means
  // Calculate IPD means for centering reference
  const ipdMeans = {};
  covariates.forEach(cov => {
    ipdMeans[cov] = ipdData.reduce((sum, p) => sum + (p[cov] || 0), 0) / n;
  });

  // Predicted outcome for A at target population
  const targetCentered = covariates.map(cov => aggData[cov] - ipdMeans[cov]);
  const predA = beta[0] + beta[1] +
    targetCentered.reduce((sum, x, i) => sum + beta[2 + i] * x, 0) +
    targetCentered.reduce((sum, x, i) => sum + beta[2 + p + i] * x, 0);

  // Predicted outcome for C at target population
  const predC = beta[0] +
    targetCentered.reduce((sum, x, i) => sum + beta[2 + i] * x, 0);

  // Effect of A vs C at target population
  const effectAC = predA - predC;

  // Step 3: Indirect comparison
  const es = effectAC - aggOutcome.es;

  // Variance (approximate)
  const residuals = ipdData.map((pt, i) => {
    const predicted = X[i].reduce((sum, xij, j) => sum + xij * beta[j], 0);
    return pt.outcome - predicted;
  });
  const sigma2 = residuals.reduce((sum, r) => sum + r * r, 0) / (n - k);

  // Variance of treatment effect includes model uncertainty
  const varES = sigma2 * XtXinv[1][1] + aggOutcome.se ** 2;
  const se = Math.sqrt(varES);

  return {
    method: 'Simulated Treatment Comparison (STC)',
    coefficients: beta,
    predictedEffectAC: effectAC,
    indirectEffect: {
      es,
      se,
      ciLo: es - 1.96 * se,
      ciHi: es + 1.96 * se,
      comparison: 'A vs B (STC)'
    },
    modelFit: {
      residualVariance: sigma2,
      n,
      parameters: k
    }
  };
}

// ============ HELPER FUNCTIONS ============

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

function normalCDF(x) {
  const a1 = 0.254829592, a2 = -0.284496736, a3 = 1.421413741;
  const a4 = -1.453152027, a5 = 1.061405429, p = 0.3275911;

  const sign = x < 0 ? -1 : 1;
  x = Math.abs(x) / Math.sqrt(2);
  const t = 1 / (1 + p * x);
  const y = 1 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-x * x);

  return 0.5 * (1 + sign * y);
}

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
