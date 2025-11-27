/**
 * Transportability and Generalizability Analysis
 *
 * Methods for assessing and adjusting trial results for target populations:
 * - Inverse odds weighting (IOW) for transportability
 * - Inverse probability of selection weighting (IPSW) for generalizability
 * - Calibration weighting
 * - Population comparison and diagnostics
 *
 * Based on LFA (Local-to-Federated Analysis) and FDA guidance on transportability
 */

import { jStat } from 'jstat';

/**
 * Calculate population similarity metrics
 * Compares trial population characteristics to target population
 *
 * @param {Array} trialChars - Trial population characteristics [{name, mean, sd}]
 * @param {Array} targetChars - Target population characteristics [{name, mean, sd}]
 * @returns {Object} Similarity metrics
 */
export function calculatePopulationSimilarity(trialChars, targetChars) {
  if (!trialChars.length || !targetChars.length) {
    return { error: 'Population characteristics required' };
  }

  const metrics = [];
  let totalStdDiff = 0;
  let nComparisons = 0;

  for (const trialChar of trialChars) {
    const targetChar = targetChars.find(t => t.name === trialChar.name);
    if (!targetChar) continue;

    // Standardized mean difference (SMD)
    const pooledSD = Math.sqrt(
      (trialChar.sd * trialChar.sd + targetChar.sd * targetChar.sd) / 2
    );

    const smd = pooledSD > 0
      ? (trialChar.mean - targetChar.mean) / pooledSD
      : 0;

    // Variance ratio
    const varianceRatio = targetChar.sd > 0
      ? (trialChar.sd * trialChar.sd) / (targetChar.sd * targetChar.sd)
      : 1;

    metrics.push({
      name: trialChar.name,
      trialMean: trialChar.mean,
      targetMean: targetChar.mean,
      smd: smd,
      absSmd: Math.abs(smd),
      varianceRatio,
      imbalanced: Math.abs(smd) > 0.1  // Common threshold
    });

    totalStdDiff += Math.abs(smd);
    nComparisons++;
  }

  // Overall similarity index
  const avgAbsSmd = nComparisons > 0 ? totalStdDiff / nComparisons : 0;

  // Mahalanobis-style distance (simplified)
  const mahalanobisApprox = Math.sqrt(
    metrics.reduce((sum, m) => sum + m.smd * m.smd, 0)
  );

  return {
    characteristics: metrics,
    nComparisons,
    avgAbsSmd,
    mahalanobisDistance: mahalanobisApprox,
    overallSimilarity: avgAbsSmd < 0.1 ? 'high' : avgAbsSmd < 0.25 ? 'moderate' : 'low',
    imbalancedCount: metrics.filter(m => m.imbalanced).length,
    recommendation: getTransportabilityRecommendation(avgAbsSmd, metrics)
  };
}

/**
 * Inverse Odds Weighting (IOW) for Transportability
 * Weights trial participants to match target population
 *
 * @param {Array} trialData - Individual participant data [{id, treated, outcome, covariates: {...}}]
 * @param {Array} targetSummary - Target population summary [{name, mean, sd}]
 * @param {Object} options - Weighting options
 * @returns {Object} Weighted analysis results
 */
export function inverseOddsWeighting(trialData, targetSummary, options = {}) {
  const { covariateNames = [], truncate = 0.01 } = options;

  if (!trialData.length) {
    return { error: 'Trial data required' };
  }

  const n = trialData.length;

  // Use all covariates if none specified
  const covNames = covariateNames.length > 0
    ? covariateNames
    : Object.keys(trialData[0].covariates || {});

  if (covNames.length === 0) {
    return { error: 'Covariates required for weighting' };
  }

  // Calculate propensity scores (simplified logistic model)
  // In practice, this would use actual logistic regression
  // Here we use a simplified distance-based approach

  const weights = trialData.map(participant => {
    let logOdds = 0;

    for (const covName of covNames) {
      const targetCov = targetSummary.find(t => t.name === covName);
      if (!targetCov) continue;

      const participantValue = participant.covariates[covName];
      if (participantValue == null) continue;

      // Standardized distance from target mean
      const zScore = targetCov.sd > 0
        ? (participantValue - targetCov.mean) / targetCov.sd
        : 0;

      // Add to log-odds (participants closer to target mean get higher weight)
      logOdds -= zScore * zScore * 0.5;
    }

    // Convert to probability and then to weight
    const prob = 1 / (1 + Math.exp(-logOdds));
    let weight = prob / (1 - Math.max(prob, 0.001));

    return weight;
  });

  // Normalize weights
  const sumWeights = weights.reduce((a, b) => a + b, 0);
  const normalizedWeights = weights.map(w => w * n / sumWeights);

  // Truncate extreme weights
  const lowerBound = jStat.percentile(normalizedWeights, truncate * 100);
  const upperBound = jStat.percentile(normalizedWeights, (1 - truncate) * 100);

  const truncatedWeights = normalizedWeights.map(w =>
    Math.max(lowerBound, Math.min(upperBound, w))
  );

  // Re-normalize after truncation
  const sumTrunc = truncatedWeights.reduce((a, b) => a + b, 0);
  const finalWeights = truncatedWeights.map(w => w * n / sumTrunc);

  // Calculate weighted outcomes
  const treated = trialData.filter((p, i) => p.treated);
  const control = trialData.filter((p, i) => !p.treated);

  // Weighted mean outcomes
  let sumWtY1 = 0, sumWt1 = 0;
  let sumWtY0 = 0, sumWt0 = 0;

  trialData.forEach((p, i) => {
    if (p.treated) {
      sumWtY1 += finalWeights[i] * p.outcome;
      sumWt1 += finalWeights[i];
    } else {
      sumWtY0 += finalWeights[i] * p.outcome;
      sumWt0 += finalWeights[i];
    }
  });

  const weightedMeanTreated = sumWt1 > 0 ? sumWtY1 / sumWt1 : 0;
  const weightedMeanControl = sumWt0 > 0 ? sumWtY0 / sumWt0 : 0;
  const weightedATE = weightedMeanTreated - weightedMeanControl;

  // Unweighted comparison
  const unweightedTreated = treated.reduce((sum, p) => sum + p.outcome, 0) / treated.length;
  const unweightedControl = control.reduce((sum, p) => sum + p.outcome, 0) / control.length;
  const unweightedATE = unweightedTreated - unweightedControl;

  // Effective sample size
  const ESS = (finalWeights.reduce((sum, w) => sum + w, 0) ** 2) /
    finalWeights.reduce((sum, w) => sum + w * w, 0);

  return {
    method: 'Inverse Odds Weighting',
    n,
    nTreated: treated.length,
    nControl: control.length,
    covariates: covNames,
    unweighted: {
      treatedMean: unweightedTreated,
      controlMean: unweightedControl,
      ATE: unweightedATE
    },
    weighted: {
      treatedMean: weightedMeanTreated,
      controlMean: weightedMeanControl,
      ATE: weightedATE
    },
    effectChange: weightedATE - unweightedATE,
    ESS,
    ESSreduction: 1 - ESS / n,
    weights: {
      mean: finalWeights.reduce((a, b) => a + b, 0) / n,
      min: Math.min(...finalWeights),
      max: Math.max(...finalWeights),
      cv: calculateCV(finalWeights)
    },
    diagnostics: assessWeightingDiagnostics(finalWeights, ESS, n)
  };
}

/**
 * Calculate coefficient of variation
 */
function calculateCV(arr) {
  const mean = arr.reduce((a, b) => a + b, 0) / arr.length;
  const variance = arr.reduce((sum, x) => sum + (x - mean) ** 2, 0) / arr.length;
  return mean > 0 ? Math.sqrt(variance) / mean : 0;
}

/**
 * Assess weighting diagnostics
 */
function assessWeightingDiagnostics(weights, ESS, n) {
  const issues = [];

  // Check ESS reduction
  const essRatio = ESS / n;
  if (essRatio < 0.5) {
    issues.push('ESS reduced by more than 50% - consider model simplification');
  }

  // Check weight extremity
  const cv = calculateCV(weights);
  if (cv > 1) {
    issues.push('High coefficient of variation in weights - positivity may be violated');
  }

  const maxW = Math.max(...weights);
  const meanW = weights.reduce((a, b) => a + b, 0) / weights.length;
  if (maxW > 10 * meanW) {
    issues.push('Some weights are extremely large - check for near-positivity violations');
  }

  return {
    essRatio,
    cv,
    issues,
    acceptable: issues.length === 0
  };
}

/**
 * Transportability Analysis for Meta-Analysis
 * Assesses whether each study's population is transportable to target
 *
 * @param {Array} studies - Studies with population characteristics
 * @param {Array} targetChars - Target population characteristics
 * @returns {Object} Study-level transportability assessment
 */
export function assessStudyTransportability(studies, targetChars) {
  const assessments = [];

  for (const study of studies) {
    if (!study.populationChars) {
      assessments.push({
        id: study.id,
        transportable: 'unknown',
        reason: 'No population characteristics available'
      });
      continue;
    }

    const similarity = calculatePopulationSimilarity(study.populationChars, targetChars);

    if (similarity.error) {
      assessments.push({
        id: study.id,
        transportable: 'unknown',
        reason: similarity.error
      });
      continue;
    }

    assessments.push({
      id: study.id,
      avgAbsSmd: similarity.avgAbsSmd,
      transportable: similarity.overallSimilarity === 'high' ? 'yes' :
        similarity.overallSimilarity === 'moderate' ? 'conditional' : 'limited',
      imbalancedCovariates: similarity.characteristics
        .filter(c => c.imbalanced)
        .map(c => c.name),
      recommendation: similarity.recommendation
    });
  }

  // Summary statistics
  const transportable = assessments.filter(a => a.transportable === 'yes').length;
  const conditional = assessments.filter(a => a.transportable === 'conditional').length;
  const limited = assessments.filter(a => a.transportable === 'limited').length;
  const unknown = assessments.filter(a => a.transportable === 'unknown').length;

  return {
    studies: assessments,
    summary: {
      total: studies.length,
      transportable,
      conditional,
      limited,
      unknown
    },
    overallAssessment: transportable / studies.length > 0.7
      ? 'Good transportability'
      : conditional / studies.length > 0.5
        ? 'Moderate transportability - consider subgroup analysis'
        : 'Limited transportability - results may not generalize'
  };
}

/**
 * Calibration Weighting for Target Population
 * Uses calibration constraints to match target moments
 *
 * @param {Array} trialData - IPD from trial
 * @param {Array} targetMoments - Target population moments [{name, mean}]
 * @param {Object} options - Options
 * @returns {Object} Calibration weights
 */
export function calibrationWeighting(trialData, targetMoments, options = {}) {
  const { maxIterations = 100, tolerance = 1e-6 } = options;

  if (!trialData.length || !targetMoments.length) {
    return { error: 'Trial data and target moments required' };
  }

  const n = trialData.length;
  const m = targetMoments.length;

  // Extract covariate values
  const X = targetMoments.map(tm => {
    return trialData.map(p => p.covariates[tm.name] || 0);
  });

  const targets = targetMoments.map(tm => tm.mean);

  // Initialize weights uniformly
  let weights = new Array(n).fill(1);

  // Iterative proportional fitting (raking)
  for (let iter = 0; iter < maxIterations; iter++) {
    let maxChange = 0;

    for (let j = 0; j < m; j++) {
      // Current weighted mean for covariate j
      const sumWX = X[j].reduce((sum, x, i) => sum + weights[i] * x, 0);
      const sumW = weights.reduce((a, b) => a + b, 0);
      const currentMean = sumWX / sumW;

      // Adjustment factor
      const ratio = targets[j] / (currentMean || 1);

      // Update weights
      for (let i = 0; i < n; i++) {
        const oldWeight = weights[i];
        // Adjust based on covariate value relative to mean
        const adjustment = 1 + (ratio - 1) * (X[j][i] / (currentMean || 1) - 1) * 0.1;
        weights[i] *= Math.max(0.1, Math.min(10, adjustment));
        maxChange = Math.max(maxChange, Math.abs(weights[i] - oldWeight));
      }

      // Normalize weights
      const newSumW = weights.reduce((a, b) => a + b, 0);
      weights = weights.map(w => w * n / newSumW);
    }

    if (maxChange < tolerance) break;
  }

  // Check calibration achieved
  const calibrationCheck = targetMoments.map(tm => {
    const idx = targetMoments.indexOf(tm);
    const sumWX = X[idx].reduce((sum, x, i) => sum + weights[i] * x, 0);
    const sumW = weights.reduce((a, b) => a + b, 0);
    const achieved = sumWX / sumW;

    return {
      name: tm.name,
      target: tm.mean,
      achieved,
      difference: achieved - tm.mean,
      calibrated: Math.abs(achieved - tm.mean) < 0.01
    };
  });

  const ESS = (weights.reduce((sum, w) => sum + w, 0) ** 2) /
    weights.reduce((sum, w) => sum + w * w, 0);

  return {
    method: 'Calibration Weighting',
    n,
    weights: {
      values: weights,
      mean: weights.reduce((a, b) => a + b, 0) / n,
      min: Math.min(...weights),
      max: Math.max(...weights),
      cv: calculateCV(weights)
    },
    ESS,
    calibrationCheck,
    allCalibrated: calibrationCheck.every(c => c.calibrated)
  };
}

/**
 * Generalizability Index
 * Calculates how well trial results generalize to target population
 *
 * @param {Object} trialResult - Trial treatment effect {es, se}
 * @param {Object} transportedResult - Transported treatment effect {es, se}
 * @param {number} ESS - Effective sample size after weighting
 * @param {number} originalN - Original sample size
 * @returns {Object} Generalizability assessment
 */
export function calculateGeneralizabilityIndex(trialResult, transportedResult, ESS, originalN) {
  // Effect change
  const effectChange = transportedResult.es - trialResult.es;
  const relativeChange = trialResult.es !== 0
    ? effectChange / Math.abs(trialResult.es)
    : effectChange;

  // Precision loss
  const precisionRatio = transportedResult.se / trialResult.se;

  // ESS ratio
  const essRatio = ESS / originalN;

  // Composite generalizability index (0-100)
  // Higher is better (more generalizable, less change needed)
  const stabilityScore = Math.max(0, 100 - Math.abs(relativeChange) * 100);
  const precisionScore = Math.max(0, 100 * (2 - precisionRatio)); // Penalize SE inflation
  const essScore = essRatio * 100;

  const generalizabilityIndex = (stabilityScore * 0.4 + precisionScore * 0.3 + essScore * 0.3);

  return {
    effectChange,
    relativeChange: relativeChange * 100,  // As percentage
    precisionRatio,
    essRatio,
    generalizabilityIndex: Math.max(0, Math.min(100, generalizabilityIndex)),
    interpretation: generalizabilityIndex > 80
      ? 'High generalizability - trial results likely apply to target population'
      : generalizabilityIndex > 60
        ? 'Moderate generalizability - some caution warranted'
        : generalizabilityIndex > 40
          ? 'Limited generalizability - results may differ in target population'
          : 'Poor generalizability - significant differences expected'
  };
}

function getTransportabilityRecommendation(avgAbsSmd, metrics) {
  if (avgAbsSmd < 0.1) {
    return 'Populations are similar. Direct application of results is reasonable.';
  } else if (avgAbsSmd < 0.25) {
    const imbalanced = metrics.filter(m => m.imbalanced).map(m => m.name);
    return `Moderate differences in ${imbalanced.join(', ')}. Consider adjustment or subgroup analysis.`;
  } else {
    return 'Substantial population differences. Transportability analysis strongly recommended before applying results.';
  }
}

/**
 * Target Population Effect Estimation
 * Estimates treatment effect for target population using multiple methods
 *
 * @param {Array} studies - Studies with effects and population characteristics
 * @param {Array} targetChars - Target population characteristics
 * @returns {Object} Target population effect estimates
 */
export function estimateTargetPopulationEffect(studies, targetChars) {
  const transportability = assessStudyTransportability(studies, targetChars);

  // Method 1: Simple subset (transportable studies only)
  const transportableStudies = studies.filter((s, i) =>
    transportability.studies[i]?.transportable === 'yes'
  );

  let subsetEstimate = null;
  if (transportableStudies.length >= 2) {
    const w = transportableStudies.map(s => 1 / s.vi);
    const sumW = w.reduce((a, b) => a + b, 0);
    const pooledES = w.reduce((sum, wi, i) => sum + wi * transportableStudies[i].es, 0) / sumW;
    const pooledSE = Math.sqrt(1 / sumW);

    subsetEstimate = {
      method: 'Transportable Subset',
      k: transportableStudies.length,
      es: pooledES,
      se: pooledSE,
      ci: {
        lo: pooledES - 1.96 * pooledSE,
        hi: pooledES + 1.96 * pooledSE
      }
    };
  }

  // Method 2: Weighted by population similarity
  const similarityWeights = transportability.studies.map(a => {
    if (a.transportable === 'unknown') return 0.5;
    return Math.max(0.1, 1 - (a.avgAbsSmd || 0));
  });

  const validStudies = studies.filter((s, i) =>
    !s.excluded && s.es != null && s.vi > 0
  );

  if (validStudies.length >= 2) {
    const combinedWeights = validStudies.map((s, i) => {
      const baseWeight = 1 / s.vi;
      const simWeight = similarityWeights[i] || 0.5;
      return baseWeight * simWeight;
    });

    const sumCW = combinedWeights.reduce((a, b) => a + b, 0);
    const weightedES = combinedWeights.reduce((sum, w, i) =>
      sum + w * validStudies[i].es, 0) / sumCW;

    // Approximate SE (conservative)
    const baseW = validStudies.map(s => 1 / s.vi);
    const sumBaseW = baseW.reduce((a, b) => a + b, 0);
    const baseSE = Math.sqrt(1 / sumBaseW);
    const avgSimWeight = similarityWeights.reduce((a, b) => a + b, 0) / similarityWeights.length;
    const adjustedSE = baseSE / avgSimWeight;  // Inflate for uncertainty

    var similarityWeightedEstimate = {
      method: 'Similarity-Weighted',
      k: validStudies.length,
      es: weightedES,
      se: adjustedSE,
      ci: {
        lo: weightedES - 1.96 * adjustedSE,
        hi: weightedES + 1.96 * adjustedSE
      }
    };
  }

  return {
    transportabilityAssessment: transportability,
    estimates: {
      subset: subsetEstimate,
      similarityWeighted: similarityWeightedEstimate || null
    },
    recommendation: subsetEstimate && similarityWeightedEstimate
      ? Math.abs(subsetEstimate.es - similarityWeightedEstimate.es) < 0.1
        ? 'Estimates are consistent - reasonable confidence in transportability'
        : 'Estimates differ - interpret target population effect with caution'
      : 'Limited data for target population estimation'
  };
}

export default {
  calculatePopulationSimilarity,
  inverseOddsWeighting,
  assessStudyTransportability,
  calibrationWeighting,
  calculateGeneralizabilityIndex,
  estimateTargetPopulationEffect
};
