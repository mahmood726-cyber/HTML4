/**
 * PET-PEESE Bias Correction Methods
 *
 * Implements small-study effect correction using:
 * - PET (Precision-Effect Test): Regresses effect on SE
 * - PEESE (Precision-Effect Estimate with Standard Error): Regresses effect on variance
 * - Selection models for publication bias
 * - Trim-and-fill enhancement
 *
 * Based on Stanley & Doucouliagos methodology and meta-epidemiology-bias-correction approaches
 */

import { jStat } from 'jstat';

/**
 * PET (Precision-Effect Test)
 * Regresses effect sizes on their standard errors
 * Tests H0: no true effect (intercept = 0 when SE = 0)
 *
 * @param {Array} studies - Studies with es (effect size) and se (standard error)
 * @param {Object} options - Options
 * @returns {Object} PET results
 */
export function runPET(studies, options = {}) {
  const validStudies = studies.filter(s => !s.excluded && s.es != null && s.se > 0);
  const k = validStudies.length;

  if (k < 3) {
    return { error: 'At least 3 studies required for PET' };
  }

  // Extract data
  const y = validStudies.map(s => s.es);
  const se = validStudies.map(s => s.se);

  // Weights: inverse variance
  const w = se.map(s => 1 / (s * s));

  // Weighted least squares regression: ES = β0 + β1*SE
  // Using SE as predictor (PET model)
  const sumW = w.reduce((a, b) => a + b, 0);
  const sumWX = w.reduce((sum, wi, i) => sum + wi * se[i], 0);
  const sumWY = w.reduce((sum, wi, i) => sum + wi * y[i], 0);
  const sumWXX = w.reduce((sum, wi, i) => sum + wi * se[i] * se[i], 0);
  const sumWXY = w.reduce((sum, wi, i) => sum + wi * se[i] * y[i], 0);

  const meanX = sumWX / sumW;
  const meanY = sumWY / sumW;

  // Calculate slope and intercept
  const Sxx = sumWXX - sumWX * sumWX / sumW;
  const Sxy = sumWXY - sumWX * sumWY / sumW;

  if (Math.abs(Sxx) < 1e-10) {
    return { error: 'Insufficient variance in standard errors' };
  }

  const slope = Sxy / Sxx;
  const intercept = meanY - slope * meanX;

  // Residuals and MSE
  let SSres = 0;
  for (let i = 0; i < k; i++) {
    const fitted = intercept + slope * se[i];
    const resid = y[i] - fitted;
    SSres += w[i] * resid * resid;
  }

  const df = k - 2;
  const MSE = SSres / df;

  // Standard errors of coefficients
  const seIntercept = Math.sqrt(MSE * sumWXX / (sumW * Sxx));
  const seSlope = Math.sqrt(MSE * sumW / (sumW * Sxx));

  // t-statistics and p-values
  const tIntercept = intercept / seIntercept;
  const tSlope = slope / seSlope;

  const pIntercept = 2 * (1 - jStat.studentt.cdf(Math.abs(tIntercept), df));
  const pSlope = 2 * (1 - jStat.studentt.cdf(Math.abs(tSlope), df));

  // Confidence intervals
  const tCrit = jStat.studentt.inv(0.975, df);
  const interceptCI = {
    lo: intercept - tCrit * seIntercept,
    hi: intercept + tCrit * seIntercept
  };
  const slopeCI = {
    lo: slope - tCrit * seSlope,
    hi: slope + tCrit * seSlope
  };

  // PET-corrected estimate is the intercept (effect when SE = 0)
  const correctedES = intercept;
  const correctedSE = seIntercept;

  // Test for small-study effects (is slope significant?)
  const hasSmallStudyEffect = pSlope < 0.10;

  return {
    method: 'PET',
    k,
    intercept: {
      estimate: intercept,
      se: seIntercept,
      t: tIntercept,
      p: pIntercept,
      ci: interceptCI
    },
    slope: {
      estimate: slope,
      se: seSlope,
      t: tSlope,
      p: pSlope,
      ci: slopeCI
    },
    correctedES,
    correctedSE,
    correctedCI: interceptCI,
    hasSmallStudyEffect,
    interpretation: hasSmallStudyEffect
      ? 'Significant small-study effect detected. PET-corrected estimate may be more accurate.'
      : 'No significant small-study effect. Standard meta-analysis may be appropriate.',
    studies: validStudies.map((s, i) => ({
      id: s.id,
      es: y[i],
      se: se[i],
      fitted: intercept + slope * se[i],
      residual: y[i] - (intercept + slope * se[i])
    }))
  };
}

/**
 * PEESE (Precision-Effect Estimate with Standard Error)
 * Regresses effect sizes on variance (SE²)
 * Better for estimating true effect when there IS publication bias
 *
 * @param {Array} studies - Studies with es and se
 * @param {Object} options - Options
 * @returns {Object} PEESE results
 */
export function runPEESE(studies, options = {}) {
  const validStudies = studies.filter(s => !s.excluded && s.es != null && s.se > 0);
  const k = validStudies.length;

  if (k < 3) {
    return { error: 'At least 3 studies required for PEESE' };
  }

  const y = validStudies.map(s => s.es);
  const se = validStudies.map(s => s.se);
  const variance = se.map(s => s * s);

  // Weights
  const w = se.map(s => 1 / (s * s));

  // WLS regression: ES = β0 + β1*Variance
  const sumW = w.reduce((a, b) => a + b, 0);
  const sumWX = w.reduce((sum, wi, i) => sum + wi * variance[i], 0);
  const sumWY = w.reduce((sum, wi, i) => sum + wi * y[i], 0);
  const sumWXX = w.reduce((sum, wi, i) => sum + wi * variance[i] * variance[i], 0);
  const sumWXY = w.reduce((sum, wi, i) => sum + wi * variance[i] * y[i], 0);

  const meanX = sumWX / sumW;
  const meanY = sumWY / sumW;

  const Sxx = sumWXX - sumWX * sumWX / sumW;
  const Sxy = sumWXY - sumWX * sumWY / sumW;

  if (Math.abs(Sxx) < 1e-10) {
    return { error: 'Insufficient variance in study variances' };
  }

  const slope = Sxy / Sxx;
  const intercept = meanY - slope * meanX;

  // Residuals
  let SSres = 0;
  for (let i = 0; i < k; i++) {
    const fitted = intercept + slope * variance[i];
    const resid = y[i] - fitted;
    SSres += w[i] * resid * resid;
  }

  const df = k - 2;
  const MSE = SSres / df;

  const seIntercept = Math.sqrt(MSE * sumWXX / (sumW * Sxx));
  const seSlope = Math.sqrt(MSE * sumW / (sumW * Sxx));

  const tIntercept = intercept / seIntercept;
  const tSlope = slope / seSlope;

  const pIntercept = 2 * (1 - jStat.studentt.cdf(Math.abs(tIntercept), df));
  const pSlope = 2 * (1 - jStat.studentt.cdf(Math.abs(tSlope), df));

  const tCrit = jStat.studentt.inv(0.975, df);
  const interceptCI = {
    lo: intercept - tCrit * seIntercept,
    hi: intercept + tCrit * seIntercept
  };

  return {
    method: 'PEESE',
    k,
    intercept: {
      estimate: intercept,
      se: seIntercept,
      t: tIntercept,
      p: pIntercept,
      ci: interceptCI
    },
    slope: {
      estimate: slope,
      se: seSlope,
      t: tSlope,
      p: pSlope
    },
    correctedES: intercept,
    correctedSE: seIntercept,
    correctedCI: interceptCI,
    studies: validStudies.map((s, i) => ({
      id: s.id,
      es: y[i],
      variance: variance[i],
      fitted: intercept + slope * variance[i],
      residual: y[i] - (intercept + slope * variance[i])
    }))
  };
}

/**
 * PET-PEESE Conditional Estimator
 * Uses PET to test for bias, then:
 * - If PET intercept is significant: use PEESE estimate
 * - If PET intercept is not significant: conclude no effect
 *
 * @param {Array} studies - Studies with es and se
 * @param {Object} options - Options including alpha threshold
 * @returns {Object} Combined PET-PEESE results
 */
export function petPeese(studies, options = {}) {
  const { alpha = 0.10 } = options;

  const pet = runPET(studies);
  if (pet.error) return pet;

  const peese = runPEESE(studies);
  if (peese.error) return peese;

  // Decision rule:
  // If PET intercept is NOT significant → conclude no effect
  // If PET intercept IS significant → use PEESE estimate
  const petSignificant = pet.intercept.p < alpha;

  let finalEstimate, finalSE, finalCI, decision;

  if (!petSignificant) {
    // PET not significant: conclude no effect after bias correction
    finalEstimate = 0;
    finalSE = pet.correctedSE;
    finalCI = { lo: 0, hi: 0 };
    decision = 'PET-insignificant';
  } else {
    // PET significant: use PEESE
    finalEstimate = peese.correctedES;
    finalSE = peese.correctedSE;
    finalCI = peese.correctedCI;
    decision = 'PEESE';
  }

  return {
    method: 'PET-PEESE',
    k: pet.k,
    pet,
    peese,
    decision,
    correctedES: finalEstimate,
    correctedSE: finalSE,
    correctedCI: finalCI,
    interpretation: decision === 'PET-insignificant'
      ? 'After accounting for publication bias, no significant effect remains.'
      : `Using PEESE estimate: ${finalEstimate.toFixed(3)} [${finalCI.lo.toFixed(3)}, ${finalCI.hi.toFixed(3)}]`
  };
}

/**
 * Selection Model for Publication Bias
 * Implements a simple step-function selection model
 * based on p-value cutoffs
 *
 * @param {Array} studies - Studies with es and se
 * @param {Object} options - Options including p-value cutoffs
 * @returns {Object} Selection model results
 */
export function selectionModel(studies, options = {}) {
  const {
    cutoffs = [0.025, 0.5, 1.0],  // Two-tailed: 0.05, NS
    weights = null  // If null, weights are estimated
  } = options;

  const validStudies = studies.filter(s => !s.excluded && s.es != null && s.se > 0);
  const k = validStudies.length;

  if (k < 5) {
    return { error: 'At least 5 studies required for selection model' };
  }

  const y = validStudies.map(s => s.es);
  const se = validStudies.map(s => s.se);

  // Calculate z-scores and one-sided p-values
  const z = y.map((yi, i) => yi / se[i]);
  const pvals = z.map(zi => 1 - jStat.normal.cdf(Math.abs(zi), 0, 1));

  // Categorize studies by significance
  const categories = pvals.map(p => {
    for (let i = 0; i < cutoffs.length; i++) {
      if (p <= cutoffs[i]) return i;
    }
    return cutoffs.length - 1;
  });

  // Count studies in each category
  const counts = cutoffs.map((_, i) => categories.filter(c => c === i).length);

  // Simple estimate: weighted by inverse of assumed selection probability
  // More sophisticated models would use maximum likelihood

  // For now, use reweighting based on observed proportions
  const totalStudies = k;
  const expectedProp = cutoffs.map((c, i) => {
    const prev = i === 0 ? 0 : cutoffs[i - 1];
    return c - prev;
  });

  const observedProp = counts.map(c => c / totalStudies);

  // Selection weights (inverse of selection probability ratio)
  const selectionWeights = expectedProp.map((exp, i) => {
    const obs = observedProp[i];
    if (obs < 0.01) return 1;  // Avoid division by zero
    return exp / obs;
  });

  // Apply selection weights to studies
  const adjustedWeights = validStudies.map((s, i) => {
    const cat = categories[i];
    const baseW = 1 / (se[i] * se[i]);
    return baseW * selectionWeights[cat];
  });

  // Weighted mean with adjusted weights
  const sumAdjW = adjustedWeights.reduce((a, b) => a + b, 0);
  const adjustedES = adjustedWeights.reduce((sum, w, i) => sum + w * y[i], 0) / sumAdjW;
  const adjustedSE = Math.sqrt(1 / sumAdjW);

  // Calculate unadjusted for comparison
  const baseW = se.map(s => 1 / (s * s));
  const sumBaseW = baseW.reduce((a, b) => a + b, 0);
  const unadjustedES = baseW.reduce((sum, w, i) => sum + w * y[i], 0) / sumBaseW;

  return {
    method: 'Selection Model',
    k,
    cutoffs,
    categoryCounts: counts,
    expectedProportions: expectedProp,
    observedProportions: observedProp,
    selectionWeights,
    unadjustedES,
    adjustedES,
    adjustedSE,
    adjustedCI: {
      lo: adjustedES - 1.96 * adjustedSE,
      hi: adjustedES + 1.96 * adjustedSE
    },
    biasEstimate: unadjustedES - adjustedES,
    interpretation: Math.abs(unadjustedES - adjustedES) > 0.1
      ? 'Substantial publication bias detected. Consider adjusted estimate.'
      : 'Limited evidence of publication bias from selection model.'
  };
}

/**
 * Enhanced Trim-and-Fill with Effect Direction
 * Extends basic trim-and-fill with proper direction estimation
 *
 * @param {Array} studies - Studies with es and se
 * @param {Object} options - Options
 * @returns {Object} Enhanced trim-and-fill results
 */
export function enhancedTrimFill(studies, options = {}) {
  const {
    side = 'auto',  // 'left', 'right', or 'auto'
    estimator = 'R0',  // 'L0', 'R0', or 'Q0'
    maxIter = 100
  } = options;

  const validStudies = studies.filter(s => !s.excluded && s.es != null && s.se > 0);
  const k = validStudies.length;

  if (k < 3) {
    return { error: 'At least 3 studies required for trim-and-fill' };
  }

  const y = validStudies.map(s => s.es);
  const se = validStudies.map(s => s.se);
  const w = se.map(s => 1 / (s * s));

  // Initial pooled estimate
  const sumW = w.reduce((a, b) => a + b, 0);
  let pooledES = w.reduce((sum, wi, i) => sum + wi * y[i], 0) / sumW;

  // Determine side
  let fillSide = side;
  if (side === 'auto') {
    // Check asymmetry direction using Egger's test slope
    const meanSE = se.reduce((a, b) => a + b, 0) / k;
    const correlation = calculateCorrelation(y, se);
    fillSide = correlation > 0 ? 'left' : 'right';
  }

  // Iterative trim-and-fill
  let k0 = 0;
  let imputed = [];

  for (let iter = 0; iter < maxIter; iter++) {
    // Rank studies by deviation from pooled estimate
    const deviations = y.map((yi, i) => ({
      idx: i,
      y: yi,
      se: se[i],
      dev: yi - pooledES,
      rank: 0
    }));

    // Sort by deviation
    if (fillSide === 'right') {
      deviations.sort((a, b) => b.dev - a.dev);  // Most positive first
    } else {
      deviations.sort((a, b) => a.dev - b.dev);  // Most negative first
    }

    // Assign ranks
    deviations.forEach((d, i) => d.rank = i + 1);

    // Estimate k0 using R0 estimator
    const n = k;
    let sumRanks = 0;
    deviations.forEach(d => {
      if ((fillSide === 'right' && d.dev > 0) || (fillSide === 'left' && d.dev < 0)) {
        sumRanks += d.rank;
      }
    });

    // R0 estimator
    const newK0 = Math.round(Math.max(0, 4 * sumRanks / n - n));

    if (newK0 === k0) break;
    k0 = newK0;

    // Create imputed studies (mirror around pooled estimate)
    imputed = [];
    const extremeStudies = deviations.slice(0, k0);

    for (const study of extremeStudies) {
      const mirroredES = 2 * pooledES - study.y;
      imputed.push({
        id: `Imputed_${study.idx}`,
        es: mirroredES,
        se: study.se,
        imputed: true
      });
    }

    // Recalculate pooled estimate with imputed studies
    const allY = [...y, ...imputed.map(s => s.es)];
    const allSE = [...se, ...imputed.map(s => s.se)];
    const allW = allSE.map(s => 1 / (s * s));
    const sumAllW = allW.reduce((a, b) => a + b, 0);
    pooledES = allW.reduce((sum, wi, i) => sum + wi * allY[i], 0) / sumAllW;
  }

  // Final calculation
  const allY = [...y, ...imputed.map(s => s.es)];
  const allSE = [...se, ...imputed.map(s => s.se)];
  const allW = allSE.map(s => 1 / (s * s));
  const sumAllW = allW.reduce((a, b) => a + b, 0);
  const adjustedES = allW.reduce((sum, wi, i) => sum + wi * allY[i], 0) / sumAllW;
  const adjustedSE = Math.sqrt(1 / sumAllW);

  // Unadjusted
  const unadjustedES = w.reduce((sum, wi, i) => sum + wi * y[i], 0) / sumW;

  return {
    method: 'Trim-and-Fill',
    k,
    k0,
    side: fillSide,
    unadjustedES,
    adjustedES,
    adjustedSE,
    adjustedCI: {
      lo: adjustedES - 1.96 * adjustedSE,
      hi: adjustedES + 1.96 * adjustedSE
    },
    imputedStudies: imputed,
    biasEstimate: unadjustedES - adjustedES,
    interpretation: k0 > 0
      ? `${k0} studies imputed to restore symmetry. Adjusted estimate: ${adjustedES.toFixed(3)}`
      : 'No missing studies detected. Funnel plot appears symmetric.'
  };
}

/**
 * Calculate correlation between two arrays
 */
function calculateCorrelation(x, y) {
  const n = x.length;
  const meanX = x.reduce((a, b) => a + b, 0) / n;
  const meanY = y.reduce((a, b) => a + b, 0) / n;

  let sumXY = 0, sumX2 = 0, sumY2 = 0;
  for (let i = 0; i < n; i++) {
    const dx = x[i] - meanX;
    const dy = y[i] - meanY;
    sumXY += dx * dy;
    sumX2 += dx * dx;
    sumY2 += dy * dy;
  }

  const denom = Math.sqrt(sumX2 * sumY2);
  return denom < 1e-10 ? 0 : sumXY / denom;
}

/**
 * Comprehensive Publication Bias Assessment
 * Runs multiple methods and provides overall assessment
 *
 * @param {Array} studies - Studies with es and se
 * @returns {Object} Comprehensive bias assessment
 */
export function comprehensiveBiasAssessment(studies) {
  const pet = runPET(studies);
  const peese = runPEESE(studies);
  const petPeeseResult = petPeese(studies);
  const trimFill = enhancedTrimFill(studies);

  // Count indicators of bias
  let biasIndicators = 0;
  const indicators = [];

  if (!pet.error && pet.hasSmallStudyEffect) {
    biasIndicators++;
    indicators.push('PET slope significant');
  }

  if (!trimFill.error && trimFill.k0 > 0) {
    biasIndicators++;
    indicators.push(`Trim-and-fill imputed ${trimFill.k0} studies`);
  }

  // Calculate bias severity
  const estimates = [];
  if (!pet.error) estimates.push({ method: 'PET', es: pet.correctedES });
  if (!peese.error) estimates.push({ method: 'PEESE', es: peese.correctedES });
  if (!trimFill.error) estimates.push({ method: 'Trim-Fill', es: trimFill.adjustedES });

  const unadjusted = !pet.error ? pet.studies.reduce((sum, s) => {
    const w = 1 / (s.se * s.se);
    return { sumW: sum.sumW + w, sumWY: sum.sumWY + w * s.es };
  }, { sumW: 0, sumWY: 0 }) : null;

  const unadjustedES = unadjusted ? unadjusted.sumWY / unadjusted.sumW : null;

  let severity = 'none';
  if (biasIndicators >= 2 && estimates.length > 0) {
    const maxDiff = Math.max(...estimates.map(e => Math.abs(e.es - (unadjustedES || 0))));
    if (maxDiff > 0.3) severity = 'severe';
    else if (maxDiff > 0.1) severity = 'moderate';
    else severity = 'mild';
  } else if (biasIndicators === 1) {
    severity = 'possible';
  }

  return {
    k: pet.error ? 0 : pet.k,
    unadjustedES,
    pet: pet.error ? { error: pet.error } : {
      correctedES: pet.correctedES,
      slopeP: pet.slope.p,
      hasSmallStudyEffect: pet.hasSmallStudyEffect
    },
    peese: peese.error ? { error: peese.error } : {
      correctedES: peese.correctedES
    },
    petPeese: petPeeseResult.error ? { error: petPeeseResult.error } : {
      decision: petPeeseResult.decision,
      correctedES: petPeeseResult.correctedES
    },
    trimFill: trimFill.error ? { error: trimFill.error } : {
      k0: trimFill.k0,
      adjustedES: trimFill.adjustedES
    },
    biasIndicators,
    indicators,
    severity,
    recommendation: getRecommendation(severity, petPeeseResult)
  };
}

function getRecommendation(severity, petPeese) {
  switch (severity) {
    case 'severe':
      return 'Strong evidence of publication bias. Report bias-corrected estimates (PET-PEESE). Interpret with caution.';
    case 'moderate':
      return 'Moderate evidence of publication bias. Report both unadjusted and adjusted estimates. Consider sensitivity analyses.';
    case 'mild':
      return 'Mild evidence of publication bias. Results are likely robust but report bias assessment.';
    case 'possible':
      return 'Some indication of publication bias. Report bias tests and consider adjusted estimates.';
    default:
      return 'No strong evidence of publication bias. Standard estimates are likely appropriate.';
  }
}

export default {
  runPET,
  runPEESE,
  petPeese,
  selectionModel,
  enhancedTrimFill,
  comprehensiveBiasAssessment
};
