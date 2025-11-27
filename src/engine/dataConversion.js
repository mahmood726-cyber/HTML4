/**
 * Data Conversion Utilities for Meta-Analysis
 *
 * Methods for converting summary statistics:
 * - Median/IQR to Mean/SD (Wan, Luo, McGrath methods)
 * - Risk measures (OR, RR, HR conversions)
 * - Correlation coefficient transformations
 * - SMD variations (Cohen's d, Hedges' g, Glass's Δ)
 * - Sample size estimation from confidence intervals
 *
 * Based on Cochrane Handbook and methodological papers
 */

import { jStat } from 'jstat';

/**
 * Convert Median and IQR to Mean and SD
 * Using Wan et al. (2014) method
 *
 * @param {number} median - Median value
 * @param {number} q1 - First quartile (25th percentile)
 * @param {number} q3 - Third quartile (75th percentile)
 * @param {number} n - Sample size
 * @returns {Object} Estimated mean and SD
 */
export function medianIQRToMeanSD(median, q1, q3, n) {
  if (n < 1) {
    return { error: 'Sample size must be at least 1' };
  }

  // Wan et al. (2014) formula for mean estimation from median and IQR
  // Mean ≈ (q1 + median + q3) / 3 for normal distribution

  // More refined estimate
  const mean = (q1 + median + q3) / 3;

  // SD from IQR: SD ≈ IQR / 1.35 for normal distribution
  // Wan et al. provide more accurate formula based on sample size
  const iqr = q3 - q1;

  // Adjustment factor based on sample size
  let sdEstimate;
  if (n <= 15) {
    // Small sample adjustment
    sdEstimate = iqr / (2 * jStat.normal.inv(0.75, 0, 1));
  } else if (n <= 70) {
    // Medium sample
    sdEstimate = iqr / (2 * jStat.normal.inv((0.75 * n - 0.125) / (n + 0.25), 0, 1));
  } else {
    // Large sample: IQR / 1.35
    sdEstimate = iqr / 1.35;
  }

  return {
    method: 'Wan et al. (2014)',
    input: { median, q1, q3, n },
    mean: mean,
    sd: sdEstimate,
    se: sdEstimate / Math.sqrt(n),
    note: 'Assumes approximately normal distribution'
  };
}

/**
 * Convert Median and Range to Mean and SD
 * Using Hozo et al. (2005) and Wan et al. (2014) methods
 *
 * @param {number} median - Median value
 * @param {number} min - Minimum value
 * @param {number} max - Maximum value
 * @param {number} n - Sample size
 * @returns {Object} Estimated mean and SD
 */
export function medianRangeToMeanSD(median, min, max, n) {
  if (n < 1) {
    return { error: 'Sample size must be at least 1' };
  }

  const range = max - min;

  // Wan et al. (2014) improved estimators
  let mean, sd;

  if (n <= 15) {
    // Small sample: use Hozo formula
    mean = (min + 2 * median + max) / 4;
    // SD from range for small samples
    sd = range / (2 * jStat.normal.inv((n - 0.375) / (n + 0.25), 0, 1));
  } else if (n <= 70) {
    // Medium sample
    mean = (min + 2 * median + max) / 4 + (min - 2 * median + max) / (4 * n);
    sd = range / (2 * jStat.normal.inv((n - 0.375) / (n + 0.25), 0, 1));
  } else {
    // Large sample (n > 70)
    mean = (min + 2 * median + max) / 4;
    // For large n, range/4 approximates SD
    sd = range / 4;
  }

  return {
    method: 'Wan et al. (2014)',
    input: { median, min, max, n },
    mean: mean,
    sd: Math.max(sd, 0.001),  // Ensure positive
    se: sd / Math.sqrt(n),
    note: 'Assumes approximately normal distribution'
  };
}

/**
 * Convert Median, IQR, and Range to Mean and SD
 * Combining both pieces of information (Luo et al. 2018)
 *
 * @param {Object} data - {median, q1, q3, min, max, n}
 * @returns {Object} Estimated mean and SD
 */
export function combinedToMeanSD(data) {
  const { median, q1, q3, min, max, n } = data;

  if (n < 1) {
    return { error: 'Sample size must be at least 1' };
  }

  // Luo et al. (2018) optimal linear combination
  // Use all five numbers: min, q1, median, q3, max

  let mean, sd;

  if (n <= 50) {
    // Small to medium sample
    // Mean estimation using weighted combination
    const w1 = 2.2 / (2.2 + n * 0.05);  // Weight for range-based
    const w2 = 1 - w1;  // Weight for IQR-based

    const meanRange = (min + 2 * median + max) / 4;
    const meanIQR = (q1 + median + q3) / 3;

    mean = w1 * meanRange + w2 * meanIQR;

    // SD using both IQR and range
    const iqr = q3 - q1;
    const range = max - min;

    const sdIQR = iqr / 1.35;
    const sdRange = range / (2 * jStat.normal.inv((n - 0.375) / (n + 0.25), 0, 1));

    // Optimal weighting (Luo et al.)
    sd = 0.7 * sdIQR + 0.3 * sdRange;
  } else {
    // Large sample
    mean = (q1 + median + q3) / 3;
    sd = (q3 - q1) / 1.35;
  }

  return {
    method: 'Luo et al. (2018)',
    input: data,
    mean,
    sd: Math.max(sd, 0.001),
    se: sd / Math.sqrt(n),
    note: 'Optimal combination of IQR and range'
  };
}

/**
 * McGrath et al. (2020) quantile estimation method
 * More robust for skewed distributions
 *
 * @param {Object} data - Summary statistics
 * @param {string} distribution - Assumed distribution ('normal', 'lognormal', 'unknown')
 * @returns {Object} Mean and SD estimates
 */
export function mcgrathQuantileMethod(data, distribution = 'unknown') {
  const { median, q1, q3, min, max, n } = data;

  if (!median || !n) {
    return { error: 'Median and sample size required' };
  }

  // Check for skewness
  let skewIndicator = 0;
  if (q1 && q3) {
    const iqr = q3 - q1;
    const lowerHalf = median - q1;
    const upperHalf = q3 - median;
    skewIndicator = (upperHalf - lowerHalf) / iqr;
  }

  const isSkewed = Math.abs(skewIndicator) > 0.2;

  if (distribution === 'unknown' && isSkewed) {
    // Use log-normal approximation for skewed data
    return estimateFromLognormal(data);
  }

  // Normal distribution assumption
  if (q1 && q3 && min && max) {
    return combinedToMeanSD(data);
  } else if (q1 && q3) {
    return medianIQRToMeanSD(median, q1, q3, n);
  } else if (min && max) {
    return medianRangeToMeanSD(median, min, max, n);
  }

  return { error: 'Insufficient data for estimation' };
}

/**
 * Estimate mean and SD assuming log-normal distribution
 */
function estimateFromLognormal(data) {
  const { median, q1, q3, n } = data;

  if (!median || !q1 || !q3 || median <= 0) {
    return { error: 'Positive median and quartiles required for log-normal' };
  }

  // For log-normal: median = exp(μ)
  const mu = Math.log(median);

  // IQR on log scale corresponds to 2 * 0.6745 * σ
  const logIQR = Math.log(q3) - Math.log(q1);
  const sigma = logIQR / (2 * 0.6745);

  // Mean of log-normal = exp(μ + σ²/2)
  const mean = Math.exp(mu + sigma * sigma / 2);

  // SD of log-normal = mean * sqrt(exp(σ²) - 1)
  const sd = mean * Math.sqrt(Math.exp(sigma * sigma) - 1);

  return {
    method: 'Log-normal estimation',
    input: data,
    mean,
    sd,
    se: sd / Math.sqrt(n),
    logParams: { mu, sigma },
    note: 'Assumes log-normal distribution due to skewness'
  };
}

/**
 * Convert between effect size measures
 */

/**
 * Convert Odds Ratio to Risk Ratio
 *
 * @param {number} or - Odds ratio
 * @param {number} p0 - Baseline risk (control group event rate)
 * @returns {Object} Risk ratio and confidence interval
 */
export function oddsRatioToRiskRatio(or, p0, options = {}) {
  const { orCI = null, confLevel = 0.95 } = options;

  if (p0 < 0 || p0 > 1) {
    return { error: 'Baseline risk must be between 0 and 1' };
  }

  // RR = OR / (1 - p0 + p0 * OR)
  const rr = or / (1 - p0 + p0 * or);

  let result = {
    method: 'Zhang & Yu (1998)',
    or,
    baselineRisk: p0,
    rr
  };

  if (orCI) {
    // Transform CI using same formula
    const rrLo = orCI.lo / (1 - p0 + p0 * orCI.lo);
    const rrHi = orCI.hi / (1 - p0 + p0 * orCI.hi);
    result.rrCI = { lo: rrLo, hi: rrHi };
  }

  return result;
}

/**
 * Convert Risk Ratio to Odds Ratio
 *
 * @param {number} rr - Risk ratio
 * @param {number} p0 - Baseline risk
 * @returns {Object} Odds ratio
 */
export function riskRatioToOddsRatio(rr, p0) {
  if (p0 < 0 || p0 > 1) {
    return { error: 'Baseline risk must be between 0 and 1' };
  }

  // OR = RR * (1 - p0) / (1 - p0 * RR)
  const denominator = 1 - p0 * rr;
  if (denominator <= 0) {
    return { error: 'Invalid combination: RR and p0 imply p1 > 1' };
  }

  const or = rr * (1 - p0) / denominator;

  return {
    method: 'Algebraic conversion',
    rr,
    baselineRisk: p0,
    or
  };
}

/**
 * Convert Hazard Ratio to Odds Ratio (approximate)
 * Valid for short follow-up and rare events
 *
 * @param {number} hr - Hazard ratio
 * @param {number} eventRate - Overall event rate
 * @returns {Object} Approximate OR
 */
export function hazardRatioToOddsRatio(hr, eventRate) {
  // For rare events: OR ≈ HR
  // More generally: OR ≈ (1 - (1-eventRate)^HR) / (1 - (1-eventRate))

  if (eventRate < 0.1) {
    // Rare events approximation
    return {
      method: 'Rare events approximation',
      hr,
      or: hr,
      note: 'HR ≈ OR for rare events'
    };
  }

  // General approximation assuming exponential survival
  const p0 = eventRate;
  const p1 = 1 - Math.pow(1 - p0, hr);

  const or = (p1 / (1 - p1)) / (p0 / (1 - p0));

  return {
    method: 'Exponential survival approximation',
    hr,
    eventRate,
    or,
    note: 'Approximate conversion assuming constant hazards'
  };
}

/**
 * Convert between standardized mean difference variants
 */

/**
 * Convert Cohen's d to Hedges' g (bias-corrected SMD)
 *
 * @param {number} d - Cohen's d
 * @param {number} n1 - Group 1 sample size
 * @param {number} n2 - Group 2 sample size
 * @returns {Object} Hedges' g and correction factor
 */
export function cohensD_to_HedgesG(d, n1, n2) {
  const df = n1 + n2 - 2;

  // Correction factor J (approximation)
  const J = 1 - 3 / (4 * df - 1);

  const g = d * J;

  // Variance of g
  const vd = (n1 + n2) / (n1 * n2) + d * d / (2 * (n1 + n2));
  const vg = J * J * vd;

  return {
    d,
    g,
    J,
    seG: Math.sqrt(vg),
    n1,
    n2
  };
}

/**
 * Convert Point-Biserial Correlation to Cohen's d
 *
 * @param {number} r - Point-biserial correlation
 * @param {number} p - Proportion in one group (for unequal n)
 * @returns {Object} Cohen's d
 */
export function correlationToD(r, p = 0.5) {
  // d = 2r / sqrt(1 - r²) for equal groups
  // With correction for unequal groups
  const correction = Math.sqrt(p * (1 - p));

  const d = r / (correction * Math.sqrt(1 - r * r));

  return {
    r,
    d,
    groupProportion: p
  };
}

/**
 * Convert Cohen's d to Point-Biserial Correlation
 *
 * @param {number} d - Cohen's d
 * @param {number} n1 - Group 1 sample size
 * @param {number} n2 - Group 2 sample size
 * @returns {Object} Correlation
 */
export function dToCorrelation(d, n1, n2) {
  const a = (n1 + n2) * (n1 + n2) / (n1 * n2);

  // r = d / sqrt(d² + a)
  const r = d / Math.sqrt(d * d + a);

  return {
    d,
    r,
    n1,
    n2
  };
}

/**
 * Fisher's z transformation for correlations
 */
export function fisherZ(r) {
  const z = 0.5 * Math.log((1 + r) / (1 - r));
  return {
    r,
    z,
    se: null  // Needs n for SE
  };
}

export function fisherZWithSE(r, n) {
  const z = 0.5 * Math.log((1 + r) / (1 - r));
  const se = 1 / Math.sqrt(n - 3);

  return {
    r,
    z,
    se,
    ci: {
      lo: z - 1.96 * se,
      hi: z + 1.96 * se
    }
  };
}

export function inverseFisherZ(z) {
  const r = (Math.exp(2 * z) - 1) / (Math.exp(2 * z) + 1);
  return { z, r };
}

/**
 * Estimate sample size from confidence interval
 *
 * @param {number} effect - Point estimate
 * @param {number} ciLo - Lower CI bound
 * @param {number} ciHi - Upper CI bound
 * @param {string} metric - Type of effect ('smd', 'or', 'rr', 'rd')
 * @param {number} confLevel - Confidence level (default 0.95)
 * @returns {Object} Estimated sample size components
 */
export function estimateSampleSizeFromCI(effect, ciLo, ciHi, metric = 'smd', confLevel = 0.95) {
  const z = jStat.normal.inv((1 + confLevel) / 2, 0, 1);

  // CI width
  const ciWidth = ciHi - ciLo;

  // SE from CI
  const se = ciWidth / (2 * z);

  let nEstimate;

  switch (metric) {
    case 'smd':
    case 'd':
    case 'g':
      // For SMD: SE² ≈ 2(1 + d²/8) * 2/n assuming equal groups
      // n ≈ 4(1 + d²/8) / SE²
      const d = Math.abs(effect);
      nEstimate = 4 * (1 + d * d / 8) / (se * se);
      break;

    case 'or':
    case 'logor':
      // For log(OR): SE² ≈ 1/a + 1/b + 1/c + 1/d
      // Very rough: n ≈ 4 / SE² for balanced design with moderate event rates
      nEstimate = 4 / (se * se);
      break;

    case 'rr':
    case 'logrr':
      // Similar approximation
      nEstimate = 4 / (se * se);
      break;

    case 'rd':
      // For RD: SE² ≈ p(1-p)/n1 + p(1-p)/n2
      // Approximate with p ≈ 0.5
      nEstimate = 0.5 / (se * se);
      break;

    case 'r':
    case 'correlation':
      // For Fisher's z: SE = 1/sqrt(n-3)
      nEstimate = 1 / (se * se) + 3;
      break;

    default:
      nEstimate = 4 / (se * se);
  }

  return {
    effect,
    ciLo,
    ciHi,
    metric,
    se,
    estimatedTotalN: Math.round(nEstimate),
    estimatedPerGroup: Math.round(nEstimate / 2),
    note: 'Rough estimate assuming balanced groups'
  };
}

/**
 * Convert t-statistic to effect size
 *
 * @param {number} t - t-statistic
 * @param {number} n1 - Group 1 sample size
 * @param {number} n2 - Group 2 sample size
 * @returns {Object} Cohen's d and Hedges' g
 */
export function tToEffectSize(t, n1, n2) {
  // d = t * sqrt(1/n1 + 1/n2)
  const d = t * Math.sqrt(1 / n1 + 1 / n2);

  const hedges = cohensD_to_HedgesG(d, n1, n2);

  return {
    t,
    n1,
    n2,
    d,
    g: hedges.g,
    seD: Math.sqrt((n1 + n2) / (n1 * n2) + d * d / (2 * (n1 + n2))),
    seG: hedges.seG
  };
}

/**
 * Convert F-statistic to effect size (for 2-group comparison)
 *
 * @param {number} F - F-statistic
 * @param {number} n1 - Group 1 sample size
 * @param {number} n2 - Group 2 sample size
 * @returns {Object} Effect sizes
 */
export function fToEffectSize(F, n1, n2) {
  // For 2 groups: F = t², so t = sqrt(F)
  const t = Math.sqrt(F);
  return tToEffectSize(t, n1, n2);
}

/**
 * Convert p-value and sample size to approximate effect size
 * (last resort when only p-value is reported)
 *
 * @param {number} p - p-value (two-tailed)
 * @param {number} n - Total sample size
 * @param {string} direction - 'positive' or 'negative'
 * @returns {Object} Approximate effect size
 */
export function pValueToEffectSize(p, n, direction = 'positive') {
  // Convert p to z-score
  const z = jStat.normal.inv(1 - p / 2, 0, 1);

  // Approximate d = z / sqrt(n/4) for equal groups
  const d = (direction === 'negative' ? -1 : 1) * z / Math.sqrt(n / 4);

  return {
    p,
    n,
    z,
    approximateD: d,
    note: 'Very rough approximation - use with caution'
  };
}

/**
 * Batch convert studies with various reporting formats
 *
 * @param {Array} studies - Studies with heterogeneous formats
 * @returns {Array} Studies with standardized es and vi
 */
export function standardizeStudyData(studies) {
  return studies.map(study => {
    const result = { ...study };

    // Already has es and vi
    if (study.es != null && study.vi != null) {
      return result;
    }

    // Convert from mean/SD
    if (study.mean1 != null && study.mean2 != null && study.sd1 != null && study.sd2 != null) {
      const pooledSD = Math.sqrt(
        ((study.n1 - 1) * study.sd1 * study.sd1 + (study.n2 - 1) * study.sd2 * study.sd2) /
        (study.n1 + study.n2 - 2)
      );
      result.es = (study.mean1 - study.mean2) / pooledSD;
      result.vi = (study.n1 + study.n2) / (study.n1 * study.n2) +
        result.es * result.es / (2 * (study.n1 + study.n2));
      result.conversionMethod = 'mean/SD to SMD';
    }

    // Convert from median/IQR
    else if (study.median1 != null && study.median2 != null &&
      study.q1_1 != null && study.q3_1 != null &&
      study.q1_2 != null && study.q3_2 != null) {
      const conv1 = medianIQRToMeanSD(study.median1, study.q1_1, study.q3_1, study.n1);
      const conv2 = medianIQRToMeanSD(study.median2, study.q1_2, study.q3_2, study.n2);

      if (!conv1.error && !conv2.error) {
        const pooledSD = Math.sqrt(
          ((study.n1 - 1) * conv1.sd * conv1.sd + (study.n2 - 1) * conv2.sd * conv2.sd) /
          (study.n1 + study.n2 - 2)
        );
        result.es = (conv1.mean - conv2.mean) / pooledSD;
        result.vi = (study.n1 + study.n2) / (study.n1 * study.n2) +
          result.es * result.es / (2 * (study.n1 + study.n2));
        result.conversionMethod = 'median/IQR to SMD';
      }
    }

    // Convert from t-statistic
    else if (study.t != null && study.n1 != null && study.n2 != null) {
      const conv = tToEffectSize(study.t, study.n1, study.n2);
      result.es = conv.g;  // Use Hedges' g
      result.vi = conv.seG * conv.seG;
      result.conversionMethod = 't-statistic to SMD';
    }

    // Convert from OR
    else if (study.or != null && study.orSE != null) {
      result.es = Math.log(study.or);
      result.vi = study.orSE * study.orSE;
      result.conversionMethod = 'OR to log(OR)';
    }

    // Convert from correlation
    else if (study.r != null && study.n != null) {
      const conv = fisherZWithSE(study.r, study.n);
      result.es = conv.z;
      result.vi = conv.se * conv.se;
      result.conversionMethod = 'r to Fisher z';
    }

    return result;
  });
}

export default {
  medianIQRToMeanSD,
  medianRangeToMeanSD,
  combinedToMeanSD,
  mcgrathQuantileMethod,
  oddsRatioToRiskRatio,
  riskRatioToOddsRatio,
  hazardRatioToOddsRatio,
  cohensD_to_HedgesG,
  correlationToD,
  dToCorrelation,
  fisherZ,
  fisherZWithSE,
  inverseFisherZ,
  estimateSampleSizeFromCI,
  tToEffectSize,
  fToEffectSize,
  pValueToEffectSize,
  standardizeStudyData
};
