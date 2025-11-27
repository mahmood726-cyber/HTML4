/**
 * Diagnostic Test Accuracy (DTA) Meta-Analysis
 *
 * Implements methods for meta-analysis of diagnostic accuracy studies:
 * - Bivariate model (Reitsma et al.)
 * - HSROC model (Rutter & Gatsonis)
 * - Summary ROC curves
 * - Heterogeneity assessment
 * - QUADAS-2 quality integration
 *
 * Based on Cochrane DTA Handbook and NICE diagnostics guidelines
 */

import { jStat } from 'jstat';

/**
 * Calculate diagnostic accuracy measures from 2x2 table
 *
 * @param {Object} data - 2x2 table {tp, fp, fn, tn}
 * @returns {Object} Accuracy measures with CIs
 */
export function calculateAccuracyMeasures(data) {
  const { tp, fp, fn, tn } = data;

  // Validate
  if (tp < 0 || fp < 0 || fn < 0 || tn < 0) {
    return { error: 'All cells must be non-negative' };
  }

  const totalDiseased = tp + fn;
  const totalHealthy = fp + tn;
  const totalPositive = tp + fp;
  const totalNegative = fn + tn;
  const total = tp + fp + fn + tn;

  if (totalDiseased === 0 || totalHealthy === 0) {
    return { error: 'Both diseased and healthy groups must have observations' };
  }

  // Sensitivity and Specificity
  const sens = tp / totalDiseased;
  const spec = tn / totalHealthy;

  // Standard errors using Wilson score interval components
  const seSens = Math.sqrt((sens * (1 - sens)) / totalDiseased);
  const seSpec = Math.sqrt((spec * (1 - spec)) / totalHealthy);

  // Positive and Negative Predictive Values
  const ppv = totalPositive > 0 ? tp / totalPositive : 0;
  const npv = totalNegative > 0 ? tn / totalNegative : 0;

  // Likelihood Ratios
  const lrPlus = spec < 1 ? sens / (1 - spec) : Infinity;
  const lrMinus = sens < 1 ? (1 - sens) / spec : 0;

  // Diagnostic Odds Ratio
  const dor = (fp === 0 || fn === 0)
    ? (tp * tn + 0.5) / ((fp + 0.5) * (fn + 0.5))
    : (tp * tn) / (fp * fn);

  // Log DOR and its SE
  const logDOR = Math.log(dor);
  const seLogDOR = Math.sqrt(1 / tp + 1 / fp + 1 / fn + 1 / tn);

  // Youden's Index
  const youden = sens + spec - 1;

  // Area under ROC (single study approximation)
  const auc = (sens + spec) / 2;

  // Wilson score CIs for proportions
  const sensCI = wilsonCI(tp, totalDiseased, 0.95);
  const specCI = wilsonCI(tn, totalHealthy, 0.95);

  return {
    tp, fp, fn, tn,
    n: total,
    prevalence: totalDiseased / total,
    sensitivity: {
      estimate: sens,
      se: seSens,
      ci: sensCI
    },
    specificity: {
      estimate: spec,
      se: seSpec,
      ci: specCI
    },
    ppv: {
      estimate: ppv,
      ci: wilsonCI(tp, totalPositive, 0.95)
    },
    npv: {
      estimate: npv,
      ci: wilsonCI(tn, totalNegative, 0.95)
    },
    lrPlus: {
      estimate: lrPlus,
      ci: lrCI(sens, spec, totalDiseased, totalHealthy, 'positive')
    },
    lrMinus: {
      estimate: lrMinus,
      ci: lrCI(sens, spec, totalDiseased, totalHealthy, 'negative')
    },
    dor: {
      estimate: dor,
      logDOR,
      seLogDOR,
      ci: {
        lo: Math.exp(logDOR - 1.96 * seLogDOR),
        hi: Math.exp(logDOR + 1.96 * seLogDOR)
      }
    },
    youden,
    auc
  };
}

/**
 * Wilson score confidence interval for proportions
 */
function wilsonCI(x, n, conf) {
  if (n === 0) return { lo: 0, hi: 1 };

  const p = x / n;
  const z = jStat.normal.inv((1 + conf) / 2, 0, 1);
  const z2 = z * z;

  const denom = 1 + z2 / n;
  const center = (p + z2 / (2 * n)) / denom;
  const margin = (z / denom) * Math.sqrt(p * (1 - p) / n + z2 / (4 * n * n));

  return {
    lo: Math.max(0, center - margin),
    hi: Math.min(1, center + margin)
  };
}

/**
 * Likelihood ratio confidence interval
 */
function lrCI(sens, spec, nDiseased, nHealthy, type) {
  let lr, seLnLR;

  if (type === 'positive') {
    lr = sens / (1 - spec);
    if (spec === 1 || sens === 0) return { lo: 0, hi: Infinity };
    seLnLR = Math.sqrt(
      (1 - sens) / (sens * nDiseased) +
      spec / ((1 - spec) * nHealthy)
    );
  } else {
    lr = (1 - sens) / spec;
    if (sens === 1 || spec === 0) return { lo: 0, hi: Infinity };
    seLnLR = Math.sqrt(
      sens / ((1 - sens) * nDiseased) +
      (1 - spec) / (spec * nHealthy)
    );
  }

  const lnLR = Math.log(lr);
  return {
    lo: Math.exp(lnLR - 1.96 * seLnLR),
    hi: Math.exp(lnLR + 1.96 * seLnLR)
  };
}

/**
 * Transform sensitivity/specificity to logit scale
 */
function logitTransform(studies) {
  return studies.map(s => {
    const acc = calculateAccuracyMeasures(s);
    if (acc.error) return null;

    const sens = acc.sensitivity.estimate;
    const spec = acc.specificity.estimate;

    // Add continuity correction for 0/1
    const adjSens = Math.max(0.001, Math.min(0.999, sens));
    const adjSpec = Math.max(0.001, Math.min(0.999, spec));

    const logitSens = Math.log(adjSens / (1 - adjSens));
    const logitSpec = Math.log(adjSpec / (1 - adjSpec));

    // Variance of logit
    const nDis = s.tp + s.fn;
    const nHealthy = s.fp + s.tn;
    const varLogitSens = 1 / (sens * (1 - sens) * nDis);
    const varLogitSpec = 1 / (spec * (1 - spec) * nHealthy);

    return {
      ...s,
      sens, spec,
      logitSens, logitSpec,
      varLogitSens, varLogitSpec,
      seLogitSens: Math.sqrt(varLogitSens),
      seLogitSpec: Math.sqrt(varLogitSpec)
    };
  }).filter(s => s !== null);
}

/**
 * Bivariate Random Effects Model for DTA Meta-Analysis
 * (Reitsma et al. 2005)
 *
 * @param {Array} studies - Studies with {id, tp, fp, fn, tn}
 * @param {Object} options - Options
 * @returns {Object} Bivariate model results
 */
export function bivariateModel(studies, options = {}) {
  const { maxIter = 100, tolerance = 1e-6 } = options;

  const validStudies = studies.filter(s =>
    !s.excluded && s.tp != null && s.fp != null && s.fn != null && s.tn != null
  );

  if (validStudies.length < 4) {
    return { error: 'At least 4 studies required for bivariate model' };
  }

  // Transform to logit scale
  const transformed = logitTransform(validStudies);
  const k = transformed.length;

  // Method of moments estimation for bivariate model
  // Estimate mean logit sensitivity and specificity

  // Initial estimates (simple pooling)
  const totalTP = transformed.reduce((sum, s) => sum + s.tp, 0);
  const totalFN = transformed.reduce((sum, s) => sum + s.fn, 0);
  const totalFP = transformed.reduce((sum, s) => sum + s.fp, 0);
  const totalTN = transformed.reduce((sum, s) => sum + s.tn, 0);

  const pooledSens = totalTP / (totalTP + totalFN);
  const pooledSpec = totalTN / (totalTN + totalFP);

  // Weighted means of logit values
  const wSens = transformed.map(s => 1 / s.varLogitSens);
  const wSpec = transformed.map(s => 1 / s.varLogitSpec);

  const sumWSens = wSens.reduce((a, b) => a + b, 0);
  const sumWSpec = wSpec.reduce((a, b) => a + b, 0);

  let muSens = wSens.reduce((sum, w, i) => sum + w * transformed[i].logitSens, 0) / sumWSens;
  let muSpec = wSpec.reduce((sum, w, i) => sum + w * transformed[i].logitSpec, 0) / sumWSpec;

  // Estimate between-study variances (DerSimonian-Laird type)
  let QSens = 0, QSpec = 0;
  for (let i = 0; i < k; i++) {
    QSens += wSens[i] * (transformed[i].logitSens - muSens) ** 2;
    QSpec += wSpec[i] * (transformed[i].logitSpec - muSpec) ** 2;
  }

  const cSens = sumWSens - wSens.reduce((sum, w) => sum + w * w, 0) / sumWSens;
  const cSpec = sumWSpec - wSpec.reduce((sum, w) => sum + w * w, 0) / sumWSpec;

  let tau2Sens = Math.max(0, (QSens - (k - 1)) / cSens);
  let tau2Spec = Math.max(0, (QSpec - (k - 1)) / cSpec);

  // Correlation between logit sens and logit spec (typically negative)
  let sumCov = 0;
  for (let i = 0; i < k; i++) {
    sumCov += (transformed[i].logitSens - muSens) * (transformed[i].logitSpec - muSpec);
  }
  const covariance = sumCov / (k - 1);
  const correlation = covariance / (Math.sqrt(tau2Sens + 0.001) * Math.sqrt(tau2Spec + 0.001));

  // Iterative re-estimation with random effects
  for (let iter = 0; iter < maxIter; iter++) {
    const wSensRE = transformed.map(s => 1 / (s.varLogitSens + tau2Sens));
    const wSpecRE = transformed.map(s => 1 / (s.varLogitSpec + tau2Spec));

    const sumWSensRE = wSensRE.reduce((a, b) => a + b, 0);
    const sumWSpecRE = wSpecRE.reduce((a, b) => a + b, 0);

    const newMuSens = wSensRE.reduce((sum, w, i) => sum + w * transformed[i].logitSens, 0) / sumWSensRE;
    const newMuSpec = wSpecRE.reduce((sum, w, i) => sum + w * transformed[i].logitSpec, 0) / sumWSpecRE;

    if (Math.abs(newMuSens - muSens) < tolerance && Math.abs(newMuSpec - muSpec) < tolerance) {
      break;
    }

    muSens = newMuSens;
    muSpec = newMuSpec;
  }

  // Final weights and variances
  const finalWSens = transformed.map(s => 1 / (s.varLogitSens + tau2Sens));
  const finalWSpec = transformed.map(s => 1 / (s.varLogitSpec + tau2Spec));
  const sumFinalWSens = finalWSens.reduce((a, b) => a + b, 0);
  const sumFinalWSpec = finalWSpec.reduce((a, b) => a + b, 0);

  const varMuSens = 1 / sumFinalWSens;
  const varMuSpec = 1 / sumFinalWSpec;

  // Transform back to probability scale
  const summarySens = 1 / (1 + Math.exp(-muSens));
  const summarySpec = 1 / (1 + Math.exp(-muSpec));

  // Delta method for SE on probability scale
  const seSens = Math.sqrt(varMuSens) * summarySens * (1 - summarySens);
  const seSpec = Math.sqrt(varMuSpec) * summarySpec * (1 - summarySpec);

  // I² statistics
  const I2Sens = tau2Sens / (tau2Sens + 1 / sumWSens) * 100;
  const I2Spec = tau2Spec / (tau2Spec + 1 / sumWSpec) * 100;

  return {
    method: 'Bivariate Random Effects',
    k,
    sensitivity: {
      logit: muSens,
      estimate: summarySens,
      se: seSens,
      ci: {
        lo: Math.max(0, summarySens - 1.96 * seSens),
        hi: Math.min(1, summarySens + 1.96 * seSens)
      }
    },
    specificity: {
      logit: muSpec,
      estimate: summarySpec,
      se: seSpec,
      ci: {
        lo: Math.max(0, summarySpec - 1.96 * seSpec),
        hi: Math.min(1, summarySpec + 1.96 * seSpec)
      }
    },
    heterogeneity: {
      tau2Sens,
      tau2Spec,
      correlation,
      I2Sens,
      I2Spec
    },
    summaryLRPlus: summarySens / (1 - summarySpec),
    summaryLRMinus: (1 - summarySens) / summarySpec,
    summaryDOR: (summarySens * summarySpec) / ((1 - summarySens) * (1 - summarySpec)),
    studies: transformed.map(s => ({
      id: s.id,
      sens: s.sens,
      spec: s.spec,
      logitSens: s.logitSens,
      logitSpec: s.logitSpec
    }))
  };
}

/**
 * Generate Summary ROC (SROC) Curve Data
 *
 * @param {Object} bivarResult - Result from bivariateModel
 * @param {Object} options - Options including number of points
 * @returns {Object} SROC curve data for plotting
 */
export function generateSROCCurve(bivarResult, options = {}) {
  const { nPoints = 100 } = options;

  if (bivarResult.error) {
    return { error: bivarResult.error };
  }

  const muSens = bivarResult.sensitivity.logit;
  const muSpec = bivarResult.specificity.logit;
  const tau2Sens = bivarResult.heterogeneity.tau2Sens;
  const tau2Spec = bivarResult.heterogeneity.tau2Spec;
  const rho = bivarResult.heterogeneity.correlation;

  // Generate SROC curve points
  // Using the bivariate normal on logit scale
  const curve = [];
  const confRegion = [];

  // FPR from 0 to 1
  for (let i = 0; i <= nPoints; i++) {
    const fpr = i / nPoints;
    const logitFPR = Math.log((fpr + 0.001) / (1 - fpr + 0.001));

    // Conditional mean of logit sensitivity given logit(1-specificity)
    // Assumes bivariate normal relationship
    const expectedLogitSens = muSens + rho * Math.sqrt(tau2Sens / tau2Spec) *
      (logitFPR - (-muSpec));

    const sens = 1 / (1 + Math.exp(-expectedLogitSens));

    curve.push({
      fpr: fpr,
      sens: Math.max(0, Math.min(1, sens)),
      spec: 1 - fpr
    });

    // Confidence region (approximate)
    const seSens = Math.sqrt(tau2Sens) * sens * (1 - sens);
    confRegion.push({
      fpr,
      sensLo: Math.max(0, sens - 1.96 * seSens),
      sensHi: Math.min(1, sens + 1.96 * seSens)
    });
  }

  // Summary operating point
  const summaryPoint = {
    fpr: 1 - bivarResult.specificity.estimate,
    sens: bivarResult.sensitivity.estimate,
    spec: bivarResult.specificity.estimate
  };

  // AUC calculation (trapezoidal)
  let auc = 0;
  for (let i = 1; i < curve.length; i++) {
    const width = curve[i].fpr - curve[i - 1].fpr;
    const avgHeight = (curve[i].sens + curve[i - 1].sens) / 2;
    auc += width * avgHeight;
  }

  return {
    curve,
    confRegion,
    summaryPoint,
    summaryPointCI: {
      sensCI: bivarResult.sensitivity.ci,
      specCI: bivarResult.specificity.ci
    },
    auc,
    studies: bivarResult.studies.map(s => ({
      id: s.id,
      fpr: 1 - s.spec,
      sens: s.sens
    }))
  };
}

/**
 * HSROC Model (Hierarchical Summary ROC)
 * Rutter & Gatsonis parameterization
 *
 * @param {Array} studies - Studies with 2x2 data
 * @param {Object} options - Options
 * @returns {Object} HSROC model results
 */
export function hsrocModel(studies, options = {}) {
  const validStudies = studies.filter(s =>
    !s.excluded && s.tp != null && s.fp != null && s.fn != null && s.tn != null
  );

  if (validStudies.length < 4) {
    return { error: 'At least 4 studies required for HSROC model' };
  }

  const transformed = logitTransform(validStudies);
  const k = transformed.length;

  // HSROC parameterization:
  // D = logit(sens) + logit(spec) (accuracy)
  // S = logit(sens) - logit(spec) (threshold)

  const D = transformed.map(s => s.logitSens + s.logitSpec);
  const S = transformed.map(s => s.logitSens - s.logitSpec);

  // Estimate HSROC parameters
  const meanD = D.reduce((a, b) => a + b, 0) / k;
  const meanS = S.reduce((a, b) => a + b, 0) / k;

  // Variances
  const varD = D.reduce((sum, d) => sum + (d - meanD) ** 2, 0) / (k - 1);
  const varS = S.reduce((sum, s) => sum + (s - meanS) ** 2, 0) / (k - 1);

  // Shape parameter (beta in HSROC)
  // Simplified: regress D on S
  let sumDS = 0, sumS2 = 0;
  for (let i = 0; i < k; i++) {
    sumDS += (D[i] - meanD) * (S[i] - meanS);
    sumS2 += (S[i] - meanS) ** 2;
  }

  const beta = sumS2 > 0 ? sumDS / sumS2 : 0;

  // Lambda (accuracy parameter)
  // exp(Lambda/2) is the DOR when threshold = mean
  const Lambda = meanD - beta * meanS;

  // Theta* (mean threshold)
  const ThetaStar = meanS;

  // Convert to summary sensitivity and specificity at mean threshold
  const summaryLogitSpec = (Lambda - ThetaStar) / 2;
  const summaryLogitSens = (Lambda + ThetaStar) / 2;

  const summarySens = 1 / (1 + Math.exp(-summaryLogitSens));
  const summarySpec = 1 / (1 + Math.exp(-summaryLogitSpec));

  // Summary DOR
  const summaryLogDOR = Lambda;
  const summaryDOR = Math.exp(summaryLogDOR);

  return {
    method: 'HSROC',
    k,
    parameters: {
      Lambda,          // Accuracy
      ThetaStar,       // Mean threshold
      beta,            // Shape (asymmetry)
      varD,
      varS
    },
    sensitivity: {
      estimate: summarySens,
      logit: summaryLogitSens
    },
    specificity: {
      estimate: summarySpec,
      logit: summaryLogitSpec
    },
    summaryDOR: {
      estimate: summaryDOR,
      logDOR: summaryLogDOR
    },
    interpretation: {
      symmetry: Math.abs(beta) < 0.5 ? 'Symmetric SROC curve' : 'Asymmetric SROC curve',
      accuracy: Lambda > 3 ? 'High accuracy' : Lambda > 2 ? 'Moderate accuracy' : 'Low accuracy'
    },
    studies: transformed.map((s, i) => ({
      id: s.id,
      D: D[i],
      S: S[i],
      sens: s.sens,
      spec: s.spec
    }))
  };
}

/**
 * Meta-Regression for DTA (Covariate analysis)
 *
 * @param {Array} studies - Studies with 2x2 data and covariates
 * @param {string} covariateName - Name of covariate to analyze
 * @returns {Object} Meta-regression results
 */
export function dtaMetaRegression(studies, covariateName) {
  const validStudies = studies.filter(s =>
    !s.excluded && s.tp != null && s[covariateName] != null
  );

  if (validStudies.length < 5) {
    return { error: 'At least 5 studies required for DTA meta-regression' };
  }

  const transformed = logitTransform(validStudies);
  const k = transformed.length;

  const x = validStudies.map(s => s[covariateName]);
  const yS = transformed.map(s => s.logitSens);
  const yP = transformed.map(s => s.logitSpec);

  // Simple weighted regression for sensitivity
  const wS = transformed.map(s => 1 / s.varLogitSens);
  const sumWS = wS.reduce((a, b) => a + b, 0);
  const meanX = wS.reduce((sum, w, i) => sum + w * x[i], 0) / sumWS;
  const meanYS = wS.reduce((sum, w, i) => sum + w * yS[i], 0) / sumWS;

  let SxxS = 0, SxyS = 0;
  for (let i = 0; i < k; i++) {
    SxxS += wS[i] * (x[i] - meanX) ** 2;
    SxyS += wS[i] * (x[i] - meanX) * (yS[i] - meanYS);
  }

  const slopeSens = SxxS > 0 ? SxyS / SxxS : 0;
  const interceptSens = meanYS - slopeSens * meanX;

  // Same for specificity
  const wP = transformed.map(s => 1 / s.varLogitSpec);
  const sumWP = wP.reduce((a, b) => a + b, 0);
  const meanYP = wP.reduce((sum, w, i) => sum + w * yP[i], 0) / sumWP;

  let SxxP = 0, SxyP = 0;
  for (let i = 0; i < k; i++) {
    SxxP += wP[i] * (x[i] - meanX) ** 2;
    SxyP += wP[i] * (x[i] - meanX) * (yP[i] - meanYP);
  }

  const slopeSpec = SxxP > 0 ? SxyP / SxxP : 0;
  const interceptSpec = meanYP - slopeSpec * meanX;

  // R² calculations
  let SStotS = 0, SSresS = 0;
  let SStotP = 0, SSresP = 0;
  for (let i = 0; i < k; i++) {
    SStotS += wS[i] * (yS[i] - meanYS) ** 2;
    SSresS += wS[i] * (yS[i] - (interceptSens + slopeSens * x[i])) ** 2;
    SStotP += wP[i] * (yP[i] - meanYP) ** 2;
    SSresP += wP[i] * (yP[i] - (interceptSpec + slopeSpec * x[i])) ** 2;
  }

  const R2Sens = SStotS > 0 ? 1 - SSresS / SStotS : 0;
  const R2Spec = SStotP > 0 ? 1 - SSresP / SStotP : 0;

  return {
    method: 'DTA Meta-Regression',
    covariate: covariateName,
    k,
    sensitivity: {
      intercept: interceptSens,
      slope: slopeSens,
      R2: R2Sens * 100,
      significant: Math.abs(slopeSens) > 0.1
    },
    specificity: {
      intercept: interceptSpec,
      slope: slopeSpec,
      R2: R2Spec * 100,
      significant: Math.abs(slopeSpec) > 0.1
    },
    interpretation: generateDTARegressionInterpretation(slopeSens, slopeSpec, covariateName)
  };
}

function generateDTARegressionInterpretation(slopeSens, slopeSpec, covariateName) {
  const parts = [];

  if (Math.abs(slopeSens) > 0.1) {
    parts.push(`${covariateName} ${slopeSens > 0 ? 'increases' : 'decreases'} sensitivity`);
  }
  if (Math.abs(slopeSpec) > 0.1) {
    parts.push(`${covariateName} ${slopeSpec > 0 ? 'increases' : 'decreases'} specificity`);
  }

  if (parts.length === 0) {
    return `No significant association between ${covariateName} and test accuracy`;
  }

  return parts.join('; ');
}

/**
 * QUADAS-2 Weighted Analysis
 * Adjusts for study quality in DTA meta-analysis
 *
 * @param {Array} studies - Studies with quality scores
 * @param {Object} options - Options
 * @returns {Object} Quality-weighted results
 */
export function quadasWeightedAnalysis(studies, options = {}) {
  const validStudies = studies.filter(s =>
    !s.excluded && s.tp != null && s.quadas != null
  );

  if (validStudies.length < 3) {
    return { error: 'At least 3 studies with QUADAS scores required' };
  }

  // QUADAS-2 domains: patient selection, index test, reference standard, flow and timing
  // Score each domain as: low risk = 1, unclear = 0.5, high risk = 0

  const qualityWeights = validStudies.map(s => {
    const q = s.quadas;
    const domains = ['patientSelection', 'indexTest', 'referenceStandard', 'flowTiming'];
    let score = 0;

    for (const domain of domains) {
      if (q[domain] === 'low') score += 1;
      else if (q[domain] === 'unclear') score += 0.5;
      // high risk adds 0
    }

    return score / domains.length;  // 0 to 1
  });

  // Run bivariate model with quality weights
  const bivar = bivariateModel(validStudies);

  if (bivar.error) {
    return { error: bivar.error };
  }

  // Quality-weighted pooling
  const transformed = logitTransform(validStudies);
  const k = transformed.length;

  // Combine inverse-variance weights with quality weights
  const combinedWSens = transformed.map((s, i) =>
    (1 / s.varLogitSens) * qualityWeights[i]
  );
  const combinedWSpec = transformed.map((s, i) =>
    (1 / s.varLogitSpec) * qualityWeights[i]
  );

  const sumCWSens = combinedWSens.reduce((a, b) => a + b, 0);
  const sumCWSpec = combinedWSpec.reduce((a, b) => a + b, 0);

  const qwLogitSens = combinedWSens.reduce((sum, w, i) =>
    sum + w * transformed[i].logitSens, 0) / sumCWSens;
  const qwLogitSpec = combinedWSpec.reduce((sum, w, i) =>
    sum + w * transformed[i].logitSpec, 0) / sumCWSpec;

  const qwSens = 1 / (1 + Math.exp(-qwLogitSens));
  const qwSpec = 1 / (1 + Math.exp(-qwLogitSpec));

  return {
    method: 'QUADAS-2 Weighted',
    k,
    qualityScores: validStudies.map((s, i) => ({
      id: s.id,
      score: qualityWeights[i],
      domains: s.quadas
    })),
    avgQualityScore: qualityWeights.reduce((a, b) => a + b, 0) / k,
    unweighted: {
      sensitivity: bivar.sensitivity.estimate,
      specificity: bivar.specificity.estimate
    },
    qualityWeighted: {
      sensitivity: qwSens,
      specificity: qwSpec
    },
    sensitivityChange: qwSens - bivar.sensitivity.estimate,
    specificityChange: qwSpec - bivar.specificity.estimate,
    interpretation: Math.abs(qwSens - bivar.sensitivity.estimate) > 0.05 ||
      Math.abs(qwSpec - bivar.specificity.estimate) > 0.05
        ? 'Quality weighting substantially changes estimates - interpret with caution'
        : 'Quality weighting has limited impact - results appear robust'
  };
}

export default {
  calculateAccuracyMeasures,
  bivariateModel,
  generateSROCCurve,
  hsrocModel,
  dtaMetaRegression,
  quadasWeightedAnalysis
};
