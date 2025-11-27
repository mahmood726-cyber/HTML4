/**
 * Robust Bayesian Meta-Analysis (RoBMA)
 *
 * Implements Bayesian model averaging across:
 * - Presence/absence of effect
 * - Presence/absence of heterogeneity
 * - Presence/absence of publication bias
 *
 * Based on Bartoš & Maier (2020), Maier et al. (2022)
 */

import { jStat } from 'jstat';

/**
 * Robust Bayesian Meta-Analysis
 * Model-averages across different assumptions
 *
 * @param {Array} studies - Studies with es and se
 * @param {Object} options - Options including priors
 * @returns {Object} RoBMA results
 */
export function robustBayesianMA(studies, options = {}) {
  const {
    effectPrior = { type: 'normal', mean: 0, sd: 1 },
    heterogeneityPrior = { type: 'invGamma', shape: 1, scale: 0.15 },
    biasPrior = { type: 'twoSided' },  // Publication bias model
    nSamples = 5000,
    seed = null
  } = options;

  const validStudies = studies.filter(s => !s.excluded && s.es != null && s.se > 0);
  const k = validStudies.length;

  if (k < 3) {
    return { error: 'At least 3 studies required for RoBMA' };
  }

  const y = validStudies.map(s => s.es);
  const se = validStudies.map(s => s.se);

  // Define model space
  // Each model is defined by (hasEffect, hasHeterogeneity, hasBias)
  const models = [
    { hasEffect: false, hasHeterogeneity: false, hasBias: false, name: 'Null' },
    { hasEffect: true, hasHeterogeneity: false, hasBias: false, name: 'FE' },
    { hasEffect: false, hasHeterogeneity: true, hasBias: false, name: 'RE-Null' },
    { hasEffect: true, hasHeterogeneity: true, hasBias: false, name: 'RE' },
    { hasEffect: false, hasHeterogeneity: false, hasBias: true, name: 'Null+Bias' },
    { hasEffect: true, hasHeterogeneity: false, hasBias: true, name: 'FE+Bias' },
    { hasEffect: false, hasHeterogeneity: true, hasBias: true, name: 'RE-Null+Bias' },
    { hasEffect: true, hasHeterogeneity: true, hasBias: true, name: 'RE+Bias' }
  ];

  // Calculate marginal likelihood for each model (via Laplace approximation)
  const modelResults = models.map(model => {
    const result = fitSingleModel(y, se, model, effectPrior, heterogeneityPrior, nSamples);
    return { ...model, ...result };
  });

  // Prior model probabilities (equal by default)
  const priorProbs = models.map(() => 1 / models.length);

  // Calculate posterior model probabilities
  const logMarginalLikelihoods = modelResults.map(m => m.logMarginalLikelihood);
  const maxLogML = Math.max(...logMarginalLikelihoods);

  const unnormalizedPosterior = modelResults.map((m, i) =>
    priorProbs[i] * Math.exp(m.logMarginalLikelihood - maxLogML)
  );
  const sumPosterior = unnormalizedPosterior.reduce((a, b) => a + b, 0);
  const posteriorProbs = unnormalizedPosterior.map(p => p / sumPosterior);

  // Model-averaged estimates
  let maEffect = 0;
  let maEffectVar = 0;
  let maTau2 = 0;

  for (let i = 0; i < models.length; i++) {
    maEffect += posteriorProbs[i] * modelResults[i].effectEstimate;
    maTau2 += posteriorProbs[i] * modelResults[i].tau2Estimate;
  }

  // Variance includes both within-model and between-model variance
  for (let i = 0; i < models.length; i++) {
    const withinVar = modelResults[i].effectSE ** 2;
    const betweenVar = (modelResults[i].effectEstimate - maEffect) ** 2;
    maEffectVar += posteriorProbs[i] * (withinVar + betweenVar);
  }

  const maEffectSE = Math.sqrt(maEffectVar);

  // Inclusion Bayes Factors
  // P(effect | data) / P(no effect | data) divided by prior odds
  const posteriorEffect = posteriorProbs.filter((_, i) => models[i].hasEffect)
    .reduce((a, b) => a + b, 0);
  const posteriorNoEffect = 1 - posteriorEffect;
  const priorEffect = 0.5;  // Assuming equal prior

  const bfEffect = (posteriorEffect / posteriorNoEffect) / (priorEffect / (1 - priorEffect));

  const posteriorHet = posteriorProbs.filter((_, i) => models[i].hasHeterogeneity)
    .reduce((a, b) => a + b, 0);
  const bfHeterogeneity = posteriorHet > 0.001 && posteriorHet < 0.999
    ? (posteriorHet / (1 - posteriorHet)) / 1 : posteriorHet > 0.5 ? Infinity : 0;

  const posteriorBias = posteriorProbs.filter((_, i) => models[i].hasBias)
    .reduce((a, b) => a + b, 0);
  const bfBias = posteriorBias > 0.001 && posteriorBias < 0.999
    ? (posteriorBias / (1 - posteriorBias)) / 1 : posteriorBias > 0.5 ? Infinity : 0;

  // Credible interval for model-averaged effect
  // Approximate using normal distribution
  const ci95Lo = maEffect - 1.96 * maEffectSE;
  const ci95Hi = maEffect + 1.96 * maEffectSE;

  return {
    method: 'Robust Bayesian Meta-Analysis',
    k,
    modelAveraged: {
      effect: maEffect,
      se: maEffectSE,
      ci95: { lo: ci95Lo, hi: ci95Hi },
      tau2: maTau2,
      tau: Math.sqrt(maTau2)
    },
    inclusionBayesFactors: {
      effect: {
        BF10: bfEffect,
        interpretation: interpretBF(bfEffect)
      },
      heterogeneity: {
        BF10: bfHeterogeneity,
        interpretation: interpretBF(bfHeterogeneity)
      },
      publicationBias: {
        BF10: bfBias,
        interpretation: interpretBF(bfBias)
      }
    },
    posteriorProbabilities: {
      effect: posteriorEffect,
      heterogeneity: posteriorHet,
      publicationBias: posteriorBias
    },
    models: modelResults.map((m, i) => ({
      name: m.name,
      priorProb: priorProbs[i],
      posteriorProb: posteriorProbs[i],
      logML: m.logMarginalLikelihood,
      effect: m.effectEstimate,
      tau2: m.tau2Estimate
    })),
    interpretation: generateRoBMAInterpretation(bfEffect, posteriorEffect, posteriorBias)
  };
}

/**
 * Fit a single model component
 */
function fitSingleModel(y, se, model, effectPrior, heterogeneityPrior, nSamples) {
  const k = y.length;
  const v = se.map(s => s * s);

  let effectEstimate = 0;
  let effectSE = 0.1;
  let tau2Estimate = 0;

  if (!model.hasEffect) {
    // Null model: effect = 0
    effectEstimate = 0;
    effectSE = 0.001;
  } else if (!model.hasHeterogeneity) {
    // Fixed effect model
    const w = v.map(vi => 1 / vi);
    const sumW = w.reduce((a, b) => a + b, 0);
    effectEstimate = w.reduce((sum, wi, i) => sum + wi * y[i], 0) / sumW;
    effectSE = Math.sqrt(1 / sumW);
  } else {
    // Random effects model
    // Estimate tau² using DL
    const w = v.map(vi => 1 / vi);
    const sumW = w.reduce((a, b) => a + b, 0);
    const muFE = w.reduce((sum, wi, i) => sum + wi * y[i], 0) / sumW;
    const Q = w.reduce((sum, wi, i) => sum + wi * (y[i] - muFE) ** 2, 0);
    const C = sumW - w.reduce((sum, wi) => sum + wi * wi, 0) / sumW;
    tau2Estimate = Math.max(0, (Q - (k - 1)) / C);

    // RE estimate
    const wRE = v.map(vi => 1 / (vi + tau2Estimate));
    const sumWRE = wRE.reduce((a, b) => a + b, 0);
    effectEstimate = wRE.reduce((sum, wi, i) => sum + wi * y[i], 0) / sumWRE;
    effectSE = Math.sqrt(1 / sumWRE);
  }

  // Apply selection model if hasBias
  if (model.hasBias) {
    // Simple selection adjustment: weight by inverse of selection probability
    // Assuming studies with p < 0.05 are selected with probability 1
    // and others with probability 0.5

    const zScores = y.map((yi, i) => Math.abs(yi / se[i]));
    const selectionWeights = zScores.map(z =>
      z > 1.96 ? 1 : 0.5  // Significant vs not
    );

    const adjustedW = v.map((vi, i) => selectionWeights[i] / (vi + tau2Estimate));
    const sumAdjW = adjustedW.reduce((a, b) => a + b, 0);

    if (model.hasEffect && sumAdjW > 0) {
      effectEstimate = adjustedW.reduce((sum, wi, i) => sum + wi * y[i], 0) / sumAdjW;
      effectSE = Math.sqrt(1 / sumAdjW) * 1.2;  // Inflate SE for bias adjustment
    }
  }

  // Calculate log marginal likelihood (Laplace approximation)
  let logML = 0;

  // Likelihood contribution
  for (let i = 0; i < k; i++) {
    const totalVar = v[i] + tau2Estimate;
    const resid = y[i] - effectEstimate;
    logML += -0.5 * Math.log(2 * Math.PI * totalVar) - 0.5 * resid * resid / totalVar;
  }

  // Prior contribution
  if (model.hasEffect) {
    // Normal prior on effect
    const priorVar = effectPrior.sd ** 2;
    logML += -0.5 * Math.log(2 * Math.PI * priorVar) -
      0.5 * (effectEstimate - effectPrior.mean) ** 2 / priorVar;
  }

  if (model.hasHeterogeneity && tau2Estimate > 0) {
    // Inverse gamma prior on tau²
    const shape = heterogeneityPrior.shape;
    const scale = heterogeneityPrior.scale;
    logML += shape * Math.log(scale) - jStat.gammafn(shape) -
      (shape + 1) * Math.log(tau2Estimate) - scale / tau2Estimate;
  }

  // Penalty for publication bias model complexity
  if (model.hasBias) {
    logML -= 1;  // BIC-like penalty
  }

  return {
    effectEstimate,
    effectSE,
    tau2Estimate,
    logMarginalLikelihood: logML
  };
}

/**
 * Interpret Bayes Factor
 */
function interpretBF(bf) {
  if (!isFinite(bf)) return 'Extreme evidence';
  if (bf > 100) return 'Extreme evidence for H1';
  if (bf > 30) return 'Very strong evidence for H1';
  if (bf > 10) return 'Strong evidence for H1';
  if (bf > 3) return 'Moderate evidence for H1';
  if (bf > 1) return 'Anecdotal evidence for H1';
  if (bf > 1 / 3) return 'Anecdotal evidence for H0';
  if (bf > 1 / 10) return 'Moderate evidence for H0';
  if (bf > 1 / 30) return 'Strong evidence for H0';
  if (bf > 1 / 100) return 'Very strong evidence for H0';
  return 'Extreme evidence for H0';
}

function generateRoBMAInterpretation(bfEffect, posteriorEffect, posteriorBias) {
  const parts = [];

  if (bfEffect > 3) {
    parts.push(`Evidence for effect presence (BF=${bfEffect.toFixed(1)}).`);
  } else if (bfEffect < 1 / 3) {
    parts.push(`Evidence against effect (BF=${bfEffect.toFixed(2)}).`);
  } else {
    parts.push('Inconclusive evidence for effect.');
  }

  parts.push(`P(effect|data) = ${(posteriorEffect * 100).toFixed(1)}%.`);

  if (posteriorBias > 0.5) {
    parts.push('Publication bias likely present.');
  }

  return parts.join(' ');
}

/**
 * Simplified Bayesian Model Averaging for effect/no-effect
 *
 * @param {Array} studies - Studies
 * @param {Object} options - Options
 * @returns {Object} BMA results
 */
export function bayesianModelAveraging(studies, options = {}) {
  const {
    effectPriorSD = 0.5,
    priorProbEffect = 0.5
  } = options;

  const validStudies = studies.filter(s => !s.excluded && s.es != null && s.se > 0);
  const k = validStudies.length;

  if (k < 2) {
    return { error: 'At least 2 studies required' };
  }

  const y = validStudies.map(s => s.es);
  const v = validStudies.map(s => s.se * s.se);

  // Model 0: No effect (mu = 0)
  let logML0 = 0;
  for (let i = 0; i < k; i++) {
    logML0 += -0.5 * Math.log(2 * Math.PI * v[i]) - 0.5 * y[i] * y[i] / v[i];
  }

  // Model 1: Effect present
  // Marginal likelihood with normal prior on mu
  const w = v.map(vi => 1 / vi);
  const sumW = w.reduce((a, b) => a + b, 0);
  const muHat = w.reduce((sum, wi, i) => sum + wi * y[i], 0) / sumW;

  const priorVar = effectPriorSD ** 2;
  const posteriorVar = 1 / (sumW + 1 / priorVar);
  const posteriorMean = posteriorVar * (sumW * muHat);

  let logML1 = -0.5 * Math.log(priorVar) + 0.5 * Math.log(posteriorVar);
  logML1 += 0.5 * posteriorMean * posteriorMean / posteriorVar;
  for (let i = 0; i < k; i++) {
    logML1 += -0.5 * Math.log(2 * Math.PI * v[i]) - 0.5 * y[i] * y[i] / v[i];
  }

  // Bayes factor
  const logBF10 = logML1 - logML0;
  const BF10 = Math.exp(Math.min(logBF10, 100));  // Cap to avoid overflow

  // Posterior probability
  const priorOdds = priorProbEffect / (1 - priorProbEffect);
  const posteriorOdds = priorOdds * BF10;
  const posteriorProbEffect = posteriorOdds / (1 + posteriorOdds);

  // Model-averaged estimate
  const maEffect = posteriorProbEffect * posteriorMean + (1 - posteriorProbEffect) * 0;

  return {
    method: 'Bayesian Model Averaging',
    k,
    BF10,
    logBF10,
    interpretation: interpretBF(BF10),
    posteriorProbEffect,
    posteriorProbNull: 1 - posteriorProbEffect,
    modelAveragedEffect: maEffect,
    conditionalEffect: {
      mean: posteriorMean,
      sd: Math.sqrt(posteriorVar),
      ci95: {
        lo: posteriorMean - 1.96 * Math.sqrt(posteriorVar),
        hi: posteriorMean + 1.96 * Math.sqrt(posteriorVar)
      }
    }
  };
}

export default {
  robustBayesianMA,
  bayesianModelAveraging
};
