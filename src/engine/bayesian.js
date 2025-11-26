/**
 * Bayesian Meta-Analysis Methods
 * Provides Bayesian random-effects meta-analysis using MCMC approximation
 * References:
 * - Sutton & Abrams (2001) - Bayesian methods in meta-analysis
 * - Rhodes et al. (2015) - Predictive distributions for tau
 * - Turner et al. (2012) - Informative priors for heterogeneity
 * - Lilienthal et al. (2024) - Empirical priors for HTA
 */

import { chiSqP } from './distributions.js';

/**
 * Bayesian random-effects meta-analysis using Gibbs sampling approximation
 *
 * Model:
 *   yi | θi ~ N(θi, σi²)
 *   θi | μ, τ² ~ N(μ, τ²)
 *   μ ~ N(μ0, σ0²)      [prior for overall effect]
 *   τ² ~ Half-Cauchy(0, scale) or InverseGamma(a, b)
 *
 * @param {Array} effects - Array of { es, vi } effect estimates
 * @param {Object} options - Bayesian analysis options
 * @returns {Object} Posterior summaries
 */
export function bayesianMetaAnalysis(effects, options = {}) {
  const {
    // Prior for μ (overall effect)
    muPrior = { mean: 0, sd: 10 },  // Weakly informative

    // Prior for τ² (between-study variance)
    tauPrior = {
      type: 'half-cauchy',  // 'half-cauchy', 'inverse-gamma', 'uniform', 'informative'
      scale: 0.5,           // For half-Cauchy
      shape: 0.001,         // For inverse-gamma (a)
      rate: 0.001,          // For inverse-gamma (b)
      // For informative priors (Turner et al.)
      category: null        // 'pharmacological-placebo', 'pharmacological-active', etc.
    },

    // MCMC settings
    nIterations = 10000,
    nBurnin = 2000,
    nThin = 1,
    nChains = 2,

    // Output
    confLevel = 0.95,
    seed = null
  } = options;

  const active = effects.filter(e => !e.excluded && e.es !== null && e.vi > 0);

  if (active.length < 2) {
    return { error: 'Need at least 2 studies for meta-analysis' };
  }

  const k = active.length;
  const y = active.map(e => e.es);
  const sigma2 = active.map(e => e.vi);

  // Initialize random number generator
  const rng = createRNG(seed);

  // Get informative prior if specified
  let tauPriorParams = tauPrior;
  if (tauPrior.type === 'informative' && tauPrior.category) {
    tauPriorParams = getInformativePrior(tauPrior.category);
  }

  // Run multiple chains
  const chains = [];
  for (let chain = 0; chain < nChains; chain++) {
    const samples = runGibbsSampler(
      y, sigma2, muPrior, tauPriorParams,
      nIterations, nBurnin, nThin, rng
    );
    chains.push(samples);
  }

  // Combine chains and calculate summaries
  const allMu = chains.flatMap(c => c.mu);
  const allTau2 = chains.flatMap(c => c.tau2);
  const allTau = allTau2.map(t => Math.sqrt(t));

  // Posterior summaries
  const alpha = 1 - confLevel;

  const muSummary = summarizePosterior(allMu, alpha);
  const tau2Summary = summarizePosterior(allTau2, alpha);
  const tauSummary = summarizePosterior(allTau, alpha);

  // Study-specific shrinkage estimates
  const theta = active.map((e, i) => {
    // Posterior mean of θi given μ and τ²
    // E[θi | y, μ, τ²] = (yi/σi² + μ/τ²) / (1/σi² + 1/τ²)
    const shrinkageEstimates = [];

    for (let s = 0; s < allMu.length; s++) {
      const mu = allMu[s];
      const tau2 = allTau2[s];

      if (tau2 > 0) {
        const w1 = 1 / sigma2[i];
        const w2 = 1 / tau2;
        const thetaS = (y[i] * w1 + mu * w2) / (w1 + w2);
        shrinkageEstimates.push(thetaS);
      } else {
        shrinkageEstimates.push(mu);
      }
    }

    return summarizePosterior(shrinkageEstimates, alpha);
  });

  // Prediction interval for new study
  const predictive = [];
  for (let s = 0; s < allMu.length; s++) {
    const mu = allMu[s];
    const tau2 = allTau2[s];
    // Sample from predictive: θ_new ~ N(μ, τ²)
    predictive.push(mu + Math.sqrt(tau2) * rng.normal());
  }
  const predictiveSummary = summarizePosterior(predictive, alpha);

  // Model diagnostics
  const diagnostics = calculateDiagnostics(chains, allMu, allTau2);

  // Probability statements
  const probabilities = {
    // P(μ > 0) - probability of positive effect
    muPositive: allMu.filter(m => m > 0).length / allMu.length,
    // P(μ < 0) - probability of negative effect
    muNegative: allMu.filter(m => m < 0).length / allMu.length,
    // P(τ < 0.1) - probability of low heterogeneity
    tauLow: allTau.filter(t => t < 0.1).length / allTau.length,
    // P(τ > 0.5) - probability of substantial heterogeneity
    tauHigh: allTau.filter(t => t > 0.5).length / allTau.length
  };

  // Calculate DIC (Deviance Information Criterion)
  const dic = calculateDIC(y, sigma2, allMu, allTau2);

  return {
    method: 'Bayesian Random Effects',

    // Posterior summaries
    mu: muSummary,
    tau2: tau2Summary,
    tau: tauSummary,
    predictive: predictiveSummary,

    // Study-specific estimates
    studies: active.map((e, i) => ({
      id: e.id,
      observed: { es: e.es, se: Math.sqrt(e.vi) },
      posterior: theta[i],
      shrinkage: 1 - (theta[i].sd / Math.sqrt(e.vi))
    })),

    // Probability statements
    probabilities,

    // Model comparison
    dic,

    // Diagnostics
    diagnostics,

    // Settings
    settings: {
      k,
      nIterations,
      nBurnin,
      nThin,
      nChains,
      muPrior,
      tauPrior: tauPriorParams,
      confLevel
    },

    // Raw samples (for further analysis)
    samples: {
      mu: allMu,
      tau2: allTau2
    }
  };
}

/**
 * Gibbs sampler for random-effects model
 */
function runGibbsSampler(y, sigma2, muPrior, tauPrior, nIter, nBurnin, thin, rng) {
  const k = y.length;

  // Initialize
  let mu = muPrior.mean;
  let tau2 = 0.1;

  const samples = { mu: [], tau2: [] };

  for (let iter = 0; iter < nIter + nBurnin; iter++) {
    // Sample μ | τ², y
    // Posterior: N(μ*, V*)
    // V* = 1 / (1/σ0² + Σ 1/(σi² + τ²))
    // μ* = V* * (μ0/σ0² + Σ yi/(σi² + τ²))

    const priorPrec = 1 / (muPrior.sd * muPrior.sd);
    let postPrec = priorPrec;
    let postMean = muPrior.mean * priorPrec;

    for (let i = 0; i < k; i++) {
      const w = 1 / (sigma2[i] + tau2);
      postPrec += w;
      postMean += y[i] * w;
    }

    postMean /= postPrec;
    const postSD = Math.sqrt(1 / postPrec);
    mu = postMean + postSD * rng.normal();

    // Sample τ² | μ, y
    // Depends on prior type
    tau2 = sampleTau2(y, sigma2, mu, tau2, tauPrior, rng);

    // Store sample (after burnin, with thinning)
    if (iter >= nBurnin && (iter - nBurnin) % thin === 0) {
      samples.mu.push(mu);
      samples.tau2.push(tau2);
    }
  }

  return samples;
}

/**
 * Sample τ² using Metropolis-Hastings
 */
function sampleTau2(y, sigma2, mu, tau2Current, tauPrior, rng) {
  const k = y.length;

  // Proposal: log-normal random walk
  const proposalSD = 0.5;
  const logTau2Prop = Math.log(Math.max(1e-10, tau2Current)) + proposalSD * rng.normal();
  const tau2Prop = Math.exp(logTau2Prop);

  // Log-likelihood
  const logLikCurrent = logLikelihood(y, sigma2, mu, tau2Current);
  const logLikProp = logLikelihood(y, sigma2, mu, tau2Prop);

  // Log-prior
  const logPriorCurrent = logTau2Prior(tau2Current, tauPrior);
  const logPriorProp = logTau2Prior(tau2Prop, tauPrior);

  // Jacobian for log-transform
  const logJacobianCurrent = Math.log(tau2Current);
  const logJacobianProp = Math.log(tau2Prop);

  // Acceptance ratio
  const logAlpha = (logLikProp + logPriorProp + logJacobianProp) -
                   (logLikCurrent + logPriorCurrent + logJacobianCurrent);

  if (Math.log(rng.uniform()) < logAlpha) {
    return tau2Prop;
  }
  return tau2Current;
}

/**
 * Log-likelihood for random effects model
 */
function logLikelihood(y, sigma2, mu, tau2) {
  let ll = 0;
  for (let i = 0; i < y.length; i++) {
    const v = sigma2[i] + tau2;
    ll -= 0.5 * Math.log(2 * Math.PI * v);
    ll -= 0.5 * (y[i] - mu) ** 2 / v;
  }
  return ll;
}

/**
 * Log-prior for τ²
 */
function logTau2Prior(tau2, prior) {
  if (tau2 <= 0) return -Infinity;

  switch (prior.type) {
    case 'half-cauchy': {
      // Half-Cauchy(0, scale) for τ (not τ²)
      const tau = Math.sqrt(tau2);
      const scale = prior.scale || 0.5;
      // p(τ) ∝ 1/(1 + (τ/scale)²)
      // Jacobian for τ → τ²: 1/(2τ)
      return -Math.log(1 + (tau / scale) ** 2) - 0.5 * Math.log(tau2);
    }

    case 'inverse-gamma': {
      // InverseGamma(shape, rate)
      const a = prior.shape || 0.001;
      const b = prior.rate || 0.001;
      return a * Math.log(b) - logGamma(a) - (a + 1) * Math.log(tau2) - b / tau2;
    }

    case 'uniform': {
      // Uniform(0, max) for τ
      const max = prior.max || 2;
      const tau = Math.sqrt(tau2);
      if (tau > max) return -Infinity;
      return -0.5 * Math.log(tau2) - Math.log(max);  // Jacobian included
    }

    case 'log-normal': {
      // Log-normal for τ²
      const logMean = prior.logMean || -2;
      const logSD = prior.logSD || 1;
      const logTau2 = Math.log(tau2);
      return -0.5 * ((logTau2 - logMean) / logSD) ** 2 - logTau2 - Math.log(logSD);
    }

    default:
      // Improper uniform (reference prior)
      return 0;
  }
}

/**
 * Get informative prior based on Turner et al. (2012) taxonomy
 */
function getInformativePrior(category) {
  // Log-normal priors for τ² based on empirical meta-epidemiological research
  const priors = {
    // Pharmacological vs placebo
    'pharmacological-placebo': { type: 'log-normal', logMean: -2.56, logSD: 1.74 },
    // Pharmacological vs active control
    'pharmacological-active': { type: 'log-normal', logMean: -2.34, logSD: 1.62 },
    // Non-pharmacological vs inactive control
    'non-pharmacological-inactive': { type: 'log-normal', logMean: -1.81, logSD: 1.51 },
    // Non-pharmacological vs active control
    'non-pharmacological-active': { type: 'log-normal', logMean: -1.67, logSD: 1.39 },
    // Surgery
    'surgery': { type: 'log-normal', logMean: -2.01, logSD: 1.64 },
    // All-cause mortality (Rhodes et al.)
    'mortality': { type: 'log-normal', logMean: -4.18, logSD: 1.41 },
    // Semi-objective outcomes
    'semi-objective': { type: 'log-normal', logMean: -2.13, logSD: 1.58 },
    // Subjective outcomes
    'subjective': { type: 'log-normal', logMean: -1.38, logSD: 1.48 }
  };

  return priors[category] || { type: 'half-cauchy', scale: 0.5 };
}

/**
 * Summarize posterior distribution
 */
function summarizePosterior(samples, alpha) {
  const sorted = [...samples].sort((a, b) => a - b);
  const n = sorted.length;

  const mean = samples.reduce((a, b) => a + b, 0) / n;
  const variance = samples.reduce((acc, x) => acc + (x - mean) ** 2, 0) / (n - 1);
  const sd = Math.sqrt(variance);

  // Quantiles
  const q = (p) => {
    const idx = Math.floor(p * (n - 1));
    const frac = p * (n - 1) - idx;
    if (idx >= n - 1) return sorted[n - 1];
    return sorted[idx] * (1 - frac) + sorted[idx + 1] * frac;
  };

  return {
    mean,
    sd,
    median: q(0.5),
    ciLo: q(alpha / 2),
    ciHi: q(1 - alpha / 2),
    q025: q(0.025),
    q975: q(0.975)
  };
}

/**
 * Calculate convergence diagnostics
 */
function calculateDiagnostics(chains, allMu, allTau2) {
  // Gelman-Rubin R-hat (simplified)
  const rhatMu = calculateRhat(chains.map(c => c.mu));
  const rhatTau2 = calculateRhat(chains.map(c => c.tau2));

  // Effective sample size (simplified)
  const essMu = calculateESS(allMu);
  const essTau2 = calculateESS(allTau2);

  return {
    rhat: { mu: rhatMu, tau2: rhatTau2 },
    ess: { mu: essMu, tau2: essTau2 },
    converged: rhatMu < 1.1 && rhatTau2 < 1.1
  };
}

/**
 * Gelman-Rubin R-hat statistic
 */
function calculateRhat(chainSamples) {
  const nChains = chainSamples.length;
  if (nChains < 2) return 1;

  const n = chainSamples[0].length;

  // Within-chain variance
  const chainMeans = chainSamples.map(c => c.reduce((a, b) => a + b, 0) / n);
  const chainVars = chainSamples.map((c, j) =>
    c.reduce((acc, x) => acc + (x - chainMeans[j]) ** 2, 0) / (n - 1)
  );
  const W = chainVars.reduce((a, b) => a + b, 0) / nChains;

  // Between-chain variance
  const overallMean = chainMeans.reduce((a, b) => a + b, 0) / nChains;
  const B = chainMeans.reduce((acc, m) => acc + (m - overallMean) ** 2, 0) * n / (nChains - 1);

  // R-hat
  const varPlus = (n - 1) / n * W + B / n;
  return Math.sqrt(varPlus / W);
}

/**
 * Effective sample size
 */
function calculateESS(samples) {
  const n = samples.length;
  const mean = samples.reduce((a, b) => a + b, 0) / n;

  // Autocorrelation at lag 1
  let rho1 = 0;
  for (let i = 0; i < n - 1; i++) {
    rho1 += (samples[i] - mean) * (samples[i + 1] - mean);
  }
  rho1 /= samples.reduce((acc, x) => acc + (x - mean) ** 2, 0);

  // ESS approximation
  return Math.max(1, n * (1 - rho1) / (1 + rho1));
}

/**
 * Calculate DIC (Deviance Information Criterion)
 */
function calculateDIC(y, sigma2, muSamples, tau2Samples) {
  const n = muSamples.length;

  // Mean deviance
  let Dbar = 0;
  for (let s = 0; s < n; s++) {
    Dbar += -2 * logLikelihood(y, sigma2, muSamples[s], tau2Samples[s]);
  }
  Dbar /= n;

  // Deviance at posterior mean
  const muBar = muSamples.reduce((a, b) => a + b, 0) / n;
  const tau2Bar = tau2Samples.reduce((a, b) => a + b, 0) / n;
  const Dhat = -2 * logLikelihood(y, sigma2, muBar, tau2Bar);

  // Effective number of parameters
  const pD = Dbar - Dhat;

  // DIC
  const DIC = Dbar + pD;

  return { DIC, pD, Dbar, Dhat };
}

/**
 * Simple pseudo-random number generator with seed support
 */
function createRNG(seed = null) {
  let state = seed || Date.now();

  const uniform = () => {
    state = (state * 1103515245 + 12345) & 0x7fffffff;
    return state / 0x7fffffff;
  };

  const normal = () => {
    // Box-Muller transform
    const u1 = uniform();
    const u2 = uniform();
    return Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
  };

  return { uniform, normal };
}

/**
 * Log gamma function (Lanczos approximation)
 */
function logGamma(x) {
  const g = 7;
  const c = [
    0.99999999999980993, 676.5203681218851, -1259.1392167224028,
    771.32342877765313, -176.61502916214059, 12.507343278686905,
    -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7
  ];

  if (x < 0.5) {
    return Math.log(Math.PI / Math.sin(Math.PI * x)) - logGamma(1 - x);
  }

  x -= 1;
  let a = c[0];
  for (let i = 1; i < g + 2; i++) {
    a += c[i] / (x + i);
  }
  const t = x + g + 0.5;

  return 0.5 * Math.log(2 * Math.PI) + (x + 0.5) * Math.log(t) - t + Math.log(a);
}

/**
 * Bayesian model comparison between fixed and random effects
 */
export function bayesianModelComparison(effects, options = {}) {
  // Run random effects model
  const reResult = bayesianMetaAnalysis(effects, { ...options, tauPrior: { type: 'half-cauchy', scale: 0.5 } });

  // Run with very tight prior on tau (approximates fixed effects)
  const feResult = bayesianMetaAnalysis(effects, { ...options, tauPrior: { type: 'half-cauchy', scale: 0.001 } });

  if (reResult.error || feResult.error) {
    return { error: reResult.error || feResult.error };
  }

  // Compare DICs
  const deltaDIC = feResult.dic.DIC - reResult.dic.DIC;

  return {
    randomEffects: {
      mu: reResult.mu,
      tau: reResult.tau,
      dic: reResult.dic
    },
    fixedEffects: {
      mu: feResult.mu,
      dic: feResult.dic
    },
    comparison: {
      deltaDIC,
      preferredModel: deltaDIC > 0 ? 'Random Effects' : 'Fixed Effects',
      evidence: Math.abs(deltaDIC) < 2 ? 'Weak' :
                Math.abs(deltaDIC) < 5 ? 'Moderate' :
                Math.abs(deltaDIC) < 10 ? 'Strong' : 'Very Strong'
    }
  };
}
