/**
 * Statistical Engine - Main exports
 */

// Distribution functions
export {
  normCDF,
  pFromZ,
  tCritical,
  tCDF,
  pFromT,
  chiSqP,
  normQuantile
} from './distributions.js';

// Effect size calculations
export {
  calculateEffect,
  transformEffect,
  isRatioMetric
} from './effects.js';

// Tau² estimators
export {
  tauDL,
  tauREML,
  tauPM,
  tauSJ,
  tauHE,
  tauHS,
  tauEB,
  getTauEstimator
} from './tau.js';

// Pooling methods
export {
  poolInverseVariance,
  leaveOneOut,
  cumulativeMetaAnalysis,
  calculateInfluenceDiagnostics
} from './pooling.js';

// Publication bias tests
export {
  eggerTest,
  beggTest,
  petersTest,
  failSafeN,
  orwinFailSafeN,
  trimAndFill,
  getFunnelPlotData,
  getContourFunnelData
} from './bias.js';

// Trial Sequential Analysis
export {
  calculateTSA,
  getTSAPlotData
} from './tsa.js';
