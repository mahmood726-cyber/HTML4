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

// Mantel-Haenszel and Peto methods (sparse data)
export {
  mantelHaenszel,
  petoOddsRatio,
  recommendBinaryMethod
} from './mhPeto.js';

// Meta-regression
export {
  univariateMetaRegression,
  multivariateMetaRegression,
  subgroupMetaRegression,
  createBubblePlotData
} from './metaRegression.js';

// Network Meta-Analysis
export {
  createNetwork,
  isNetworkConnected,
  bucherIndirect,
  allBucherComparisons,
  frequentistNMA,
  nodeSplitting,
  globalInconsistencyTest,
  generateLeagueTable,
  getNetworkDiagramData
} from './nma.js';

// MAIC and Population-Adjusted Methods
export {
  performMAIC,
  anchoredMAIC,
  unanchoredMAIC,
  simulatedTreatmentComparison
} from './maic.js';

// Component NMA
export {
  additiveCNMA,
  interactionCNMA,
  prepareCNMAUpsetData
} from './componentNMA.js';

// Bayesian Meta-Analysis
export {
  bayesianMetaAnalysis,
  bayesianModelComparison
} from './bayesian.js';

// Survival Analysis and IPD Reconstruction
export {
  reconstructIPD,
  calculateKMFromIPD,
  calculateRMST,
  extractHazardRatio,
  fitParametricSurvival,
  compareSurvivalDistributions
} from './survival.js';

// PET-PEESE Bias Correction
export {
  runPET,
  runPEESE,
  petPeese,
  selectionModel,
  enhancedTrimFill,
  comprehensiveBiasAssessment
} from './petPeese.js';

// Transportability and Generalizability
export {
  calculatePopulationSimilarity,
  inverseOddsWeighting,
  assessStudyTransportability,
  calibrationWeighting,
  calculateGeneralizabilityIndex,
  estimateTargetPopulationEffect
} from './transportability.js';

// Diagnostic Test Accuracy (DTA) Meta-Analysis
export {
  calculateAccuracyMeasures,
  bivariateModel,
  generateSROCCurve,
  hsrocModel,
  dtaMetaRegression,
  quadasWeightedAnalysis
} from './dta.js';

// Data Conversion Utilities
export {
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
} from './dataConversion.js';

// Dose-Response Meta-Analysis
export {
  calculateDoseCovariance,
  linearDoseResponse,
  quadraticDoseResponse,
  splineDoseResponse,
  compareDoseResponseModels,
  getDoseResponsePlotData
} from './doseResponse.js';

// Cross-Validation and Overfitting Detection
export {
  leaveOneOutCV,
  kFoldCV,
  influenceDiagnostics,
  predictionInterval,
  assessOverfittingRisk,
  bootstrapValidation
} from './crossValidation.js';
