/**
 * Statistical constants and lookup tables
 */

// T-distribution critical values for common df (alpha = 0.05, two-tailed)
export const T_CRIT_TABLE = {
  1: 12.706,
  2: 4.303,
  3: 3.182,
  4: 2.776,
  5: 2.571,
  6: 2.447,
  7: 2.365,
  8: 2.306,
  9: 2.262,
  10: 2.228,
  15: 2.131,
  20: 2.086,
  25: 2.060,
  30: 2.042,
  40: 2.021,
  50: 2.009,
  60: 2.000,
  80: 1.990,
  100: 1.984,
  Infinity: 1.96
};

// Normal distribution approximation constants (Abramowitz & Stegun)
export const NORM_CDF_COEFFICIENTS = {
  a1: 0.254829592,
  a2: -0.284496736,
  a3: 1.421413741,
  a4: -1.453152027,
  a5: 1.061405429,
  p: 0.3275911
};

// Default continuity correction
export const DEFAULT_CONTINUITY_CORRECTION = 0.5;

// REML iteration parameters
export const REML_CONFIG = {
  maxIterations: 50,
  tolerance: 1e-6
};

// Paule-Mandel iteration parameters
export const PM_CONFIG = {
  maxIterations: 100,
  tolerance: 1e-6
};

// GOSH analysis limits
export const GOSH_CONFIG = {
  minStudies: 3,
  maxStudies: 15,
  batchSize: 100
};

// Default confidence level
export const DEFAULT_CONF_LEVEL = 0.95;

// Z-scores for common alpha levels
export const Z_SCORES = {
  0.10: 1.645,
  0.05: 1.96,
  0.01: 2.576
};

// GRADE certainty levels
export const GRADE_LEVELS = {
  HIGH: 'High',
  MODERATE: 'Moderate',
  LOW: 'Low',
  VERY_LOW: 'Very Low'
};

// Risk of bias thresholds
export const ROB_THRESHOLDS = {
  highRiskRatio: 0.3 // Downgrade if >30% of studies are high risk
};

// Heterogeneity thresholds (I²)
export const I2_THRESHOLDS = {
  low: 25,
  moderate: 50,
  high: 75
};

// NICE willingness-to-pay threshold (GBP per QALY)
export const WTP_THRESHOLD = 30000;

// Column definitions for data types
export const DATA_COLUMNS = {
  binary: ['study', 'year', 'e1', 'n1', 'e2', 'n2', 'rob', 'subgroup', 'exclude'],
  continuous: ['study', 'year', 'm1', 's1', 'n1', 'm2', 's2', 'n2', 'rob', 'subgroup', 'exclude'],
  proportion: ['study', 'year', 'e', 'n', 'rob', 'subgroup', 'exclude'],
  survival: ['study', 'year', 'hr', 'll', 'ul', 'rob', 'subgroup', 'exclude'],
  correlation: ['study', 'year', 'r', 'n', 'rob', 'subgroup', 'exclude'],
  generic: ['study', 'year', 'es', 'se', 'rob', 'subgroup', 'exclude']
};

// Column labels for display
export const COLUMN_LABELS = {
  study: 'Study',
  year: 'Year',
  e1: 'E\u2081',
  n1: 'N\u2081',
  e2: 'E\u2082',
  n2: 'N\u2082',
  m1: 'M\u2081',
  s1: 'SD\u2081',
  m2: 'M\u2082',
  s2: 'SD\u2082',
  e: 'Events',
  n: 'N',
  hr: 'HR',
  ll: 'CI Low',
  ul: 'CI High',
  r: 'r',
  es: 'Effect',
  se: 'SE',
  rob: 'RoB',
  subgroup: 'Group',
  exclude: '\u2717'
};

// Numeric columns (for input type determination)
export const NUMERIC_COLUMNS = [
  'e1', 'n1', 'e2', 'n2', 'm1', 's1', 'm2', 's2',
  'n', 'e', 'hr', 'll', 'ul', 'es', 'se', 'r', 'year'
];

// Ratio-based effect measures
export const RATIO_METRICS = ['OR', 'RR', 'HR'];

// Plotly common layout settings
export const PLOTLY_LAYOUT_DEFAULTS = {
  margin: { t: 30, r: 30, b: 50, l: 60 },
  paper_bgcolor: 'rgba(0,0,0,0)',
  plot_bgcolor: 'rgba(250,249,247,1)'
};

// Plotly color palette
export const COLORS = {
  primary: '#0d7377',
  accent: '#c43e3e',
  muted: '#94a3b8',
  success: '#22c55e',
  warning: '#ef4444',
  info: '#3b82f6'
};
