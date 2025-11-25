/**
 * DOM element caching and management
 */

/**
 * Cached DOM element references
 */
let cachedElements = null;

/**
 * DOM element ID mappings
 */
const ELEMENT_IDS = {
  // Data inputs
  dataType: 'data-type',
  metric: 'metric',
  metricGroup: 'metric-group',
  method: 'method',
  model: 'model',
  tauMethod: 'tau-method',
  confLevel: 'conf-level',
  contCorr: 'cont-corr',
  hksj: 'hksj',
  predictionInterval: 'prediction-interval',
  teachingMode: 'teaching-mode',
  subgroupVar: 'subgroup-var',
  robFilter: 'rob-filter',

  // Economic inputs
  cer: 'cer',
  costEvent: 'cost-event',
  costTreatment: 'cost-treatment',
  qalyLoss: 'qaly-loss',

  // TSA inputs
  tsaAlpha: 'tsa-alpha',
  tsaBeta: 'tsa-beta',
  tsaRrr: 'tsa-rrr',

  // Table elements
  tableHead: 'table-head',
  tableBody: 'table-body',
  studyCount: 'study-count',

  // Results display
  resEs: 'res-es',
  resCi: 'res-ci',
  resP: 'res-p',
  resZ: 'res-z',
  resI2: 'res-i2',
  resTau: 'res-tau',
  resQ: 'res-q',
  resQp: 'res-qp',
  resK: 'res-k',
  resN: 'res-n',
  resPi: 'res-pi',
  resH2: 'res-h2',
  resH2Ci: 'res-h2ci',
  resTauVal: 'res-tau-val',
  resTauCi: 'res-tau-ci',

  // Bias results
  eggerP: 'egger-p',
  eggerInt: 'egger-int',
  beggP: 'begg-p',
  beggTau: 'begg-tau',
  fsnVal: 'fsn-val',
  tfMissing: 'tf-missing',
  tfAdjusted: 'tf-adjusted',

  // TSA results
  tsaRis: 'tsa-ris',
  tsaAccrued: 'tsa-accrued',
  tsaPct: 'tsa-pct',
  tsaD2: 'tsa-d2',
  tsaConclusion: 'tsa-conclusion',
  tsaBoundary: 'tsa-boundary',

  // Economic results
  nntValue: 'nnt-value',
  nntLabel: 'nnt-label',
  arrValue: 'arr-value',
  arrCi: 'arr-ci',
  costAvoided: 'cost-avoided',
  icerValue: 'icer-value',
  icerThreshold: 'icer-threshold',
  nntGrid: 'nnt-grid',
  nntDisplay: 'nnt-display',
  nntType: 'nnt-type',
  nntDesc: 'nnt-desc',

  // Plot containers
  forestPlot: 'forest-plot',
  funnelPlot: 'funnel-plot',
  contourPlot: 'contour-plot',
  sunsetPlot: 'sunset-plot',
  baujatPlot: 'baujat-plot',
  labbePlot: 'labbe-plot',
  radialPlot: 'radial-plot',
  goshPlot: 'gosh-plot',
  tsaPlot: 'tsa-plot',
  looPlot: 'loo-plot',
  cumulativePlot: 'cumulative-plot',
  influencePlot: 'influence-plot',
  robPlot: 'rob-plot',

  // Other containers
  studyEffectsTable: 'study-effects-table',
  influenceStats: 'influence-stats',
  gradeProfile: 'grade-profile',
  sofTable: 'sof-table',
  robSummary: 'rob-summary',
  validationResults: 'validation-results',

  // Import
  importFile: 'import-file'
};

/**
 * Initialize and cache all DOM elements
 * @returns {Object} Cached element references
 */
export function initDOMCache() {
  cachedElements = {};

  for (const [key, id] of Object.entries(ELEMENT_IDS)) {
    const element = document.getElementById(id);
    if (element) {
      cachedElements[key] = element;
    } else {
      console.warn(`DOM element not found: #${id}`);
    }
  }

  return cachedElements;
}

/**
 * Get a cached DOM element
 * @param {string} key - Element key
 * @returns {HTMLElement|null} Cached element
 */
export function getElement(key) {
  if (!cachedElements) {
    initDOMCache();
  }
  return cachedElements[key] || null;
}

/**
 * Get multiple cached DOM elements
 * @param {...string} keys - Element keys
 * @returns {Object} Object with requested elements
 */
export function getElements(...keys) {
  if (!cachedElements) {
    initDOMCache();
  }
  return keys.reduce((acc, key) => {
    acc[key] = cachedElements[key] || null;
    return acc;
  }, {});
}

/**
 * Get value from an input element
 * @param {string} key - Element key
 * @returns {string} Input value
 */
export function getValue(key) {
  const el = getElement(key);
  if (!el) return '';

  if (el.type === 'checkbox') {
    return el.checked;
  }
  return el.value;
}

/**
 * Set value on an input element
 * @param {string} key - Element key
 * @param {*} value - Value to set
 */
export function setValue(key, value) {
  const el = getElement(key);
  if (!el) return;

  if (el.type === 'checkbox') {
    el.checked = Boolean(value);
  } else {
    el.value = value;
  }
}

/**
 * Set text content of an element
 * @param {string} key - Element key
 * @param {string} text - Text content
 */
export function setText(key, text) {
  const el = getElement(key);
  if (el) {
    el.textContent = text;
  }
}

/**
 * Set innerHTML of an element (use with caution - ensure content is sanitized)
 * @param {string} key - Element key
 * @param {string} html - HTML content
 */
export function setHTML(key, html) {
  const el = getElement(key);
  if (el) {
    el.innerHTML = html;
  }
}

/**
 * Add/remove class from element
 * @param {string} key - Element key
 * @param {string} className - Class name
 * @param {boolean} add - Add or remove
 */
export function toggleClass(key, className, add) {
  const el = getElement(key);
  if (el) {
    el.classList.toggle(className, add);
  }
}

/**
 * Set element display style
 * @param {string} key - Element key
 * @param {boolean} visible - Show or hide
 */
export function setVisible(key, visible) {
  const el = getElement(key);
  if (el) {
    el.style.display = visible ? '' : 'none';
  }
}

/**
 * Add event listener to cached element
 * @param {string} key - Element key
 * @param {string} event - Event type
 * @param {Function} handler - Event handler
 */
export function addListener(key, event, handler) {
  const el = getElement(key);
  if (el) {
    el.addEventListener(event, handler);
  }
}

/**
 * Query selector within a cached element
 * @param {string} key - Element key
 * @param {string} selector - CSS selector
 * @returns {HTMLElement|null} Found element
 */
export function queryWithin(key, selector) {
  const el = getElement(key);
  return el ? el.querySelector(selector) : null;
}

/**
 * Query all within a cached element
 * @param {string} key - Element key
 * @param {string} selector - CSS selector
 * @returns {NodeList} Found elements
 */
export function queryAllWithin(key, selector) {
  const el = getElement(key);
  return el ? el.querySelectorAll(selector) : [];
}
