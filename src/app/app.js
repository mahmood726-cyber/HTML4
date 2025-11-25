/**
 * Main Application Module
 */

import {
  initDOMCache,
  getElement,
  getValue,
  setText,
  setHTML,
  toggleClass,
  setVisible,
  addListener
} from './dom.js';

import {
  getState,
  getStateValue,
  updateState,
  setStateValue,
  exportState,
  importState,
  importCSV,
  exportCSV,
  addStudyRow,
  updateStudyRow,
  clearStudyData,
  subscribe
} from './state.js';

import {
  calculateEffect,
  getTauEstimator,
  poolInverseVariance,
  isRatioMetric
} from '../engine/index.js';

import {
  escapeHTML,
  formatNumber,
  formatPValue,
  formatEffectSize,
  validateStudyRow,
  validateConfig,
  debounce
} from '../utils/index.js';

import {
  DATA_COLUMNS,
  COLUMN_LABELS,
  NUMERIC_COLUMNS,
  RATIO_METRICS
} from '../constants/index.js';

import { renderAll } from '../renderers/index.js';

/**
 * Application controller
 */
export const App = {
  /**
   * Initialize the application
   */
  init() {
    initDOMCache();
    this.setupModules();
    this.setupTabs();
    this.setupTeaching();
    this.setupEventListeners();
    this.changeDataType();

    // Subscribe to state changes
    subscribe(debounce(() => this.onStateChange(), 100));
  },

  /**
   * Set up module navigation
   */
  setupModules() {
    document.querySelectorAll('.module-btn').forEach(btn => {
      btn.addEventListener('click', () => {
        document.querySelectorAll('.module-btn').forEach(b => b.classList.remove('active'));
        document.querySelectorAll('.module-panel').forEach(p => p.classList.remove('active'));
        btn.classList.add('active');

        const moduleId = 'module-' + btn.dataset.module;
        const moduleEl = document.getElementById(moduleId);
        if (moduleEl) {
          moduleEl.classList.add('active');
        }

        setStateValue('ui.activeModule', btn.dataset.module);
      });
    });
  },

  /**
   * Set up tab navigation within modules
   */
  setupTabs() {
    document.querySelectorAll('.tabs').forEach(tabContainer => {
      tabContainer.querySelectorAll('.tab').forEach(tab => {
        tab.addEventListener('click', () => {
          const parent = tab.closest('.module-panel');
          parent.querySelectorAll('.tab').forEach(t => t.classList.remove('active'));
          parent.querySelectorAll('.tab-content').forEach(c => c.classList.remove('active'));
          tab.classList.add('active');

          const tabContent = parent.querySelector('#tab-' + tab.dataset.tab);
          if (tabContent) {
            tabContent.classList.add('active');
          }
        });
      });
    });
  },

  /**
   * Set up teaching mode toggle
   */
  setupTeaching() {
    const teachingCheckbox = getElement('teachingMode');
    if (teachingCheckbox) {
      teachingCheckbox.addEventListener('change', (e) => {
        document.querySelectorAll('.teaching-box').forEach(box => {
          box.classList.toggle('teaching-hidden', !e.target.checked);
        });
        setStateValue('config.teachingMode', e.target.checked);
      });
    }
  },

  /**
   * Set up global event listeners
   */
  setupEventListeners() {
    // Config changes trigger recalculation
    const configInputs = [
      'dataType', 'metric', 'method', 'model', 'tauMethod',
      'confLevel', 'contCorr', 'hksj', 'predictionInterval',
      'subgroupVar', 'robFilter'
    ];

    configInputs.forEach(key => {
      addListener(key, 'change', () => this.recalc());
    });

    // Economic inputs
    ['cer', 'costEvent', 'costTreatment', 'qalyLoss'].forEach(key => {
      addListener(key, 'change', () => this.recalc());
    });

    // TSA inputs
    ['tsaAlpha', 'tsaBeta', 'tsaRrr'].forEach(key => {
      addListener(key, 'change', () => this.recalc());
    });
  },

  /**
   * Handle data type change
   */
  changeDataType() {
    const dataType = getValue('dataType');
    const cols = DATA_COLUMNS[dataType];

    // Update metric options based on data type
    const metricGroup = getElement('metricGroup');
    const metricSelect = getElement('metric');
    const methodSelect = getElement('method');

    if (dataType === 'binary') {
      setVisible('metricGroup', true);
      if (metricSelect) {
        metricSelect.innerHTML = `
          <option value="OR">Odds Ratio (OR)</option>
          <option value="RR">Risk Ratio (RR)</option>
          <option value="RD">Risk Difference (RD)</option>
        `;
      }
      if (methodSelect) methodSelect.disabled = false;
    } else if (dataType === 'survival') {
      setVisible('metricGroup', true);
      if (metricSelect) {
        metricSelect.innerHTML = `<option value="HR">Hazard Ratio (HR)</option>`;
      }
      if (methodSelect) {
        methodSelect.value = 'iv';
        methodSelect.disabled = true;
      }
    } else {
      setVisible('metricGroup', false);
      if (methodSelect) {
        methodSelect.value = 'iv';
        methodSelect.disabled = true;
      }
    }

    // Update table headers
    this.renderTableHeaders(cols);

    // Initialize data if empty
    const state = getState();
    if (state.data.length === 0) {
      const emptyRows = Array(8).fill(0).map(() => {
        const row = {};
        cols.forEach(c => row[c] = '');
        return row;
      });
      updateState({ data: emptyRows, config: { ...state.config, dataType } });
    } else {
      setStateValue('config.dataType', dataType);
    }

    this.renderTable();
  },

  /**
   * Render table headers
   * @param {Array} cols - Column names
   */
  renderTableHeaders(cols) {
    const thead = getElement('tableHead');
    if (!thead) return;

    thead.innerHTML = '<tr>' + cols.map(c => {
      const isNumeric = NUMERIC_COLUMNS.includes(c);
      return `<th class="${isNumeric ? 'num' : ''}">${escapeHTML(COLUMN_LABELS[c] || c)}</th>`;
    }).join('') + '</tr>';
  },

  /**
   * Render the data table
   */
  renderTable() {
    const state = getState();
    const dataType = state.config.dataType;
    const cols = DATA_COLUMNS[dataType];
    const tbody = getElement('tableBody');

    if (!tbody) return;

    tbody.innerHTML = state.data.map((row, i) => {
      const excluded = row.exclude === true || row.exclude === 'true';
      return `<tr class="${excluded ? 'excluded' : ''}" data-row="${i}">
        ${cols.map(c => this.renderTableCell(row, c, i)).join('')}
      </tr>`;
    }).join('');

    setText('studyCount', state.data.filter(r => r.study).length.toString());
  },

  /**
   * Render a single table cell
   * @param {Object} row - Row data
   * @param {string} col - Column name
   * @param {number} rowIndex - Row index
   * @returns {string} Cell HTML
   */
  renderTableCell(row, col, rowIndex) {
    if (col === 'exclude') {
      const checked = row[col] === true || row[col] === 'true';
      return `<td><input type="checkbox" ${checked ? 'checked' : ''}
        onchange="window.App.updateCell(${rowIndex},'${col}',this.checked)"
        aria-label="Exclude study"></td>`;
    }

    if (col === 'rob') {
      return `<td><select onchange="window.App.updateCell(${rowIndex},'${col}',this.value)"
        style="width:60px" aria-label="Risk of Bias">
        <option value="">—</option>
        <option value="low" ${row[col] === 'low' ? 'selected' : ''}>Low</option>
        <option value="some" ${row[col] === 'some' ? 'selected' : ''}>Some</option>
        <option value="high" ${row[col] === 'high' ? 'selected' : ''}>High</option>
      </select></td>`;
    }

    const isNumeric = NUMERIC_COLUMNS.includes(col);
    const inputType = isNumeric ? 'number' : 'text';
    const value = escapeHTML(row[col] || '');

    return `<td><input type="${inputType}"
      value="${value}"
      ${isNumeric ? 'step="any"' : ''}
      onchange="window.App.updateCell(${rowIndex},'${col}',this.value)"
      aria-label="${COLUMN_LABELS[col] || col}"></td>`;
  },

  /**
   * Update a cell value
   * @param {number} rowIndex - Row index
   * @param {string} col - Column name
   * @param {*} value - New value
   */
  updateCell(rowIndex, col, value) {
    updateStudyRow(rowIndex, col, value);

    // Update visual state for exclusion
    if (col === 'exclude') {
      const tr = document.querySelector(`tr[data-row="${rowIndex}"]`);
      if (tr) {
        tr.classList.toggle('excluded', Boolean(value));
      }
    }

    this.recalc();
  },

  /**
   * Add a new row to the table
   */
  addRow() {
    const state = getState();
    const cols = DATA_COLUMNS[state.config.dataType];
    const row = {};
    cols.forEach(c => row[c] = '');
    addStudyRow(row);
    this.renderTable();
  },

  /**
   * Clear all data
   */
  clearData() {
    if (confirm('Clear all study data?')) {
      clearStudyData();
      this.changeDataType();
    }
  },

  /**
   * Main recalculation function
   */
  recalc() {
    const state = getState();
    const { dataType } = state.config;
    const metric = getValue('metric') || 'OR';
    const tauMethod = getValue('tauMethod') || 'dl';
    const useHKSJ = getValue('hksj');
    const confLevel = parseFloat(getValue('confLevel')) || 0.95;

    // Validate configuration
    const configValidation = validateConfig({ confLevel });
    if (!configValidation.valid) {
      console.warn('Invalid config:', configValidation.errors);
    }

    // Calculate effects for each study
    const effects = state.data.map((row, i) => {
      const eff = calculateEffect(row, dataType, metric);
      if (!eff) return null;

      return {
        id: row.study || `Study ${i + 1}`,
        year: row.year || '',
        ...eff,
        excluded: row.exclude === true || row.exclude === 'true',
        rob: row.rob,
        subgroup: row.subgroup
      };
    }).filter(e => e !== null);

    if (effects.length === 0) {
      this.clearResults();
      return;
    }

    // Get active (non-excluded) effects
    const activeEffects = effects.filter(e => !e.excluded);

    // Estimate tau²
    let tau2 = 0;
    if (activeEffects.length >= 2) {
      const tauEstimator = getTauEstimator(tauMethod);
      tau2 = tauEstimator(activeEffects);
    }

    // Pool effects
    const result = poolInverseVariance(effects, tau2, {
      useHKSJ,
      confLevel,
      showPredictionInterval: getValue('predictionInterval')
    });

    if (!result) {
      this.clearResults();
      return;
    }

    // Store result with additional metadata
    const enrichedResult = {
      ...result,
      metric,
      type: dataType,
      isRatio: isRatioMetric(metric) || dataType === 'survival'
    };

    updateState({ currentResult: enrichedResult });
    this.updateDisplay(enrichedResult);
  },

  /**
   * Clear all results display
   */
  clearResults() {
    const ids = [
      'resEs', 'resCi', 'resP', 'resZ', 'resI2', 'resTau',
      'resQ', 'resQp', 'resK', 'resN', 'resPi', 'resH2',
      'resH2Ci', 'resTauVal', 'resTauCi'
    ];
    ids.forEach(id => setText(id, '—'));
    updateState({ currentResult: null });
  },

  /**
   * Update all display elements
   * @param {Object} result - Pooled result
   */
  updateDisplay(result) {
    const { isRatio } = result;

    // Update summary statistics
    setText('resEs', formatEffectSize(result.es, isRatio));
    setText('resCi', `95% CI: ${formatEffectSize(result.ciLo, isRatio)} – ${formatEffectSize(result.ciHi, isRatio)}`);

    const pEl = getElement('resP');
    if (pEl) {
      pEl.textContent = formatPValue(result.pVal);
      pEl.classList.toggle('sig', result.pVal < 0.05);
    }

    setText('resZ', `z = ${formatNumber(result.z)}`);
    setText('resI2', `${formatNumber(result.I2, 1)}%`);
    setText('resTau', `τ² = ${formatNumber(result.tau2, 4)}`);
    setText('resQ', formatNumber(result.Q));
    setText('resQp', `df = ${result.df}, p = ${formatPValue(result.qPVal)}`);
    setText('resK', result.k.toString());
    setText('resN', result.totalN ? `N = ${result.totalN}` : '—');
    setText('resPi', result.piLo !== null
      ? `${formatEffectSize(result.piLo, isRatio)} – ${formatEffectSize(result.piHi, isRatio)}`
      : '—');
    setText('resH2', formatNumber(result.H2));
    setText('resTauVal', formatNumber(result.tau, 4));

    // Render all visualizations
    renderAll(result, getState());
  },

  /**
   * Handle state change
   */
  onStateChange() {
    this.renderTable();
  },

  /**
   * Export state to JSON file
   */
  exportState() {
    const json = exportState();
    this.downloadFile(json, 'meta-analysis-data.json', 'application/json');
  },

  /**
   * Export data to CSV file
   */
  exportCSV() {
    const csv = exportCSV();
    this.downloadFile(csv, 'meta-analysis-data.csv', 'text/csv');
  },

  /**
   * Download a file
   * @param {string} content - File content
   * @param {string} filename - File name
   * @param {string} mimeType - MIME type
   */
  downloadFile(content, filename, mimeType) {
    const blob = new Blob([content], { type: mimeType });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = filename;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
  },

  /**
   * Import file
   * @param {Event} event - File input event
   */
  importFile(event) {
    const file = event.target.files[0];
    if (!file) return;

    const reader = new FileReader();
    reader.onload = (e) => {
      const content = e.target.result;

      if (file.name.endsWith('.json')) {
        if (importState(content)) {
          this.changeDataType();
          this.recalc();
        } else {
          alert('Failed to import JSON file');
        }
      } else if (file.name.endsWith('.csv')) {
        if (importCSV(content)) {
          this.changeDataType();
          this.recalc();
        } else {
          alert('Failed to import CSV file');
        }
      }
    };

    reader.readAsText(file);
    event.target.value = ''; // Reset input
  },

  /**
   * Load example data (BCG vaccine meta-analysis)
   */
  loadExample() {
    const bcgData = [
      { study: 'Aronson', year: 1948, e1: 4, n1: 123, e2: 11, n2: 139, rob: 'some' },
      { study: 'Ferguson & Simes', year: 1949, e1: 6, n1: 306, e2: 29, n2: 303, rob: 'low' },
      { study: 'Rosenthal', year: 1960, e1: 3, n1: 231, e2: 11, n2: 220, rob: 'low' },
      { study: 'Hart & Sutherland', year: 1977, e1: 62, n1: 13598, e2: 248, n2: 12867, rob: 'low' },
      { study: 'Frimodt-Moller', year: 1973, e1: 33, n1: 5069, e2: 47, n2: 5808, rob: 'some' },
      { study: 'Stein & Aronson', year: 1953, e1: 180, n1: 1541, e2: 372, n2: 1451, rob: 'high' },
      { study: 'Vandiviere', year: 1973, e1: 8, n1: 2545, e2: 10, n2: 629, rob: 'some' },
      { study: 'TPT Madras', year: 1980, e1: 505, n1: 88391, e2: 499, n2: 88391, rob: 'low' },
      { study: 'Coetzee & Berjak', year: 1968, e1: 29, n1: 7499, e2: 45, n2: 7277, rob: 'some' },
      { study: 'Rosenthal (Chicago)', year: 1961, e1: 17, n1: 1716, e2: 65, n2: 1665, rob: 'low' },
      { study: 'Comstock (Georgia)', year: 1974, e1: 186, n1: 50634, e2: 141, n2: 27338, rob: 'low' },
      { study: 'Comstock (P.R.)', year: 1969, e1: 141, n1: 27338, e2: 186, n2: 50634, rob: 'low' },
      { study: 'Comstock & Webster', year: 1969, e1: 5, n1: 2498, e2: 3, n2: 2341, rob: 'some' }
    ];

    updateState({ data: bcgData, config: { ...getState().config, dataType: 'binary' } });
    this.changeDataType();
    this.recalc();
  },

  /**
   * Run GOSH analysis (moved to separate function for cancellation support)
   */
  async runGOSH() {
    // Implementation moved to renderers module
    const { runGOSHAnalysis } = await import('../renderers/heterogeneity.js');
    runGOSHAnalysis(getState().currentResult);
  },

  /**
   * Run validation suite
   */
  runValidation() {
    // Implementation would go here
    alert('Validation suite not yet implemented in modular version');
  },

  /**
   * Paste data from clipboard
   */
  async pasteFromClipboard() {
    try {
      const text = await navigator.clipboard.readText();
      if (importCSV(text)) {
        this.changeDataType();
        this.recalc();
      }
    } catch (err) {
      console.error('Failed to read clipboard:', err);
      alert('Failed to paste from clipboard. Please use the Import button instead.');
    }
  }
};

// Make App available globally for inline event handlers
window.App = App;
