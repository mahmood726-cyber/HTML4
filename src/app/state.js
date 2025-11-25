/**
 * Application state management
 */

import { deepClone } from '../utils/index.js';

/**
 * Default application state
 */
const DEFAULT_STATE = {
  data: [],
  currentResult: null,
  config: {
    dataType: 'binary',
    metric: 'OR',
    method: 'iv',
    model: 're',
    tauMethod: 'dl',
    confLevel: 0.95,
    contCorr: 0.5,
    useHKSJ: true,
    showPredictionInterval: true,
    teachingMode: false
  },
  economic: {
    cer: 0.1,
    costEvent: 10000,
    costTreatment: 500,
    qalyLoss: 0.5
  },
  tsa: {
    alpha: 0.05,
    beta: 0.20,
    rrr: 20
  },
  ui: {
    activeModule: 'data',
    activeTabs: {}
  }
};

/**
 * Current application state
 */
let state = deepClone(DEFAULT_STATE);

/**
 * State change subscribers
 */
const subscribers = new Set();

/**
 * Get current state (immutable copy)
 * @returns {Object} Current state
 */
export function getState() {
  return deepClone(state);
}

/**
 * Get a specific part of state
 * @param {string} path - Dot-notation path (e.g., 'config.dataType')
 * @returns {*} Value at path
 */
export function getStateValue(path) {
  const parts = path.split('.');
  let value = state;

  for (const part of parts) {
    if (value === null || value === undefined) return undefined;
    value = value[part];
  }

  return deepClone(value);
}

/**
 * Update state
 * @param {Object} updates - Partial state updates
 */
export function updateState(updates) {
  state = mergeDeep(state, updates);
  notifySubscribers();
}

/**
 * Set a specific state value
 * @param {string} path - Dot-notation path
 * @param {*} value - New value
 */
export function setStateValue(path, value) {
  const parts = path.split('.');
  const lastPart = parts.pop();
  let target = state;

  for (const part of parts) {
    if (target[part] === undefined) {
      target[part] = {};
    }
    target = target[part];
  }

  target[lastPart] = deepClone(value);
  notifySubscribers();
}

/**
 * Reset state to defaults
 */
export function resetState() {
  state = deepClone(DEFAULT_STATE);
  notifySubscribers();
}

/**
 * Subscribe to state changes
 * @param {Function} callback - Callback function
 * @returns {Function} Unsubscribe function
 */
export function subscribe(callback) {
  subscribers.add(callback);
  return () => subscribers.delete(callback);
}

/**
 * Notify all subscribers of state change
 */
function notifySubscribers() {
  const currentState = getState();
  subscribers.forEach(callback => {
    try {
      callback(currentState);
    } catch (error) {
      console.error('State subscriber error:', error);
    }
  });
}

/**
 * Deep merge objects
 * @param {Object} target - Target object
 * @param {Object} source - Source object
 * @returns {Object} Merged object
 */
function mergeDeep(target, source) {
  const result = { ...target };

  for (const key in source) {
    if (Object.prototype.hasOwnProperty.call(source, key)) {
      if (
        source[key] &&
        typeof source[key] === 'object' &&
        !Array.isArray(source[key])
      ) {
        result[key] = mergeDeep(result[key] || {}, source[key]);
      } else {
        result[key] = deepClone(source[key]);
      }
    }
  }

  return result;
}

/**
 * Export state to JSON
 * @returns {string} JSON string
 */
export function exportState() {
  return JSON.stringify({
    version: '1.0.0',
    exportDate: new Date().toISOString(),
    data: state.data,
    config: state.config,
    economic: state.economic,
    tsa: state.tsa
  }, null, 2);
}

/**
 * Import state from JSON
 * @param {string} json - JSON string
 * @returns {boolean} Success
 */
export function importState(json) {
  try {
    const imported = JSON.parse(json);

    if (!imported.data || !Array.isArray(imported.data)) {
      throw new Error('Invalid data format');
    }

    updateState({
      data: imported.data,
      config: imported.config || state.config,
      economic: imported.economic || state.economic,
      tsa: imported.tsa || state.tsa
    });

    return true;
  } catch (error) {
    console.error('Import error:', error);
    return false;
  }
}

/**
 * Import data from CSV
 * @param {string} csv - CSV string
 * @returns {boolean} Success
 */
export function importCSV(csv) {
  try {
    const lines = csv.trim().split('\n');
    if (lines.length < 2) return false;

    const headers = lines[0].split(',').map(h => h.trim().toLowerCase());
    const data = [];

    for (let i = 1; i < lines.length; i++) {
      const values = parseCSVLine(lines[i]);
      if (values.length !== headers.length) continue;

      const row = {};
      headers.forEach((header, j) => {
        row[header] = values[j];
      });
      data.push(row);
    }

    updateState({ data });
    return true;
  } catch (error) {
    console.error('CSV import error:', error);
    return false;
  }
}

/**
 * Parse a CSV line handling quoted values
 * @param {string} line - CSV line
 * @returns {Array} Parsed values
 */
function parseCSVLine(line) {
  const values = [];
  let current = '';
  let inQuotes = false;

  for (let i = 0; i < line.length; i++) {
    const char = line[i];

    if (char === '"') {
      inQuotes = !inQuotes;
    } else if (char === ',' && !inQuotes) {
      values.push(current.trim());
      current = '';
    } else {
      current += char;
    }
  }

  values.push(current.trim());
  return values;
}

/**
 * Export data to CSV
 * @returns {string} CSV string
 */
export function exportCSV() {
  const data = state.data;
  if (!data || data.length === 0) return '';

  const headers = Object.keys(data[0]);
  const lines = [headers.join(',')];

  for (const row of data) {
    const values = headers.map(h => {
      const val = row[h];
      if (val === null || val === undefined) return '';
      if (typeof val === 'string' && (val.includes(',') || val.includes('"'))) {
        return `"${val.replace(/"/g, '""')}"`;
      }
      return String(val);
    });
    lines.push(values.join(','));
  }

  return lines.join('\n');
}

/**
 * Add a study row
 * @param {Object} row - Row data (optional)
 */
export function addStudyRow(row = {}) {
  const newData = [...state.data, row];
  updateState({ data: newData });
}

/**
 * Update a study row
 * @param {number} index - Row index
 * @param {string} column - Column name
 * @param {*} value - New value
 */
export function updateStudyRow(index, column, value) {
  if (index < 0 || index >= state.data.length) return;

  const newData = [...state.data];
  newData[index] = { ...newData[index], [column]: value };
  updateState({ data: newData });
}

/**
 * Delete a study row
 * @param {number} index - Row index
 */
export function deleteStudyRow(index) {
  if (index < 0 || index >= state.data.length) return;

  const newData = state.data.filter((_, i) => i !== index);
  updateState({ data: newData });
}

/**
 * Clear all study data
 */
export function clearStudyData() {
  updateState({ data: [], currentResult: null });
}
