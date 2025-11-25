/**
 * Utility functions for the HTA Meta-Analysis Suite
 */

/**
 * HTML entity map for escaping
 */
const HTML_ENTITIES = {
  '&': '&amp;',
  '<': '&lt;',
  '>': '&gt;',
  '"': '&quot;',
  "'": '&#39;'
};

/**
 * Escape HTML entities to prevent XSS attacks
 * @param {string} str - String to escape
 * @returns {string} Escaped string
 */
export function escapeHTML(str) {
  if (str === null || str === undefined) return '';
  return String(str).replace(/[&<>"']/g, char => HTML_ENTITIES[char]);
}

/**
 * Parse a value to number, returning null if invalid
 * @param {*} value - Value to parse
 * @returns {number|null} Parsed number or null
 */
export function parseNumber(value) {
  if (value === null || value === undefined || value === '') return null;
  const num = parseFloat(value);
  return Number.isFinite(num) ? num : null;
}

/**
 * Format a number with specified decimal places
 * @param {number} value - Number to format
 * @param {number} decimals - Decimal places (default: 2)
 * @returns {string} Formatted string or '—' if invalid
 */
export function formatNumber(value, decimals = 2) {
  if (value === null || value === undefined || !Number.isFinite(value)) {
    return '—';
  }
  return value.toFixed(decimals);
}

/**
 * Format a p-value with appropriate precision
 * @param {number} p - P-value
 * @returns {string} Formatted p-value
 */
export function formatPValue(p) {
  if (p === null || p === undefined || !Number.isFinite(p)) {
    return '—';
  }
  if (p < 0.001) return '<0.001';
  return p.toFixed(3);
}

/**
 * Format effect size (handles ratio vs difference metrics)
 * @param {number} value - Effect size on log scale for ratios
 * @param {boolean} isRatio - Whether the metric is ratio-based
 * @param {number} decimals - Decimal places
 * @returns {string} Formatted effect size
 */
export function formatEffectSize(value, isRatio = false, decimals = 2) {
  if (value === null || value === undefined || !Number.isFinite(value)) {
    return '—';
  }
  const displayValue = isRatio ? Math.exp(value) : value;
  return displayValue.toFixed(isRatio ? decimals : 3);
}

/**
 * Clamp a value between min and max
 * @param {number} value - Value to clamp
 * @param {number} min - Minimum value
 * @param {number} max - Maximum value
 * @returns {number} Clamped value
 */
export function clamp(value, min, max) {
  return Math.max(min, Math.min(max, value));
}

/**
 * Check if a value is a valid positive number
 * @param {*} value - Value to check
 * @returns {boolean} True if valid positive number
 */
export function isPositiveNumber(value) {
  const num = parseNumber(value);
  return num !== null && num > 0;
}

/**
 * Check if a value is a valid probability (0-1)
 * @param {*} value - Value to check
 * @returns {boolean} True if valid probability
 */
export function isProbability(value) {
  const num = parseNumber(value);
  return num !== null && num >= 0 && num <= 1;
}

/**
 * Debounce function execution
 * @param {Function} func - Function to debounce
 * @param {number} wait - Wait time in ms
 * @returns {Function} Debounced function
 */
export function debounce(func, wait = 250) {
  let timeout;
  return function executedFunction(...args) {
    const later = () => {
      clearTimeout(timeout);
      func(...args);
    };
    clearTimeout(timeout);
    timeout = setTimeout(later, wait);
  };
}

/**
 * Deep clone an object
 * @param {*} obj - Object to clone
 * @returns {*} Cloned object
 */
export function deepClone(obj) {
  if (obj === null || typeof obj !== 'object') return obj;
  if (Array.isArray(obj)) return obj.map(deepClone);
  return Object.fromEntries(
    Object.entries(obj).map(([key, val]) => [key, deepClone(val)])
  );
}

/**
 * Generate a unique ID
 * @param {string} prefix - Optional prefix
 * @returns {string} Unique ID
 */
export function generateId(prefix = 'id') {
  return `${prefix}_${Date.now()}_${Math.random().toString(36).slice(2, 9)}`;
}

/**
 * Validate study data row based on data type
 * @param {Object} row - Data row
 * @param {string} dataType - Type of data (binary, continuous, etc.)
 * @returns {Object} Validation result { valid: boolean, errors: string[] }
 */
export function validateStudyRow(row, dataType) {
  const errors = [];

  if (!row.study || row.study.trim() === '') {
    errors.push('Study name is required');
  }

  switch (dataType) {
    case 'binary':
      if (parseNumber(row.n1) !== null && parseNumber(row.n1) <= 0) {
        errors.push('N1 must be positive');
      }
      if (parseNumber(row.n2) !== null && parseNumber(row.n2) <= 0) {
        errors.push('N2 must be positive');
      }
      if (parseNumber(row.e1) !== null && parseNumber(row.e1) < 0) {
        errors.push('E1 cannot be negative');
      }
      if (parseNumber(row.e2) !== null && parseNumber(row.e2) < 0) {
        errors.push('E2 cannot be negative');
      }
      if (parseNumber(row.e1) > parseNumber(row.n1)) {
        errors.push('E1 cannot exceed N1');
      }
      if (parseNumber(row.e2) > parseNumber(row.n2)) {
        errors.push('E2 cannot exceed N2');
      }
      break;

    case 'continuous':
      if (parseNumber(row.s1) !== null && parseNumber(row.s1) <= 0) {
        errors.push('SD1 must be positive');
      }
      if (parseNumber(row.s2) !== null && parseNumber(row.s2) <= 0) {
        errors.push('SD2 must be positive');
      }
      if (parseNumber(row.n1) !== null && parseNumber(row.n1) <= 0) {
        errors.push('N1 must be positive');
      }
      if (parseNumber(row.n2) !== null && parseNumber(row.n2) <= 0) {
        errors.push('N2 must be positive');
      }
      break;

    case 'survival':
      if (parseNumber(row.hr) !== null && parseNumber(row.hr) <= 0) {
        errors.push('HR must be positive');
      }
      if (parseNumber(row.ll) !== null && parseNumber(row.ul) !== null) {
        if (parseNumber(row.ll) >= parseNumber(row.ul)) {
          errors.push('CI lower bound must be less than upper bound');
        }
      }
      break;

    case 'correlation':
      const r = parseNumber(row.r);
      if (r !== null && (r < -1 || r > 1)) {
        errors.push('Correlation must be between -1 and 1');
      }
      if (parseNumber(row.n) !== null && parseNumber(row.n) <= 3) {
        errors.push('N must be greater than 3 for correlation');
      }
      break;

    case 'generic':
      if (parseNumber(row.se) !== null && parseNumber(row.se) <= 0) {
        errors.push('SE must be positive');
      }
      break;
  }

  return {
    valid: errors.length === 0,
    errors
  };
}

/**
 * Validate analysis configuration
 * @param {Object} config - Configuration object
 * @returns {Object} Validation result { valid: boolean, errors: string[] }
 */
export function validateConfig(config) {
  const errors = [];

  if (config.confLevel !== undefined) {
    const cl = parseNumber(config.confLevel);
    if (cl === null || cl <= 0 || cl >= 1) {
      errors.push('Confidence level must be between 0 and 1');
    }
  }

  if (config.cer !== undefined) {
    const cer = parseNumber(config.cer);
    if (cer === null || cer <= 0 || cer >= 1) {
      errors.push('Control event rate must be between 0 and 1');
    }
  }

  if (config.tsaAlpha !== undefined) {
    const alpha = parseNumber(config.tsaAlpha);
    if (alpha === null || alpha <= 0 || alpha >= 0.5) {
      errors.push('TSA alpha must be between 0 and 0.5');
    }
  }

  if (config.tsaBeta !== undefined) {
    const beta = parseNumber(config.tsaBeta);
    if (beta === null || beta <= 0 || beta >= 1) {
      errors.push('TSA beta must be between 0 and 1');
    }
  }

  return {
    valid: errors.length === 0,
    errors
  };
}
