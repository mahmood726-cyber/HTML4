/**
 * GRADE assessment renderer
 *
 * FIXED: Imprecision logic now correctly checks null crossing for ratio metrics
 * (null is 1 on exp scale, not 0 on log scale)
 */

import { getElement, setHTML } from '../app/dom.js';
import { escapeHTML, formatEffectSize, formatNumber } from '../utils/index.js';
import { GRADE_LEVELS, I2_THRESHOLDS, ROB_THRESHOLDS } from '../constants/index.js';

/**
 * Render GRADE assessment
 * @param {Object} result - Pooled result
 */
export function renderGRADE(result) {
  const assessment = calculateGRADE(result);

  renderGradeProfile(result, assessment);
  renderSoFTable(result, assessment);
}

/**
 * Calculate GRADE certainty assessment
 * @param {Object} result - Pooled result
 * @returns {Object} GRADE assessment
 */
function calculateGRADE(result) {
  // Start with High for RCTs (assumption - could be parameterized)
  let certaintyLevel = 4; // 4=High, 3=Moderate, 2=Low, 1=Very Low
  const downgrades = [];
  const upgrades = [];

  const active = result.studies.filter(s => !s.excluded);

  // 1. Risk of Bias assessment
  const robCounts = active.reduce((acc, s) => {
    const rob = s.rob || 'unclear';
    acc[rob] = (acc[rob] || 0) + 1;
    return acc;
  }, {});

  const highRiskPct = (robCounts.high || 0) / result.k;
  if (highRiskPct > ROB_THRESHOLDS.highRiskRatio) {
    downgrades.push({
      domain: 'Risk of Bias',
      severity: highRiskPct > 0.5 ? 'Very serious' : 'Serious',
      reason: `${((robCounts.high || 0) / result.k * 100).toFixed(0)}% of studies at high risk`
    });
    certaintyLevel -= highRiskPct > 0.5 ? 2 : 1;
  }

  // 2. Inconsistency (heterogeneity)
  if (result.I2 > I2_THRESHOLDS.high) {
    downgrades.push({
      domain: 'Inconsistency',
      severity: 'Very serious',
      reason: `I² = ${formatNumber(result.I2, 1)}% (>75%)`
    });
    certaintyLevel -= 2;
  } else if (result.I2 > I2_THRESHOLDS.moderate) {
    downgrades.push({
      domain: 'Inconsistency',
      severity: 'Serious',
      reason: `I² = ${formatNumber(result.I2, 1)}% (>50%)`
    });
    certaintyLevel -= 1;
  }

  // 3. Indirectness (would need external info - placeholder)
  // This would typically be assessed based on PICO matching

  // 4. Imprecision
  // FIXED: Correct null crossing check for ratio metrics
  const imprecisionResult = assessImprecision(result);
  if (imprecisionResult.downgrade) {
    downgrades.push({
      domain: 'Imprecision',
      severity: imprecisionResult.severity,
      reason: imprecisionResult.reason
    });
    certaintyLevel -= imprecisionResult.severity === 'Very serious' ? 2 : 1;
  }

  // 5. Publication Bias
  // Would use Egger test result - simplified here
  if (result.k < 10) {
    // Can't reliably assess publication bias with few studies
  }

  // Clamp certainty level
  certaintyLevel = Math.max(1, Math.min(4, certaintyLevel));

  const certaintyLabels = {
    4: GRADE_LEVELS.HIGH,
    3: GRADE_LEVELS.MODERATE,
    2: GRADE_LEVELS.LOW,
    1: GRADE_LEVELS.VERY_LOW
  };

  return {
    certainty: certaintyLabels[certaintyLevel],
    certaintyLevel,
    downgrades,
    upgrades,
    robCounts
  };
}

/**
 * Assess imprecision for GRADE
 * FIXED: Now correctly handles ratio metrics
 * @param {Object} result - Pooled result
 * @returns {Object} Imprecision assessment
 */
function assessImprecision(result) {
  const { isRatio, es, ciLo, ciHi, totalN, k } = result;

  // Check if CI crosses the null effect
  let crossesNull;
  if (isRatio) {
    // For ratio metrics (OR, RR, HR), null is 1 on the natural scale
    // CI is on log scale, so we need to exponentiate to check
    const expLo = Math.exp(ciLo);
    const expHi = Math.exp(ciHi);
    crossesNull = expLo < 1 && expHi > 1;
  } else {
    // For difference metrics (MD, SMD, RD), null is 0
    crossesNull = ciLo < 0 && ciHi > 0;
  }

  // Check CI width (rough heuristic)
  let wideCI = false;
  if (isRatio) {
    const expLo = Math.exp(ciLo);
    const expHi = Math.exp(ciHi);
    // Wide if CI ratio > 3 (e.g., 0.5 to 1.5 vs 0.3 to 3.0)
    wideCI = (expHi / expLo) > 3;
  } else {
    // Would need clinical context to assess width
    wideCI = (ciHi - ciLo) > Math.abs(es) * 2;
  }

  // Check sample size (Optimal Information Size)
  const smallSample = totalN < 300 || (totalN < 400 && k < 5);

  // Determine downgrade
  const problems = [];
  if (crossesNull) problems.push('CI crosses null');
  if (wideCI) problems.push('wide confidence interval');
  if (smallSample) problems.push(`small sample size (N=${totalN})`);

  if (problems.length === 0) {
    return { downgrade: false };
  }

  const severity = problems.length >= 2 ? 'Very serious' : 'Serious';
  return {
    downgrade: true,
    severity,
    reason: problems.join(', ')
  };
}

/**
 * Render GRADE profile table
 * @param {Object} result - Pooled result
 * @param {Object} assessment - GRADE assessment
 */
function renderGradeProfile(result, assessment) {
  const gradeClass = {
    [GRADE_LEVELS.HIGH]: 'grade-high',
    [GRADE_LEVELS.MODERATE]: 'grade-moderate',
    [GRADE_LEVELS.LOW]: 'grade-low',
    [GRADE_LEVELS.VERY_LOW]: 'grade-very-low'
  }[assessment.certainty];

  const downgradeCells = ['Risk of Bias', 'Inconsistency', 'Indirectness', 'Imprecision', 'Publication Bias']
    .map(domain => {
      const downgrade = assessment.downgrades.find(d => d.domain === domain);
      if (downgrade) {
        return `<td title="${escapeHTML(downgrade.reason)}">
          ${downgrade.severity === 'Very serious' ? '-2' : '-1'}
        </td>`;
      }
      return '<td>0</td>';
    })
    .join('');

  const html = `
    <table class="grade-table">
      <thead>
        <tr>
          <th>Outcome</th>
          <th>Studies (N)</th>
          <th>Risk of Bias</th>
          <th>Inconsistency</th>
          <th>Indirectness</th>
          <th>Imprecision</th>
          <th>Pub. Bias</th>
          <th>Certainty</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <td>Primary Outcome</td>
          <td>${result.k} (${result.totalN || 'N/A'})</td>
          ${downgradeCells}
          <td><span class="grade-badge ${gradeClass}">${assessment.certainty}</span></td>
        </tr>
      </tbody>
    </table>
  `;

  setHTML('gradeProfile', html);
}

/**
 * Render Summary of Findings table
 * @param {Object} result - Pooled result
 * @param {Object} assessment - GRADE assessment
 */
function renderSoFTable(result, assessment) {
  const { isRatio, es, ciLo, ciHi } = result;
  const effectDisplay = formatEffectSize(es, isRatio);
  const ciDisplay = `${formatEffectSize(ciLo, isRatio)} to ${formatEffectSize(ciHi, isRatio)}`;

  const gradeClass = {
    [GRADE_LEVELS.HIGH]: 'grade-high',
    [GRADE_LEVELS.MODERATE]: 'grade-moderate',
    [GRADE_LEVELS.LOW]: 'grade-low',
    [GRADE_LEVELS.VERY_LOW]: 'grade-very-low'
  }[assessment.certainty];

  // Generate downgrade explanations
  const explanations = assessment.downgrades.length > 0
    ? assessment.downgrades.map(d => `• ${d.domain}: ${d.reason}`).join('<br>')
    : 'No serious concerns';

  const html = `
    <table class="grade-table">
      <thead>
        <tr>
          <th>Outcome</th>
          <th>Effect (95% CI)</th>
          <th>No. of Studies</th>
          <th>Certainty</th>
          <th>Interpretation</th>
        </tr>
      </thead>
      <tbody>
        <tr>
          <td>Primary Outcome</td>
          <td><strong>${effectDisplay}</strong><br><small>${ciDisplay}</small></td>
          <td>${result.k} RCTs<br><small>N = ${result.totalN || 'N/A'}</small></td>
          <td><span class="grade-badge ${gradeClass}">${assessment.certainty}</span></td>
          <td style="font-size: 0.75rem;">${explanations}</td>
        </tr>
      </tbody>
    </table>
  `;

  setHTML('sofTable', html);
}

/**
 * Render Risk of Bias summary
 * @param {Object} result - Pooled result
 */
export function renderRoB(result) {
  const active = result.studies.filter(s => !s.excluded);

  // Count RoB levels
  const counts = { low: 0, some: 0, high: 0, unclear: 0 };
  active.forEach(s => {
    const rob = s.rob || 'unclear';
    counts[rob]++;
  });

  // Traffic light table
  const rows = active.map(s => {
    const rob = s.rob || 'unclear';
    const robClass = {
      low: 'rob-low',
      some: 'rob-some',
      high: 'rob-high',
      unclear: 'rob-unclear'
    }[rob];

    return `
      <tr>
        <td>${escapeHTML(s.id)}</td>
        <td><span class="rob-cell ${robClass}" title="${rob}"></span></td>
      </tr>
    `;
  }).join('');

  const html = `
    <div style="margin-bottom: 1rem;">
      <strong>Summary:</strong>
      <span style="color: #22c55e;">● Low: ${counts.low}</span>
      <span style="color: #eab308; margin-left: 1rem;">● Some: ${counts.some}</span>
      <span style="color: #ef4444; margin-left: 1rem;">● High: ${counts.high}</span>
      <span style="color: #94a3b8; margin-left: 1rem;">● Unclear: ${counts.unclear}</span>
    </div>
    <table class="rob-table">
      <thead>
        <tr>
          <th>Study</th>
          <th>Overall RoB</th>
        </tr>
      </thead>
      <tbody>
        ${rows}
      </tbody>
    </table>
  `;

  setHTML('robSummary', html);
}
