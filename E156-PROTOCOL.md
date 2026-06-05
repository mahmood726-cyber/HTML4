# E156 Protocol — HTA Meta-Analysis Suite

- **Project**: HTML4 — HTA Meta-Analysis Suite
- **Revived**: 2026-06-05
- **Type**: Single-file offline browser tool (meta-analysis), with an extracted pure JS engine + Node tests
- **Dashboard**: https://mahmood726-cyber.github.io/HTML4/

## What changed

- Vendored Plotly 2.27.0 locally (`plotly.min.js`) and removed the Google Fonts
  link; zero external URLs remain in `index.html`.
- Extracted the pure statistical engine verbatim into `engine.js` (single source
  of truth) and deleted the inline duplicate.
- Added `tests.js` (45 pure-Node checks) with independently hand-derived
  expected values, including a fully worked 2-study OR pooling example, a
  3-study DerSimonian-Laird example with heterogeneity, and edge cases.
- Renamed the extensionless `HTML4` file to `index.html`; added `.nojekyll`,
  `.gitignore`, and `README.md`.
- No statistical methodology was changed.

## Body (E156 draft — CURRENT BODY)

How reproducibly can a browser-only tool reproduce standard random-effects meta-analysis results across the effect measures an HTA reviewer actually uses? The instrument is a single offline HTML page whose pure statistical engine was extracted verbatim into a separate module so that browser and test paths share one implementation. We pool on the log scale for ratio measures and apply inverse-variance random effects with DerSimonian-Laird tau-squared, an HKSJ Knapp-Hartung floor, and a t-based prediction interval. A worked two-study odds-ratio example reproduces a pooled OR of 0.491 (95% CI 0.287 to 0.839) and a three-study example reproduces tau-squared 0.0525 and I-squared 77.8 percent, matched to hand derivations. Forty-five unit tests pass, covering k equal to one passthrough, identical-study zero-heterogeneity, and empty-input guards. The engine reports heterogeneity on random-effects weights, which equals the Knapp-Hartung statistic but differs from fixed-effect Cochran's Q. This is a teaching and triage tool, not a regulatory replacement for metafor.

SUBMITTED: [ ]
