// Pure-Node tests for engine.js (HTA Meta-Analysis Suite).
// Every expected value below is hand-derived independently (by algebra / a
// separate calculator), NOT by running the engine. Run: node tests.js
'use strict';
const E = require('./engine.js');

let pass = 0, fail = 0;
function approx(label, got, exp, tol) {
    tol = (tol === undefined) ? 1e-6 : tol;
    const ok = (got === null && exp === null) ||
               (typeof got === 'number' && typeof exp === 'number' && Math.abs(got - exp) <= tol);
    if (ok) { pass++; console.log('  PASS  ' + label + '  (got ' + got + ')'); }
    else    { fail++; console.log('  FAIL  ' + label + '  expected ' + exp + ' got ' + got); }
}
function eq(label, got, exp) {
    const ok = got === exp;
    if (ok) { pass++; console.log('  PASS  ' + label + '  (got ' + JSON.stringify(got) + ')'); }
    else    { fail++; console.log('  FAIL  ' + label + '  expected ' + JSON.stringify(exp) + ' got ' + JSON.stringify(got)); }
}

console.log('--- Distribution helpers ---');
// normCDF(0) must be 0.5 (a value of 0 here would be a bug).
approx('normCDF(0) = 0.5', E.normCDF(0), 0.5, 1e-7);
// Phi(1.96) ~= 0.9750021; Phi(-1.96) ~= 0.0249979.
approx('normCDF(1.96) ~= 0.975', E.normCDF(1.96), 0.9750021, 1e-4);
approx('normCDF(-1.96) ~= 0.025', E.normCDF(-1.96), 0.0249979, 1e-4);
// Two-sided p for z=1.96 ~= 0.05.
approx('pFromZ(1.96) ~= 0.05', E.pFromZ(1.96), 0.04999579, 2e-4);
// Student-t critical values (table): t_{1,.975}=12.706, t_{2,.975}=4.303, t_{5,.975}=2.571.
approx('tCrit(df=1) = 12.706', E.tCrit(1, 0.05), 12.706, 1e-3);
approx('tCrit(df=2) = 4.303',  E.tCrit(2, 0.05), 4.303, 1e-3);
approx('tCrit(df=5) = 2.571',  E.tCrit(5, 0.05), 2.571, 1e-3);

console.log('\n--- Effect-size calc: binary OR (Study A: 10/100 vs 20/100) ---');
// logOR = ln((10*80)/(20*90)) = ln(800/1800) = -0.810930216
// vi    = 1/10 + 1/90 + 1/20 + 1/80 = 0.173611111
const effA = E.calcEffect({e1:10,n1:100,e2:20,n2:100}, 'binary', 'OR');
approx('logOR Study A', effA.es, -0.810930216, 1e-7);
approx('vi Study A',    effA.vi, 0.173611111, 1e-7);
approx('display OR Study A', effA.display, Math.exp(-0.810930216), 1e-7);

console.log('\n--- Effect-size calc: correlation (Fisher z), r=0.5, n=28 ---');
// z = 0.5*ln(1.5/0.5) = 0.5*ln3 = 0.549306144 ; vi = 1/(n-3) = 1/25 = 0.04
const effR = E.calcEffect({r:0.5,n:28}, 'correlation', null);
approx('Fisher z (r=0.5)', effR.es, 0.549306144, 1e-7);
approx('Fisher z var = 1/(n-3)', effR.vi, 0.04, 1e-9);

console.log('\n--- HAND-WORKED 2-study OR pooling (tau2 truncates to 0) ---');
// Study A: 10/100 vs 20/100 -> esA=-0.810930216, viA=0.173611111
// Study B: 15/100 vs 25/100 -> esB=ln((15*75)/(25*85))=ln(1125/2125)=-0.635988774
//          viB = 1/15 + 1/85 + 1/25 + 1/75 = 0.131764706
// wA=1/viA=5.76, wB=1/viB=7.589285714, sumW=13.349285714
// muFE = (wA*esA + wB*esB)/sumW = -0.711473157
// Q = wA*(esA-mu)^2 + wB*(esB-mu)^2 = 0.100219174  (< df=1) -> tau2=max(0,(Q-1)/C)=0, I2=0
// With tau2=0, RE == FE: muRE=-0.711473157, seRE=sqrt(1/sumW)=0.273697600
// OR = exp(muRE)=0.490920 ; z=muRE/seRE=-2.599486 ; CI(95%,z) OR=[0.287101,0.839435]
const A = E.calcEffect({e1:10,n1:100,e2:20,n2:100}, 'binary', 'OR');
const B = E.calcEffect({e1:15,n1:100,e2:25,n2:100}, 'binary', 'OR');
const studs = [{...A, excluded:false}, {...B, excluded:false}];
const tau2_2 = E.tauDL(studs.map(s => ({...s})));
approx('2-study tauDL = 0 (Q<df)', tau2_2, 0, 1e-12);
const pool2 = E.poolIV(studs.map(s => ({...s})), tau2_2, false, 0.95);
approx('2-study pooled logOR', pool2.es, -0.711473157, 1e-6);
approx('2-study pooled SE',    pool2.se, 0.273697600, 1e-6);
approx('2-study Q',            pool2.Q,  0.100219174, 1e-6);
approx('2-study I2 = 0',       pool2.I2, 0, 1e-9);
approx('2-study pooled OR',    Math.exp(pool2.es), 0.490920, 1e-5);
approx('2-study z',            pool2.z, -2.599486, 1e-5);
approx('2-study CI low (OR)',  Math.exp(pool2.ciLo), 0.287101, 1e-5);
approx('2-study CI high (OR)', Math.exp(pool2.ciHi), 0.839435, 1e-5);

console.log('\n--- HAND-WORKED 3-study generic pooling (DL with heterogeneity) ---');
// S1: es=0.20, se=0.10 (vi=0.01) ; S2: es=0.50, se=0.10 (vi=0.01) ; S3: es=0.80, se=0.20 (vi=0.04)
//
// tauDL uses FE (inverse-variance) weights internally:
//   w=[100,100,25], sumW=225, sumW2=20625
//   muFE = (20+50+20)/225 = 0.40 ; Q_FE = 100*.04+100*.01+25*.16 = 9
//   C = 225 - 20625/225 = 133.3333333 ; tau2 = (9-2)/133.3333333 = 0.0525  (verified below)
//
// poolIV then re-weights with RE weights w*=1/(vi+tau2):
//   w* = [16, 16, 10.810811], sumWr = 42.810811
//   muRE = (16*.2+16*.5+10.810811*.8)/42.810811 = 0.463636364
//   seRE = sqrt(1/42.810811) = 0.152835161
//
// NOTE on Q/I2/H2 reported by poolIV: poolIV recomputes these with the *RE*
// weights about muRE (this equals (k-1) * the Knapp-Hartung q statistic it uses
// for HKSJ below), NOT the FE-weighted Cochran's Q. Hand-derived RE-weighted:
//   Q_RE = 16*(0.2-0.463636)^2 + 16*(0.5-0.463636)^2 + 10.810811*(0.8-0.463636)^2
//        = 2.35636364
//   I2_RE = (2.35636364-2)/2.35636364*100 = 15.1234568% ; H2_RE = 2.35636364/2 = 1.17818182
const g = [
    {es:0.20, vi:0.01, se:0.10, excluded:false},
    {es:0.50, vi:0.01, se:0.10, excluded:false},
    {es:0.80, vi:0.04, se:0.20, excluded:false}
];
const tau2_3 = E.tauDL(g.map(s => ({...s})));
approx('3-study tauDL = 0.0525 (FE-weighted Q=9 internally)', tau2_3, 0.0525, 1e-9);
const pool3 = E.poolIV(g.map(s => ({...s})), tau2_3, false, 0.95);
approx('3-study pooled ES',  pool3.es,  0.463636364, 1e-7);
approx('3-study pooled SE',  pool3.se,  0.152835161, 1e-7);
approx('3-study Q (RE-weighted) = 2.356364', pool3.Q, 2.356363636, 1e-7);
approx('3-study I2 (RE-weighted) = 15.1235%', pool3.I2, 15.123456790, 1e-6);
approx('3-study H2 (RE-weighted) = 1.178182', pool3.H2, 1.178181818, 1e-7);

console.log('\n--- HKSJ (Knapp-Hartung) adjustment on 3-study example ---');
// Knapp-Hartung q = (1/df) * Q_RE = 2.35636364/2 = 1.17818182  (this is qStar in code)
// floor = max(1, 1.17818182) = 1.17818182 ; seHKSJ = seRE*sqrt(1.17818182) = 0.165893524
// tCrit(df=2)=4.303 ; CI = muRE +/- 4.303*seHKSJ = [-0.25020347, 1.17747620]
const pool3h = E.poolIV(g.map(s => ({...s})), tau2_3, true, 0.95);
approx('3-study HKSJ SE',     pool3h.se,   0.165893524, 1e-7);
approx('3-study HKSJ CI low', pool3h.ciLo, -0.25020347, 1e-5);
approx('3-study HKSJ CI high',pool3h.ciHi,  1.17747620, 1e-5);

console.log('\n--- Prediction interval (no-HKSJ branch) ---');
// piSe = sqrt(seRE^2 + tau2) = sqrt(0.152835161^2 + 0.0525) = sqrt(0.0758586)=0.275424
// piCrit = tCrit(df=2)=4.303 ; PI = 0.463636 +/- 4.303*0.275424 = [-0.721515, 1.648787]
approx('3-study PI low',  pool3.piLo, -0.721515, 1e-4);
approx('3-study PI high', pool3.piHi,  1.648787, 1e-4);

console.log('\n--- Edge case: k=1 passthrough ---');
// Single study es=0.3, vi=0.04 -> pooled es=0.3, se=sqrt(0.04)=0.2, Q=0, df=0, I2=0
const e1 = [{es:0.3, vi:0.04, se:0.2, excluded:false}];
const r1 = E.poolIV(e1.map(s=>({...s})), 0, false, 0.95);
approx('k=1 pooled es', r1.es, 0.3, 1e-12);
approx('k=1 pooled se', r1.se, 0.2, 1e-12);
approx('k=1 Q = 0',     r1.Q, 0, 1e-12);
eq('k=1 df = 0', r1.df, 0);
approx('k=1 I2 = 0', r1.I2, 0, 1e-12);
eq('k=1 tauDL = 0 (k<2 guard)', E.tauDL(e1), 0);

console.log('\n--- Edge case: two identical studies => tau2=0, I2=0 ---');
const ident = [
    {es:0.5, vi:0.02, se:Math.sqrt(0.02), excluded:false},
    {es:0.5, vi:0.02, se:Math.sqrt(0.02), excluded:false}
];
approx('identical tauDL = 0', E.tauDL(ident.map(s=>({...s}))), 0, 1e-12);
const rI = E.poolIV(ident.map(s=>({...s})), 0, false, 0.95);
approx('identical Q = 0',  rI.Q, 0, 1e-12);
approx('identical I2 = 0', rI.I2, 0, 1e-12);
approx('identical pooled es = 0.5', rI.es, 0.5, 1e-12);

console.log('\n--- Edge case: empty / all-excluded guard ---');
eq('poolIV([]) returns null', E.poolIV([], 0, false, 0.95), null);
const allExcl = [{es:0.3, vi:0.04, se:0.2, excluded:true}];
eq('poolIV(all-excluded) returns null', E.poolIV(allExcl, 0, false, 0.95), null);

console.log('\n--- Egger test on 3-study generic example ---');
// x=1/se=[10,10,5], y=es/se=[2,5,4]
// sumX=25, sumY=11, sumXY=10*2+10*5+5*4=90, sumX2=100+100+25=225
// slope=(n*sumXY-sumX*sumY)/(n*sumX2-sumX^2)=(3*90-25*11)/(3*225-625)=(270-275)/(675-625)=-5/50=-0.1
// intercept=(sumY-slope*sumX)/n=(11-(-0.1*25))/3=(11+2.5)/3=13.5/3=4.5
const egg = E.eggerTest(g.map(s=>({...s})));
approx('Egger slope check via intercept', egg.intercept, 4.5, 1e-9);

console.log('\n' + pass + ' passed, ' + fail + ' failed');
process.exit(fail === 0 ? 0 : 1);
