// Client-side chemistry evaluator module
// Loads CSVs, balances reactions, searches reactions/compounds, and renders into the DOM.

import StringSearcher from './searcher.js';

async function fetchText(path) {
  const res = await fetch(path);
  if (!res.ok) throw new Error(`Failed to load ${path}: ${res.status}`);
  return res.text();
}

function uniqueReactionTokens(db) {
  const rSet = new Set();
  const pSet = new Set();
  db.reactions.forEach(r => {
    (r.reactants || []).forEach(x => rSet.add(x));
    (r.products || []).forEach(x => pSet.add(x));
  });
  return { reactants: Array.from(rSet).sort(), products: Array.from(pSet).sort() };
}

function parseCSV(text) {
  const raw = text.trim().split(/\r?\n/);
  if (!raw.length) return [];
  const headers = raw[0].split(',');
  const lines = raw.slice(1).filter(l => {
    const t = l.trim();
    return t && !t.startsWith('#') && !t.startsWith('//');
  });
  return lines.map(line => {
    const vals = [];
    let cur = '', inQ = false;
    for (const ch of line) {
      if (ch === '"') inQ = !inQ;
      else if (ch === ',' && !inQ) { vals.push(cur); cur = ''; }
      else cur += ch;
    }
    vals.push(cur);
    const obj = {};
    headers.forEach((h, i) => {
      const vRaw = vals[i] ?? '';
      const v = vRaw.trim();
      const num = Number(v);
      obj[h.trim()] = v !== '' && !Number.isNaN(num) ? num : v;
    });
    return obj;
  });
}

async function loadDb() {
  const [elementsTxt, compoundsTxt, reactionsTxt] = await Promise.all([
    fetchText('elements.csv'),
    fetchText('compounds.csv'),
    fetchText('reactions.csv'),
  ]);
  const elementsArr = parseCSV(elementsTxt);
  const elements = elementsArr.reduce((acc, row) => { acc[row.symbol] = row; return acc; }, {});
  const compounds = parseCSV(compoundsTxt);
  const reactionsRaw = parseCSV(reactionsTxt);
  const compoundNameMap = compounds.reduce((acc, c) => {
    const formula = String(c.formula || '').trim();
    const name = String(c.name || '').trim();
    if (formula) acc.formulaToNames[formula.toLowerCase()] = acc.formulaToNames[formula.toLowerCase()] ? [...acc.formulaToNames[formula.toLowerCase()], name] : [name];
    if (name) acc.nameToFormula[name.toLowerCase()] = formula;
    return acc;
  }, { formulaToNames: {}, nameToFormula: {} });

  const reactions = reactionsRaw.map(r => {
    const { reactants, products } = parseReactionSides(String(r.equation || ''));
    return { ...r, reactants, products };
  });

  const tokens = uniqueReactionTokens({ reactions });
  console.log('[chemistry] reactant tokens:', tokens.reactants);
  console.log('[chemistry] product tokens:', tokens.products);

  return { elements, compounds, reactions, compoundNameMap };
}

function parseFormula(formula) {
  const counts = {};
  const stack = [{}];
  let i = 0;
  while (i < formula.length) {
    const ch = formula[i];
    if (ch === '(') { stack.push({}); i++; continue; }
    if (ch === ')') {
      i++;
      let n = '';
      while (i < formula.length && /\d/.test(formula[i])) n += formula[i++];
      const mult = n === '' ? 1 : parseInt(n);
      const popped = stack.pop();
      const top = stack[stack.length - 1];
      for (const k in popped) top[k] = (top[k] || 0) + popped[k] * mult;
      continue;
    }
    if (/[A-Z]/.test(ch)) {
      let el = ch; i++;
      if (i < formula.length && /[a-z]/.test(formula[i])) el += formula[i++];
      let n = '';
      while (i < formula.length && /\d/.test(formula[i])) n += formula[i++];
      const cnt = n === '' ? 1 : parseInt(n);
      const top = stack[stack.length - 1];
      top[el] = (top[el] || 0) + cnt;
      continue;
    }
    i++;
  }
  return stack[0];
}

function parseTerm(s) {
  const m = s.trim().match(/^(\d+)\s*(.*)$/);
  return m ? { base: parseInt(m[1], 10), formula: m[2].trim() } : { base: 1, formula: s.trim() };
}

function gcd(a, b) { return b === 0 ? a : gcd(b, a % b); }
const lcm = (a, b) => (a / gcd(a, b)) * b;
function makeFrac(n, d = 1) { const g = gcd(Math.abs(n), Math.abs(d)); const s = d < 0 ? -1 : 1; return { n: (n / g) * s, d: Math.abs(d) / g }; }
function fracIsZero(f) { return f.n === 0; }
function fAdd(a, b) { return makeFrac(a.n * b.d + b.n * a.d, a.d * b.d); }
function fSub(a, b) { return makeFrac(a.n * b.d - b.n * a.d, a.d * b.d); }
function fMul(a, b) { return makeFrac(a.n * b.n, a.d * b.d); }
function fDiv(a, b) { if (b.n === 0) throw new Error('Divide by zero'); return makeFrac(a.n * b.d, a.d * b.n); }

function nullspaceIntegerSolution(intMatrix) {
  const rows = intMatrix.length;
  const cols = intMatrix[0]?.length || 0;
  if (cols === 0) throw new Error('Empty matrix');

  const m = intMatrix.map(r => r.map(v => makeFrac(v)));
  let pivotRow = 0;
  const pivotCols = [];

  for (let col = 0; col < cols && pivotRow < rows; col++) {
    let sel = -1;
    for (let r = pivotRow; r < rows; r++) if (!fracIsZero(m[r][col])) { sel = r; break; }
    if (sel === -1) continue;
    [m[pivotRow], m[sel]] = [m[sel], m[pivotRow]];
    const pivotVal = m[pivotRow][col];
    for (let c = col; c < cols; c++) m[pivotRow][c] = fDiv(m[pivotRow][c], pivotVal);
    for (let r = 0; r < rows; r++) {
      if (r === pivotRow) continue;
      const factor = m[r][col];
      if (fracIsZero(factor)) continue;
      for (let c = col; c < cols; c++) m[r][c] = fSub(m[r][c], fMul(factor, m[pivotRow][c]));
    }
    pivotCols.push(col);
    pivotRow++;
  }

  const pivotSet = new Set(pivotCols);
  let freeCol = -1;
  for (let c = cols - 1; c >= 0; c--) { if (!pivotSet.has(c)) { freeCol = c; break; } }
  if (freeCol === -1) throw new Error('No free variable; cannot balance');

  const sol = new Array(cols).fill(makeFrac(0));
  sol[freeCol] = makeFrac(1);
  for (let i = pivotCols.length - 1; i >= 0; i--) {
    const col = pivotCols[i];
    let sum = makeFrac(0);
    for (let c = col + 1; c < cols; c++) {
      if (!fracIsZero(m[i][c])) sum = fAdd(sum, fMul(m[i][c], sol[c]));
    }
    sol[col] = makeFrac(-sum.n, sum.d);
  }

  const lcmDen = sol.reduce((acc, f) => lcm(acc, Math.abs(f.d)), 1);
  const ints = sol.map(f => (f.n * (lcmDen / f.d)));
  const gAll = ints.reduce((acc, v) => gcd(acc, Math.abs(v)), Math.abs(ints[0]) || 1);
  return ints.map(v => v / gAll);
}

function balance(reactionStr, elDb) {
  const parts = reactionStr.split('->');
  if (parts.length !== 2) throw new Error('Invalid format, expected Reactants -> Products');
  const leftTerms = parts[0].split('+').map(parseTerm).filter(t => t.formula);
  const rightTerms = parts[1].split('+').map(parseTerm).filter(t => t.formula);
  const allTerms = [...leftTerms, ...rightTerms];
  const compounds = allTerms.map(t => ({ formula: t.formula, atoms: parseFormula(t.formula), base: t.base }));
  const elements = [...new Set(compounds.flatMap(c => Object.keys(c.atoms)))];

  const matrix = elements.map(el => compounds.map((c, idx) => (idx < leftTerms.length ? 1 : -1) * (c.atoms[el] || 0) * c.base));
  const coeffs = nullspaceIntegerSolution(matrix);

  const withMass = (formula, coeff) => {
    const atoms = parseFormula(formula);
    const molarMass = Object.entries(atoms).reduce((sum, [sym, cnt]) => {
      const el = elDb[sym];
      if (!el) throw new Error(`Unknown element: ${sym}`);
      return sum + (el.molar_mass_g_mol || el.mass || 0) * cnt;
    }, 0);
    return { formula, coeff, atoms, molarMass };
  };

  const reactants = leftTerms.map((t, i) => withMass(compounds[i].formula, coeffs[i] * compounds[i].base));
  const products = rightTerms.map((t, i) => withMass(compounds[i + leftTerms.length].formula, coeffs[i + leftTerms.length] * compounds[i + leftTerms.length].base));

  const g = reactants.concat(products).reduce((acc, x) => gcd(acc, x.coeff), reactants[0]?.coeff || 1);
  reactants.forEach(r => r.coeff = r.coeff / g);
  products.forEach(p => p.coeff = p.coeff / g);

  return { reactants, products };
}

function stripCoeff(s) { return s.trim().replace(/^\d+\s+/, ''); }

function parseReactionSides(eq) {
  const parts = eq.split(/->|--.*-->/);
  if (parts.length < 2) return { reactants: [], products: [] };
  const clean = s => s.split('+').map(stripCoeff).map(t => t.trim()).filter(Boolean).map(t => t.replace(/\s+/g, ''));
  return { reactants: clean(parts[0]), products: clean(parts[1]) };
}

function termAliases(term, db) {
  const lc = term.toLowerCase();
  const aliases = new Set([lc]);

  // If query is a name, add its formula
  const byName = db.compoundNameMap.nameToFormula[lc];
  if (byName) aliases.add(String(byName).toLowerCase());

  // If query is a formula, add all known names
  const namesFromFormula = db.compoundNameMap.formulaToNames[lc];
  if (namesFromFormula) namesFromFormula.forEach(n => n && aliases.add(String(n).toLowerCase()));

  // Direct scan as fallback (covers any missing map entries)
  db.compounds.forEach(c => {
    const f = String(c.formula || '').toLowerCase();
    const n = c.name ? c.name.toLowerCase() : '';
    if (f === lc) {
      aliases.add(f);
      if (n) aliases.add(n);
    }
    if (n === lc) {
      if (f) aliases.add(f);
      aliases.add(n);
    }
  });

  return Array.from(aliases);
}

function matchTermAgainstTargets(aliases, targets, allowFuzzy = true) {
  return aliases.some(a => targets.some(t => {
    const tl = t.toLowerCase();
    if (tl === a) return true;
    return allowFuzzy && StringSearcher.compare(a, tl) >= 0.8;
  }));
}

function searchReactions(db, terms, scope = 'any', mode = 'and', allowFuzzy = true) {
  if (!terms.length) return [];
  const resolved = terms.map(t => termAliases(t, db));
  return db.reactions.filter(r => {
    const targets = scope === 'reactant' ? r.reactants : scope === 'product' ? r.products : [...r.reactants, ...r.products];
    if (mode === 'or') {
      return resolved.some(al => matchTermAgainstTargets(al, targets, allowFuzzy));
    }
    return resolved.every(al => matchTermAgainstTargets(al, targets, allowFuzzy));
  });
}

function searchCompounds(db, terms, allowFuzzy = true) {
  if (!terms.length) return [];
  return db.compounds.filter(c => {
    const f = String(c.formula || '').toLowerCase();
    const n = c.name ? c.name.toLowerCase() : '';
    return terms.every(t => {
      const lc = t.toLowerCase();
      if (!allowFuzzy) return lc === f || lc === n;
      return StringSearcher.compare(lc, f) >= 0.6 || (n && StringSearcher.compare(lc, n) >= 0.6);
    });
  });
}

function renderSearchResults({ reactions = [], compounds = [], tokens }) {
  const box = document.getElementById('searchResults');
  box.innerHTML = '';
  if (!reactions.length && !compounds.length) { box.textContent = 'No matches'; return; }
  if (reactions.length) {
    const h = document.createElement('div');
    h.className = 'results-section';
    h.innerHTML = '<h3>Reactions</h3>';
    reactions.forEach(r => {
      const div = document.createElement('div');
      div.className = 'result-item';
      div.textContent = r.equation || '(no equation)';
      div.addEventListener('click', () => renderReactionDetail(r));
      h.appendChild(div);
    });
    box.appendChild(h);
  }
  if (compounds.length) {
    const h = document.createElement('div');
    h.className = 'results-section';
    h.innerHTML = '<h3>Compounds</h3>';
    compounds.forEach(c => {
      const div = document.createElement('div');
      div.className = 'result-item';
      div.textContent = `${c.formula || '?'} — ${c.name || ''}`.trim();
      div.addEventListener('click', () => renderCompoundDetail(c));
      h.appendChild(div);
    });
    box.appendChild(h);
  }
  if (tokens) {
    const dbg = document.getElementById('searchTokens');
    if (dbg) {
      dbg.innerHTML = `<div class="results-section"><h4>All Reactants (${tokens.reactants.length})</h4><div class="result-item">${tokens.reactants.join(', ')}</div></div>
        <div class="results-section"><h4>All Products (${tokens.products.length})</h4><div class="result-item">${tokens.products.join(', ')}</div></div>`;
    }
  }
}

function showDetail(title, html) {
  const card = document.getElementById('detailCard');
  const t = document.getElementById('detailTitle');
  const body = document.getElementById('detailBody');
  if (!card || !t || !body) return;
  t.textContent = title;
  body.innerHTML = html;
  card.classList.remove('hidden');
}

function hideDetail() {
  const card = document.getElementById('detailCard');
  if (card) card.classList.add('hidden');
}

function renderReactionDetail(r) {
  const reactList = (r.reactants || []).map(x => `<li>${x}</li>`).join('');
  const prodList = (r.products || []).map(x => `<li>${x}</li>`).join('');
  const html = `
    <div class="detail-row"><strong>Equation:</strong> ${r.equation || ''}</div>
    <div class="detail-row"><strong>Difficulty:</strong> ${r.difficulty || 'N/A'}</div>
    <div class="detail-row"><strong>Temp (C):</strong> ${r.temperature_C ?? 'N/A'}</div>
    <div class="detail-row"><strong>Notes:</strong> ${r.notes || ''}</div>
    <div class="detail-row"><strong>Reactants:</strong><ul>${reactList}</ul></div>
    <div class="detail-row"><strong>Products:</strong><ul>${prodList}</ul></div>
  `;
  showDetail('Reaction details', html);
}

function renderCompoundDetail(c) {
  const html = `
    <div class="detail-row"><strong>Formula:</strong> ${c.formula || ''}</div>
    <div class="detail-row"><strong>Name:</strong> ${c.name || ''}</div>
    <div class="detail-row"><strong>Molar mass (g/mol):</strong> ${c.molar_mass_g_mol ?? ''}</div>
    <div class="detail-row"><strong>Melting point (C):</strong> ${c.melting_point_C ?? ''}</div>
    <div class="detail-row"><strong>Boiling point (C):</strong> ${c.boiling_point_C ?? ''}</div>
    <div class="detail-row"><strong>Flags:</strong> ${c.flags || ''}</div>
    <div class="detail-row"><strong>Notes:</strong> ${c.notes || ''}</div>
    <div class="detail-row"><strong>Lethal dose (mg/kg):</strong> ${c.lethal_dose_mg_kg ?? 'N/A'}</div>
  `;
  showDetail('Compound details', html);
}

function lookupReactionInfo(db, reactantsRaw, productsRaw) {
  const norm = arr => [...arr].sort().join('+');
  const sig = norm(reactantsRaw) + '->' + norm(productsRaw);
  return db.reactions.find(r => {
    const parts = String(r.equation || '').replace(/\s+/g, '').split('->');
    if (parts.length < 2) return false;
    return sig === norm(parts[0].split('+')) + '->' + norm(parts[1].split('+'));
  }) || null;
}

function renderReport(data) {
  document.getElementById('results').classList.remove('hidden');
  document.getElementById('balancedDisplay').textContent = data.equation;

  const dbCard = document.getElementById('dbInfo');
  if (data.info) {
    dbCard.classList.remove('hidden');
    document.getElementById('dbDiff').textContent = data.info.difficulty || 'N/A';
    document.getElementById('dbTemp').textContent = data.info.temperature_C || 'N/A';
    document.getElementById('dbNotes').textContent = data.info.notes || 'N/A';
    const diff = (data.info.difficulty || '').toLowerCase();
    dbCard.style.borderLeft = diff === 'easy' ? '5px solid #4caf50' : diff === 'hard' ? '5px solid #ff9800' : '5px solid #f44336';
  } else {
    dbCard.classList.add('hidden');
  }

  const tbody = document.getElementById('stoichTable');
  tbody.innerHTML = '';
  const addRows = (list, role) => {
    list.forEach(item => {
      const tr = document.createElement('tr');
      tr.className = role === 'Reactant' ? 'row-in' : 'row-out';
      let badges = '';
      if (item.flags) badges = item.flags.split(';').map(f => `<span class="badge ${f.trim().toLowerCase()}">${f.trim()}</span>`).join(' ');
      const details = item.name !== 'Unknown' ? `<div class="subtext">${item.name}</div>${badges}` : badges;
      tr.innerHTML = `
        <td>${role}</td>
        <td>${item.coeff}</td>
        <td class="formula">${item.formula}</td>
        <td>${item.molarMass}</td>
        <td>${item.totalMass}</td>
        <td>${item.massFraction}</td>
        <td class="details-cell">${details} <div class="note">${item.notes}</div></td>`;
      tbody.appendChild(tr);
    });
  };
  addRows(data.reactants, 'Reactant');
  addRows(data.products, 'Product');
}

export function setup() {
  const button = document.getElementById('evalBtn');
  const inputEl = document.getElementById('reactionInput');
  const errorDiv = document.getElementById('error');
  const resultsDiv = document.getElementById('results');
  const searchInput = document.getElementById('searchInput');
  const searchBtn = document.getElementById('searchBtn');
  const scopeSelect = document.getElementById('searchScope');
  const modeSelect = document.getElementById('searchMode');
  const logicSelect = document.getElementById('searchLogic');
  const fuzzyChk = document.getElementById('searchFuzzy');
  const tokensToggle = document.getElementById('searchTokensToggle');
  const detailClose = document.getElementById('detailClose');
  const dbPromise = loadDb();

  async function evaluateReaction() {
    errorDiv.classList.add('hidden');
    resultsDiv.classList.add('hidden');
    try {
      const db = await dbPromise;
      const balanced = balance(inputEl.value, db.elements);
      const totalMass = balanced.reactants.reduce((s, r) => s + r.molarMass * r.coeff, 0);
      const reactionInfo = lookupReactionInfo(db, balanced.reactants.map(r => r.formula), balanced.products.map(p => p.formula));
      const toRow = list => list.map(item => {
        const meta = db.compounds.find(c => c.formula === item.formula) || {};
        const total = item.molarMass * item.coeff;
        return {
          formula: item.formula,
          name: meta.name || 'Unknown',
          coeff: item.coeff,
          molarMass: item.molarMass.toFixed(2),
          totalMass: total.toFixed(2),
          massFraction: ((total / totalMass) * 100).toFixed(2) + '%',
          flags: meta.flags || '',
          notes: meta.notes || ''
        };
      });
      renderReport({
        equation: balanced.reactants.map(x => `${x.coeff > 1 ? x.coeff : ''}${x.formula}`).join(' + ') + ' ➝ ' + balanced.products.map(x => `${x.coeff > 1 ? x.coeff : ''}${x.formula}`).join(' + '),
        info: reactionInfo,
        reactants: toRow(balanced.reactants),
        products: toRow(balanced.products)
      });
    } catch (e) {
      errorDiv.textContent = 'Error: ' + e.message;
      errorDiv.classList.remove('hidden');
    }
  }

  function getScope() {
    return scopeSelect ? scopeSelect.value : 'any';
  }

  function getMode() {
    return modeSelect ? modeSelect.value : 'reactions';
  }

  function getLogic() {
    return logicSelect ? logicSelect.value : 'and';
  }

  function getFuzzy() {
    return fuzzyChk ? fuzzyChk.checked : true;
  }

  function getShowTokens() {
    return tokensToggle ? tokensToggle.checked : false;
  }

  async function runSearch() {
    const db = await dbPromise;
    const terms = searchInput.value.split(',').map(t => t.trim()).filter(Boolean);
    const scope = getScope();
    const mode = getMode();
    const logic = getLogic();
    const allowFuzzy = getFuzzy();
    const showTokens = getShowTokens();
    const tokens = uniqueReactionTokens(db);
    let reactions = [], compounds = [];
    if (mode === 'reactions' || mode === 'both') reactions = searchReactions(db, terms, scope, logic, allowFuzzy);
    if (mode === 'compounds' || mode === 'both') compounds = searchCompounds(db, terms, allowFuzzy);
    renderSearchResults({ reactions, compounds, tokens: showTokens ? tokens : null });
  }

  button.addEventListener('click', evaluateReaction);
  if (searchBtn) searchBtn.addEventListener('click', runSearch);
  if (detailClose) detailClose.addEventListener('click', hideDetail);
}

if (typeof window !== 'undefined') {
  window.addEventListener('DOMContentLoaded', () => setup());
}
