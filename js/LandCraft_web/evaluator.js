// evaluator.js (ES module)
import ELEMENT_DICT from './elements.js';

/**
 * Parses chemical formulas (e.g., "H2O", "Ca(OH)2")
 */
export function parseFormula(formula) {
    const counts = {};
    const stack = [{}];
    let i = 0;
    while (i < formula.length) {
        let char = formula[i];
        if (char === '(') {
            stack.push({});
            i++;
        } else if (char === ')') {
            i++;
            let numStr = "";
            while (i < formula.length && /\d/.test(formula[i])) numStr += formula[i++];
            const mult = numStr === "" ? 1 : parseInt(numStr);
            const popped = stack.pop();
            const top = stack[stack.length - 1];
            for (let el in popped) top[el] = (top[el] || 0) + popped[el] * mult;
        } else if (/[A-Z]/.test(char)) {
            let el = char;
            i++;
            if (i < formula.length && /[a-z]/.test(formula[i])) el += formula[i++];
            let numStr = "";
            while (i < formula.length && /\d/.test(formula[i])) numStr += formula[i++];
            const count = numStr === "" ? 1 : parseInt(numStr);
            const top = stack[stack.length - 1];
            top[el] = (top[el] || 0) + count;
        } else i++;
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

/**
 * Solves the equation via trial-and-error for small integer coefficients
 */
export function balance(reactionStr, elementDict = ELEMENT_DICT) {
    const sides = reactionStr.split('->');
    if (sides.length !== 2) throw new Error("Invalid format. Use 'Reactants -> Products'");
    
    const leftTerms = sides[0].split('+').map(parseTerm).filter(t => t.formula);
    const rightTerms = sides[1].split('+').map(parseTerm).filter(t => t.formula);
    const allTerms = [...leftTerms, ...rightTerms];
    const compounds = allTerms.map(t => ({ formula: t.formula, atoms: parseFormula(t.formula), base: t.base }));
    const elements = [...new Set(compounds.flatMap(c => Object.keys(c.atoms)))];

    const matrix = elements.map(el => compounds.map((c, idx) => (idx < leftTerms.length ? 1 : -1) * (c.atoms[el] || 0) * c.base));
    const coeffs = nullspaceIntegerSolution(matrix);

    const reactants = leftTerms.map((t, i) => ({
        formula: compounds[i].formula,
        coeff: coeffs[i] * compounds[i].base,
        atoms: compounds[i].atoms
    }));
    const products = rightTerms.map((t, i) => ({
        formula: compounds[i + leftTerms.length].formula,
        coeff: coeffs[i + leftTerms.length] * compounds[i + leftTerms.length].base,
        atoms: compounds[i + leftTerms.length].atoms
    }));

    const g = reactants.concat(products).reduce((acc, x) => gcd(acc, x.coeff), reactants[0]?.coeff || 1);
    reactants.forEach(r => r.coeff = r.coeff / g);
    products.forEach(p => p.coeff = p.coeff / g);

    return { reactants, products, elements: elementDict };
}

/**
 * Main execution function
 */
export function run(input) {
    try {
        const result = balance(input, ELEMENT_DICT);
        const all = [...result.reactants, ...result.products];
        
        console.log(`\nEquation: ${input}`);
        console.log(`Balanced: ` + 
            result.reactants.map(r => `${r.coeff > 1 ? r.coeff : ''}${r.formula}`).join(' + ') + " -> " +
            result.products.map(p => `${p.coeff > 1 ? p.coeff : ''}${p.formula}`).join(' + ')
        );

        console.log(`\nMass Analysis:`);
        let baseMass = 0;
        all.forEach((comp, idx) => {
            const molarMass = Object.entries(comp.atoms).reduce((acc, [sym, count]) => {
                if (!ELEMENT_DICT[sym]) throw new Error(`Unknown element: ${sym}`);
                return acc + (ELEMENT_DICT[sym].mass * count);
            }, 0);
            
            const totalMass = molarMass * comp.coeff;
            if (idx === 0) baseMass = totalMass;

            console.log(`[${comp.formula}]`);
            console.log(`  Molar Mass : ${molarMass.toFixed(3)} g/mol`);
            console.log(`  Stoich Mass: ${totalMass.toFixed(3)} g`);
            console.log(`  Mass Ratio : ${(totalMass / baseMass).toFixed(4)}`);
        });
    } catch (e) {
        console.error("Error:", e.message);
    }
}

if (typeof process !== 'undefined' && process.argv && process.argv[1] && import.meta && import.meta.url && process.argv[1].endsWith('evaluator.js')) {
    const input = process.argv[2] || "C3H8 + O2 -> CO2 + H2O";
    run(input);
}