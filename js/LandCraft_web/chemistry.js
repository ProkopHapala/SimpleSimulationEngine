import fs from 'fs';
import path from 'path';

const __dirname = path.dirname(new URL(import.meta.url).pathname);

// --- 1. CSV Utility ---
export class CSVManager {
    static load(filepath) {
        if (!fs.existsSync(filepath)) return [];
        const content = fs.readFileSync(filepath, 'utf-8').trim();
        const lines = content.split('\n');
        if (lines.length < 2) return [];
        const headers = lines[0].split(',').map(h => h.trim());
        
        // Handle CSV parsing considering quoted strings (e.g. "acid;oxidant")
        return lines.slice(1).map(line => {
            const values = [];
            let current = '';
            let inQuotes = false;
            for (let char of line) {
                if (char === '"') inQuotes = !inQuotes;
                else if (char === ',' && !inQuotes) { values.push(current.trim()); current = ''; }
                else current += char;
            }
            values.push(current.trim());
            
            return headers.reduce((obj, header, i) => {
                let val = values[i];
                if (!isNaN(val) && val !== '') val = parseFloat(val);
                obj[header] = val;
                return obj;
            }, {});
        });
    }
}

// ==========================================
// 1. MODIFIED UTILITIES (Browser Compatible)
// ==========================================

export class CSVManager {
    // Modified to parse string content directly (since we fetch it elsewhere)
    static parse(content) {
        if (!content) return [];
        const lines = content.trim().split('\n');
        if (lines.length < 2) return [];
        const headers = lines[0].split(',').map(h => h.trim());
        
        return lines.slice(1).map(line => {
            const values = [];
            let current = '';
            let inQuotes = false;
            for (let char of line) {
                if (char === '"') inQuotes = !inQuotes;
                else if (char === ',' && !inQuotes) { values.push(current.trim()); current = ''; }
                else current += char;
            }
            values.push(current.trim());
            
            return headers.reduce((obj, header, i) => {
                let val = values[i];
                if (!isNaN(val) && val !== '') val = parseFloat(val);
                obj[header] = val;
                return obj;
            }, {});
        });
    }
}

// ... Keep Element, Compound, and Reaction classes exactly as they were ...
// ... (Paste your existing classes here) ...


// ==========================================
// 2. NEW SEARCH FUNCTIONALITY
// ==========================================

export class ReactionSearcher {
    /**
     * Search reactions based on components (names or formulas)
     * @param {Object} db - The loaded database object
     * @param {Array<string>} terms - Array of search strings (e.g. ["CaO", "Water"])
     * @param {string} scope - "reactant", "product", or "any"
     */
    static search(db, terms, scope = 'any') {
        // 1. Resolve terms to possible formulas
        // E.g., "Water" -> ["Water", "H2O"]
        const resolvedTerms = terms.map(term => this._resolveSynonyms(term, db));

        return db.reactions.filter(reaction => {
            // Parse the equation string into two arrays of trimmed strings
            const parts = reaction.equation.split(/->|--.*-->/);
            if (parts.length < 2) return false;

            // Extract formulas/names from the reaction string (remove coefficients)
            const clean = (side) => side.split('+').map(s => {
                s = s.trim();
                // Remove leading coefficient (e.g. "2 H2O" -> "H2O")
                return s.replace(/^\d+\s+/, ''); 
            });

            const reactants = clean(parts[0]);
            const products = clean(parts[1]);

            // Define where to look based on scope
            let targetComponents = [];
            if (scope === 'reactant') targetComponents = reactants;
            else if (scope === 'product') targetComponents = products;
            else targetComponents = [...reactants, ...products];

            // 2. Check: Does the reaction contain ALL search terms?
            // We use 'every' so that searching "A, B" finds reactions having BOTH A and B.
            return resolvedTerms.every(possibilities => {
                return targetComponents.some(componentInReaction => {
                    // componentInReaction might be "H2O" or "butter"
                    // possibilities might be ["Water", "H2O"]
                    return possibilities.some(p => p.toLowerCase() === componentInReaction.toLowerCase());
                });
            });
        });
    }

    /**
     * Helper: Maps a name (e.g., "Lime") to its formula ("CaO") using the compounds DB,
     * but also keeps the original name in case the reaction uses the name directly.
     */
    static _resolveSynonyms(term, db) {
        const possibilities = [term]; // Always include the original search term
        
        // Find if this term is a name in compounds.csv
        const compoundEntry = db.compounds.find(c => 
            c.name && c.name.toLowerCase().includes(term.toLowerCase())
        );

        if (compoundEntry) {
            possibilities.push(compoundEntry.formula);
        }
        
        // Conversely, if term is a formula, try to find its name? 
        // (Optional, usually less needed for reaction search)
        
        return possibilities;
    }
}

// ==========================================
// 3. UPDATED BROWSER LOADER
// ==========================================

export async function loadDatabase() {
    // Helper to fetch and text()
    const fetchTable = async (filename) => {
        try {
            const res = await fetch(`data/${filename}`);
            if (!res.ok) throw new Error(`Failed to load ${filename}`);
            const text = await res.text();
            return CSVManager.parse(text);
        } catch (e) {
            console.error(e);
            return [];
        }
    };

    const [elementsRaw, compounds, reactions] = await Promise.all([
        fetchTable('elements.csv'),
        fetchTable('compounds.csv'),
        fetchTable('reactions.csv')
    ]);
    // Convert elements array to dict
    const elements = elementsRaw.reduce((acc, row) => { acc[row.symbol] = row; return acc; }, {});

    return { elements, compounds, reactions };
}

// --- 2. Chemical Classes ---
export class Element {
    constructor(row) {
        Object.assign(this, row);
    }
}

export class Compound {
    constructor(formula, db) {
        this.formula = formula;
        this.atoms = this._parse(formula);
        this.molarMass = this._calcMolarMass(db.elements);
        // Find metadata in compounds.csv
        this.meta = db.compounds.find(c => c.formula === formula) || {};
    }

    _parse(formula) {
        const counts = {};
        const stack = [{}];
        let i = 0;
        while (i < formula.length) {
            let char = formula[i];
            if (char === '(') { stack.push({}); i++; }
            else if (char === ')') {
                i++;
                let numStr = "";
                while (i < formula.length && /\d/.test(formula[i])) numStr += formula[i++];
                const mult = numStr === "" ? 1 : parseInt(numStr);
                const popped = stack.pop();
                const top = stack[stack.length - 1];
                for (let el in popped) top[el] = (top[el] || 0) + popped[el] * mult;
            } else if (/[A-Z]/.test(char)) {
                let el = char; i++;
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

    _calcMolarMass(elDb) {
        return Object.entries(this.atoms).reduce((sum, [sym, count]) => {
            const el = elDb[sym];
            if (!el) throw new Error(`Unknown element: ${sym}`);
            return sum + (el.mass || el.molar_mass_g_mol || 0) * count;
        }, 0);
    }
}

export class Reaction {
    constructor(str, db) {
        this.db = db;
        this.inputStr = str;
        // Clean input
        const cleanStr = str.replace(/\s+/g, '');
        const parts = cleanStr.split(/->|--.*-->/); // Handle -> or --catalyst-->
        if (parts.length < 2) throw new Error("Invalid reaction format");
        
        this.reactantsRaw = parts[0].split('+').filter(x=>x);
        this.productsRaw = parts[1].split('+').filter(x=>x);
        
        this.balanced = this._balance();
        this.info = this._findInfo();
    }

    _balance() {
        const allFormulas = [...this.reactantsRaw, ...this.productsRaw];
        const compounds = allFormulas.map(f => new Compound(f, this.db));
        const elements = [...new Set(compounds.flatMap(c => Object.keys(c.atoms)))];

        // Matrix Setup
        const matrix = elements.map(el => {
            return compounds.map((c, j) => (j < this.reactantsRaw.length ? (c.atoms[el]||0) : -(c.atoms[el]||0)));
        });

        // Simple Solver (1-20)
        const coeffs = new Array(compounds.length).fill(1);
        const solve = (idx) => {
            if (idx === compounds.length) return matrix.every(row => row.reduce((s, v, i) => s + v * coeffs[i], 0) === 0);
            for (let v = 1; v <= 20; v++) {
                coeffs[idx] = v;
                if (solve(idx + 1)) return true;
            }
            return false;
        };

        if (!solve(0)) throw new Error("Could not balance automatically (too complex or invalid)");

        return {
            reactants: this.reactantsRaw.map((f, i) => ({ 
                comp: compounds[i], 
                coeff: coeffs[i] 
            })),
            products: this.productsRaw.map((f, i) => ({ 
                comp: compounds[i + this.reactantsRaw.length], 
                coeff: coeffs[i + this.reactantsRaw.length] 
            }))
        };
    }

    _findInfo() {
        // Simple fuzzy check against reactions.csv
        // We strip spaces and compare "reactants->products"
        const normalize = (arr) => arr.sort().join('+');
        const mySig = normalize(this.reactantsRaw) + "->" + normalize(this.productsRaw);
        
        // Find best match in DB
        return this.db.reactions.find(r => {
            const rParts = r.equation.replace(/\s+/g, '').split(/->/);
            if(rParts.length < 2) return false;
            const dbSig = normalize(rParts[0].split('+')) + "->" + normalize(rParts[1].split('+'));
            return mySig === dbSig;
        }) || null;
    }

    getReport() {
        const all = [...this.balanced.reactants, ...this.balanced.products];
        // Calculate total mass of reactants to normalize fractions
        const totalMass = this.balanced.reactants.reduce((sum, item) => sum + (item.comp.molarMass * item.coeff), 0);

        const formatPart = (list) => list.map(item => ({
            formula: item.comp.formula,
            name: item.comp.meta.name || "Unknown",
            coeff: item.coeff,
            molarMass: item.comp.molarMass.toFixed(2),
            totalMass: (item.comp.molarMass * item.coeff).toFixed(2),
            massFraction: ((item.comp.molarMass * item.coeff) / totalMass * 100).toFixed(2) + "%",
            flags: item.comp.meta.flags || "",
            notes: item.comp.meta.notes || ""
        }));

        return {
            equation: this.balanced.reactants.map(x => `${x.coeff > 1 ? x.coeff : ''}${x.comp.formula}`).join(' + ') 
                      + " ➝ " + 
                      this.balanced.products.map(x => `${x.coeff > 1 ? x.coeff : ''}${x.comp.formula}`).join(' + '),
            info: this.info, // Database info if found
            reactants: formatPart(this.balanced.reactants),
            products: formatPart(this.balanced.products)
        };
    }
}

// --- 3. Database Loader ---
export async function loadDatabase() {
    const elementsRaw = CSVManager.load(path.join(__dirname, 'data/elements.csv'));
    // Convert elements array to dict
    const elements = elementsRaw.reduce((acc, row) => { acc[row.symbol] = row; return acc; }, {});
    return {
        elements,
        compounds: CSVManager.load(path.join(__dirname, 'data/compounds.csv')),
        reactions: CSVManager.load(path.join(__dirname, 'data/reactions.csv'))
    };
}

export { CSVManager, Element, Compound, Reaction, loadDatabase };