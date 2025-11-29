class ElementType {
    constructor(name) {
        this.name = name;
        this.iZ = 0;
        this.neval = 0;
        this.valence = 0;
        this.piMax = 0;
        this.color = 0xFFFFFF;
        this.Rcov = 0.0;
        this.RvdW = 0.0;
        this.EvdW = 0.0;
        this.Quff = 0.0;
        this.Uuff = 0.0;
        this.Vuff = 0.0;
        // QEq
        this.bQEq = false;
        this.Eaff = 0.0;
        this.Ehard = 0.0;
        this.Ra = 0.0;
        this.eta = 0.0;
    }
}

class AtomType {
    constructor(name) {
        this.name = name;
        this.parent_name = "*";
        this.element_name = "";
        this.epair_name = "";
        this.valence = 0;
        this.nepair = 0;
        this.npi = 0;
        this.sym = 0;
        this.Ruff = 0.0;
        this.RvdW = 0.0;
        this.EvdW = 0.0;
        this.Qbase = 0.0;
        this.Hb = 0.0;
        // MMFF
        this.bMMFF = false;
        this.Ass = 0.0;
        this.Asp = 0.0;
        this.Kss = 0.0;
        this.Ksp = 0.0;
        this.Kep = 0.0;
        this.Kpp = 0.0;

        // References
        this.element = null; // ElementType object
        this.parent = null;  // AtomType object
        this.epair = null;   // AtomType object

        // Cached from element for convenience
        this.iZ = 0;
        this.color = 0xFFFFFF;
    }
}

export class MMParams {
    constructor() {
        this.elementTypes = {}; // name -> ElementType
        this.atomTypes = {};    // name -> AtomType
        this.byAtomicNumber = {}; // iZ -> ElementType
    }

    async loadResources(elementPath, atomPath) {
        try {
            const eRes = await fetch(elementPath);
            const eText = await eRes.text();
            this.parseElementTypes(eText);

            const aRes = await fetch(atomPath);
            const aText = await aRes.text();
            this.parseAtomTypes(aText);

            window.logger.info("MMParams loaded successfully.");
        } catch (e) {
            window.logger.error("Failed to load MMParams: " + e);
        }
    }

    parseElementTypes(content) {
        const lines = content.split('\n');
        this.elementTypes = {};
        this.byAtomicNumber = {};

        for (let line of lines) {
            line = line.trim();
            if (!line || line.startsWith('#')) continue;

            const parts = line.split(/\s+/);
            if (parts.length < 12) continue;

            const name = parts[0];
            const et = new ElementType(name);

            try {
                et.iZ = parseInt(parts[1]);
                et.neval = parseInt(parts[2]);
                et.valence = parseInt(parts[3]);
                et.piMax = parseInt(parts[4]);

                let colStr = parts[5];
                if (colStr.startsWith('0x')) {
                    et.color = parseInt(colStr, 16);
                } else {
                    et.color = parseInt(colStr);
                }

                et.Rcov = parseFloat(parts[6]);
                et.RvdW = parseFloat(parts[7]);
                et.EvdW = parseFloat(parts[8]);
                et.Quff = parseFloat(parts[9]);
                et.Uuff = parseFloat(parts[10]);
                et.Vuff = parseFloat(parts[11]);

                if (parts.length >= 16) {
                    et.bQEq = true;
                    et.Eaff = parseFloat(parts[12]);
                    et.Ehard = parseFloat(parts[13]);
                    et.Ra = parseFloat(parts[14]);
                    et.eta = parseFloat(parts[15]);
                }

                this.elementTypes[name] = et;
                this.byAtomicNumber[et.iZ] = et;

            } catch (e) {
                console.warn("Error parsing element line:", line, e);
            }
        }
        window.logger.info(`Parsed ${Object.keys(this.elementTypes).length} element types.`);
        return this.elementTypes;
    }

    parseAtomTypes(content) {
        const lines = content.split('\n');
        this.atomTypes = {};

        // Pass 1: Create objects
        for (let line of lines) {
            line = line.trim();
            if (!line || line.startsWith('#')) continue;
            if (line.startsWith('*')) continue; // Skip header line starting with *

            const parts = line.split(/\s+/);
            if (parts.length < 5) continue;

            const name = parts[0];
            const at = new AtomType(name);
            at.parent_name = parts[1];
            at.element_name = parts[2];
            at.epair_name = parts[3];

            try {
                at.valence = parseInt(parts[4]);
                at.nepair = parseInt(parts[5]);
                at.npi = parseInt(parts[6]);
                at.sym = parseInt(parts[7]);
                at.Ruff = parseFloat(parts[8]);
                at.RvdW = parseFloat(parts[9]);
                at.EvdW = parseFloat(parts[10]);
                at.Qbase = parseFloat(parts[11]);
                at.Hb = parseFloat(parts[12]);

                if (parts.length >= 19) {
                    at.bMMFF = true;
                    at.Ass = parseFloat(parts[13]);
                    at.Asp = parseFloat(parts[14]);
                    at.Kss = parseFloat(parts[15]);
                    at.Ksp = parseFloat(parts[16]);
                    at.Kep = parseFloat(parts[17]);
                    at.Kpp = parseFloat(parts[18]);
                }

                this.atomTypes[name] = at;

            } catch (e) {
                console.warn("Error parsing atom type line:", line, e);
            }
        }

        // Pass 2: Resolve references
        for (const key in this.atomTypes) {
            const at = this.atomTypes[key];

            // Element
            if (this.elementTypes[at.element_name]) {
                at.element = this.elementTypes[at.element_name];
                at.iZ = at.element.iZ;
                at.color = at.element.color;
            }

            // Parent
            if (at.parent_name !== '*' && this.atomTypes[at.parent_name]) {
                at.parent = this.atomTypes[at.parent_name];
            }

            // Epair
            if (at.epair_name !== '*' && this.atomTypes[at.epair_name]) {
                at.epair = this.atomTypes[at.epair_name];
            }
        }

        window.logger.info(`Parsed ${Object.keys(this.atomTypes).length} atom types.`);
        return this.atomTypes;
    }

    // --- Helpers for Renderer ---

    getColor(iZ) {
        const et = this.byAtomicNumber[iZ];
        if (et) {
            // Convert integer color to [r, g, b] normalized
            const r = ((et.color >> 16) & 255) / 255.0;
            const g = ((et.color >> 8) & 255) / 255.0;
            const b = (et.color & 255) / 255.0;
            return [r, g, b];
        }
        return [1.0, 0.0, 1.0]; // Default Magenta
    }

    getRadius(iZ) {
        const et = this.byAtomicNumber[iZ];
        if (et) return et.RvdW;
        return 1.5; // Default
    }

    getElementName(iZ) {
        const et = this.byAtomicNumber[iZ];
        return et ? et.name : "X";
    }
}
