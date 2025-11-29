export class MoleculeSystem {
    constructor(capacity = 1000) {
        this.capacity = capacity;
        this.nAtoms = 0;

        // Structure of Arrays (SoA)
        this.pos = new Float32Array(this.capacity * 3); // x, y, z
        this.types = new Uint8Array(this.capacity);     // element type (atomic number)

        // Bonds (Simple list for now, can be optimized later)
        this.bonds = []; // Array of [id1, id2]

        // Dirty flags
        this.isDirty = true;

        // Selection
        this.selection = new Set();

        // Neighbor List (Adjacency List)
        // Array of Set<int> or Array<int>
        this.neighborList = [];
    }

    select(id, mode = 'replace') {
        if (mode === 'replace') {
            this.selection.clear();
            this.selection.add(id);
        } else if (mode === 'add') {
            this.selection.add(id);
        } else if (mode === 'subtract') {
            this.selection.delete(id);
        }
        this.isDirty = true; // Trigger redraw to update colors
    }

    clearSelection() {
        this.selection.clear();
        this.isDirty = true;
    }

    addAtom(x, y, z, type) {
        if (this.nAtoms >= this.capacity) {
            this.resize(this.capacity * 2);
        }

        const i = this.nAtoms;
        this.pos[i * 3] = x;
        this.pos[i * 3 + 1] = y;
        this.pos[i * 3 + 2] = z;
        this.types[i] = type;

        // Initialize neighbor list for new atom
        if (this.neighborList.length <= i) {
            this.neighborList.push([]);
        } else {
            this.neighborList[i] = [];
        }

        this.nAtoms++;
        this.isDirty = true;
        return i;
    }

    addBond(id1, id2) {
        this.bonds.push([id1, id2]);

        // Update neighbor list
        if (this.neighborList[id1]) this.neighborList[id1].push(id2);
        if (this.neighborList[id2]) this.neighborList[id2].push(id1);

        this.isDirty = true;
    }

    updateNeighborList() {
        this.neighborList = new Array(this.nAtoms).fill(null).map(() => []);
        for (const [id1, id2] of this.bonds) {
            if (id1 < this.nAtoms && id2 < this.nAtoms) {
                this.neighborList[id1].push(id2);
                this.neighborList[id2].push(id1);
            }
        }
    }

    recalculateBonds(mmParams = null) {
        this.bonds = [];

        // Optimization: Use Spatial Hash if N > 1000 (TODO)

        const defaultRcut = 1.6;
        const defaultRcut2 = defaultRcut * defaultRcut;
        const bondFactor = 1.3; // Tolerance factor for covalent bonds

        for (let i = 0; i < this.nAtoms; i++) {
            const type1 = this.types[i];
            let r1 = 0.7; // Default ~Carbon
            if (mmParams && mmParams.byAtomicNumber[type1]) {
                r1 = mmParams.byAtomicNumber[type1].Rcov;
            }

            for (let j = i + 1; j < this.nAtoms; j++) {
                const type2 = this.types[j];
                let r2 = 0.7;
                let rCut2 = defaultRcut2;

                if (mmParams) {
                    if (mmParams.byAtomicNumber[type2]) {
                        r2 = mmParams.byAtomicNumber[type2].Rcov;
                    }
                    const rSum = (r1 + r2) * bondFactor;
                    rCut2 = rSum * rSum;
                }

                const dx = this.pos[i * 3] - this.pos[j * 3];
                const dy = this.pos[i * 3 + 1] - this.pos[j * 3 + 1];
                const dz = this.pos[i * 3 + 2] - this.pos[j * 3 + 2];

                const dist2 = dx * dx + dy * dy + dz * dz;
                if (dist2 < rCut2) {
                    this.bonds.push([i, j]);
                }
            }
        }
        this.updateNeighborList();
        this.isDirty = true;
        window.logger.info(`Recalculated bonds. Found ${this.bonds.length} bonds.`);
    }

    deleteSelectedAtoms() {
        if (this.selection.size === 0) return;

        const toDelete = Array.from(this.selection).sort((a, b) => b - a); // Delete from end to avoid index shift issues? 
        // Actually, swap-remove changes indices, so we need to be careful.
        // Better strategy: 
        // 1. Mark atoms to keep.
        // 2. Rebuild arrays.
        // 3. Re-map selection (clear it).

        const oldToNew = new Int32Array(this.nAtoms).fill(-1);
        let newCount = 0;

        // 1. Calculate new indices
        for (let i = 0; i < this.nAtoms; i++) {
            if (!this.selection.has(i)) {
                oldToNew[i] = newCount;
                newCount++;
            }
        }

        // 2. Compact Arrays
        const newPos = new Float32Array(this.capacity * 3);
        const newTypes = new Uint8Array(this.capacity);

        for (let i = 0; i < this.nAtoms; i++) {
            if (oldToNew[i] !== -1) {
                const newIdx = oldToNew[i];
                newPos[newIdx * 3] = this.pos[i * 3];
                newPos[newIdx * 3 + 1] = this.pos[i * 3 + 1];
                newPos[newIdx * 3 + 2] = this.pos[i * 3 + 2];
                newTypes[newIdx] = this.types[i];
            }
        }

        this.pos = newPos;
        this.types = newTypes;

        // 3. Rebuild Bonds
        const newBonds = [];
        for (const [id1, id2] of this.bonds) {
            if (oldToNew[id1] !== -1 && oldToNew[id2] !== -1) {
                newBonds.push([oldToNew[id1], oldToNew[id2]]);
            }
        }
        this.bonds = newBonds;

        this.nAtoms = newCount;
        this.selection.clear();
        this.updateNeighborList();
        this.isDirty = true;

        window.logger.info(`Deleted atoms. New count: ${this.nAtoms}`);
    }

    resize(newCapacity) {
        window.logger.info(`Resizing MoleculeSystem to ${newCapacity}`);
        const newPos = new Float32Array(newCapacity * 3);
        const newTypes = new Uint8Array(newCapacity);

        newPos.set(this.pos);
        newTypes.set(this.types);

        this.pos = newPos;
        this.types = newTypes;
        this.capacity = newCapacity;
    }

    clear() {
        this.nAtoms = 0;
        this.bonds = [];
        this.isDirty = true;
    }
}
