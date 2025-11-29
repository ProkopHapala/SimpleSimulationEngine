export class Selection {
    constructor() {
        this.kind = null;      // e.g. 'vert', 'edge', 'face'
        this.bUnique = true;   // if true, avoid duplicates
        this.vec = [];         // ordered list of selected indices
        this.map = new Map();  // key -> index into vec, for fast membership
    }

    add(i) {
        if (this.bUnique) {
            if (!this.map.has(i)) {
                const idx = this.vec.length;
                this.vec.push(i);
                this.map.set(i, idx);
                return true;
            }
        } else {
            this.vec.push(i);
            return true;
        }
        return false;
    }

    remove(i) {
        const it = this.map.get(i);
        if (it === undefined) return false;
        const idx = it;
        if (idx >= 0 && idx < this.vec.length) {
            this.vec[idx] = -1; // mark as removed (keeps indices stable)
        }
        this.map.delete(i);
        return true;
    }

    toggle(i) {
        const it = this.map.get(i);
        if (it !== undefined) {
            const idx = it;
            if (idx >= 0 && idx < this.vec.length) {
                this.vec[idx] = -1;
            }
            this.map.delete(i);
            return -1;
        } else {
            const idx = this.vec.length;
            this.vec.push(i);
            this.map.set(i, idx);
            return 1;
        }
    }

    contains(i) {
        return this.map.has(i);
    }

    size() {
        return this.map.size;
    }

    data() {
        return this.vec;
    }

    sort() {
        // Remove tombstones (-1) before sorting
        this.vec = this.vec.filter((v) => v >= 0).sort((a, b) => a - b);
        this.map.clear();
        for (let idx = 0; idx < this.vec.length; idx++) {
            this.map.set(this.vec[idx], idx);
        }
    }

    clear() {
        this.vec.length = 0;
        this.map.clear();
    }

    insert(container) {
        for (const i of container) {
            this.add(i);
        }
    }

    insertArray(n, arr) {
        for (let i = 0; i < n; i++) {
            this.add(arr[i]);
        }
    }

    subtract(container, bInv = false) {
        const newVec = [];
        let nerased = 0;
        for (const i of this.vec) {
            if (i < 0) continue;
            const inOther = container.contains ? container.contains(i) : container.has(i);
            if ((inOther && !bInv) || (!inOther && bInv)) {
                this.map.delete(i);
                nerased++;
            } else {
                newVec.push(i);
            }
        }
        this.vec = newVec;
        // Rebuild map indices
        this.map.clear();
        for (let idx = 0; idx < this.vec.length; idx++) {
            this.map.set(this.vec[idx], idx);
        }
        return nerased;
    }

    intersectWith(container) {
        return this.subtract(container, true);
    }

    selectByPredicate(range, predicate) {
        this.clear();
        let index = 0;
        for (const element of range) {
            if (predicate(element)) {
                this.add(index);
            }
            index++;
        }
        return this.vec.length;
    }

    toSortedArray() {
        const arr = [];
        for (const i of this.vec) {
            if (i >= 0) arr.push(i);
        }
        arr.sort((a, b) => a - b);
        return arr;
    }
}

export class SelectionBanks {
    constructor(n) {
        this.icurSelection = 0;
        this.selections = [];
        for (let i = 0; i < n; i++) {
            this.selections.push(new Selection());
        }
        this.curSelection = this.selections[0];
    }

    pickSelection(i) {
        const n = this.selections.length;
        const idx = Math.max(0, Math.min(n - 1, i | 0));
        this.icurSelection = idx;
        this.curSelection = this.selections[idx];
        return this.curSelection;
    }

    nextSelection() {
        let i = this.icurSelection + 1;
        if (i >= this.selections.length) i = 0;
        return this.pickSelection(i);
    }
}

// Also expose on window for legacy global-script users
if (typeof window !== 'undefined') {
    window.Selection = Selection;
    window.SelectionBanks = SelectionBanks;
}
