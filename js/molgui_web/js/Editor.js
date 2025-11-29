import { Logger } from '../../common_js/Logger.js';

export class Editor {
    constructor(scene, camera, renderer, system, molRenderer) {
        this.scene = scene;
        this.camera = camera;
        this.renderer = renderer;
        this.system = system;
        this.molRenderer = molRenderer;

        this.raycaster = new THREE.Raycaster();
        this.mouse = new THREE.Vector2();

        // Reusable objects for performance
        this.tempVec = new THREE.Vector3();
        this.rayOrigin = new THREE.Vector3();
        this.rayDir = new THREE.Vector3();

        this.isDragging = false;
        this.startPos = { x: 0, y: 0 };
        this.selectionBox = document.createElement('div');
        this.selectionBox.className = 'selection-box';
        document.body.appendChild(this.selectionBox);

        // Gizmo Setup
        this.gizmo = null;
        this.gizmoMode = 'translate'; // translate, rotate, scale
        this.selectionLocked = false;

        this.selectedElement = 6; // Default Carbon
        this.selectedAtomType = 'C';

        // Gizmo Visualsection Lock flag

        this.proxyObject = new THREE.Object3D();
        this.scene.add(this.proxyObject);

        if (THREE.TransformControls) {
            this.gizmo = new THREE.TransformControls(this.camera, this.renderer.domElement);
            this.gizmo.addEventListener('change', () => this.onGizmoDragUpdate());
            this.gizmo.addEventListener('dragging-changed', (event) => {
                const isDragging = event.value;
                // Disable OrbitControls when dragging gizmo
                if (window.app && window.app.controls) {
                    window.app.controls.enabled = !isDragging;
                }

                if (isDragging) {
                    this.onGizmoDragStart();
                } else {
                    this.initialAtomStates = null; // Clear state on drag end
                }
            });
            this.scene.add(this.gizmo);
            this.gizmo.attach(this.proxyObject);
            this.gizmo.visible = false; // Hidden by default until selection
        } else {
            console.error("THREE.TransformControls not found!");
        }

        window.editorInitialized = true;
        this.initEvents();
    }

    initEvents() {
        const canvas = this.renderer.domElement;
        // console.log(`[Editor] InitEvents. Window Size: ${window.innerWidth}x${window.innerHeight}`);

        // Removed capture: true to allow Gizmo (TransformControls) to handle events first if possible
        // But we attach to document/window, while Gizmo attaches to canvas.
        // If we use capture: false (bubbling), Gizmo (on canvas) fires first, then bubble to document.
        // So removing capture: true is correct.
        canvas.addEventListener('pointerdown', this.onMouseDown.bind(this));
        window.addEventListener('pointermove', this.onMouseMove.bind(this));
        window.addEventListener('pointerup', this.onMouseUp.bind(this));

        window.editorEventsAttached = true;
    }

    // --- Gizmo / Shortcut API ---

    toggleGizmo(state) {
        if (!this.gizmo) return;
        // If state is provided, use it, otherwise toggle
        if (state !== undefined) {
            this.gizmoUserEnabled = state;
        } else {
            this.gizmoUserEnabled = !this.gizmoUserEnabled;
        }
        this.updateGizmo();
        window.logger.info(`Gizmo ${this.gizmoUserEnabled ? 'Enabled' : 'Disabled'}`);
    }

    setGizmoMode(mode) {
        if (!this.gizmo) return;
        this.gizmo.setMode(mode);
        window.logger.info(`Gizmo Mode: ${mode}`);
    }

    clearSelection() {
        if (this.selectionLocked) return; // Respect lock
        this.system.clearSelection();
        this.initialAtomStates = null; // Fix: Clear gizmo state
        this.molRenderer.update();
        this.molRenderer.update();
        this.updateGizmo();
        if (this.onSelectionChange) this.onSelectionChange();
        window.logger.info("Selection Cleared");
    }

    deleteSelection() {
        if (window.VERBOSITY_LEVEL >= Logger.DEBUG) {
            window.logger.debug(`[Editor] deleteSelection called. Selection size: ${this.system.selection.size}, Locked: ${this.selectionLocked}`);
        }

        if (this.selectionLocked) {
            window.logger.warn("[Editor] Cannot delete: Selection is locked.");
            return;
        }
        if (this.system.selection.size === 0) {
            window.logger.warn("[Editor] Cannot delete: Nothing selected.");
            return;
        }

        this.system.deleteSelectedAtoms();
        this.initialAtomStates = null;
        this.molRenderer.update();
        this.molRenderer.update(); // Double update to ensure buffers are swapped/cleared if needed
        this.updateGizmo();
        if (this.onSelectionChange) this.onSelectionChange();
    }

    recalculateBonds() {
        this.system.recalculateBonds();
        this.molRenderer.update();
        window.logger.info("Bonds Recalculated");
    }

    addAtom() {
        if (this.selectionLocked) return;
        // Add atom at origin or near selection?
        // Let's add at 0,0,0 for now, or center of view?
        // Center of view is better.

        // For now, fixed position + offset to avoid overlap
        const x = (Math.random() - 0.5) * 2;
        const y = (Math.random() - 0.5) * 2;
        const z = (Math.random() - 0.5) * 2;

        const type = this.selectedElement || 6; // Use selected element

        const id = this.system.addAtom(x, y, z, type);

        // Auto-bond to neighbors?
        // this.system.recalculateBonds(); // Maybe too expensive to do globally?
        // Just add single atom.

        // Select the new atom
        this.system.select(id, 'replace');
        this.updateGizmo();
        this.molRenderer.update();

        window.logger.info(`Added atom ${id} (Type: ${type})`);
    }

    // --- Internal Logic ---

    updateGizmo() {
        if (!this.gizmo) return;

        // Only show if user enabled AND selection exists
        if (!this.gizmoUserEnabled || this.system.selection.size === 0) {
            this.gizmo.visible = false;
            this.gizmo.enabled = false;
            return;
        }

        // Calculate Centroid
        let cx = 0, cy = 0, cz = 0;
        let count = 0;
        for (const id of this.system.selection) {
            cx += this.system.pos[id * 3];
            cy += this.system.pos[id * 3 + 1];
            cz += this.system.pos[id * 3 + 2];
            count++;
        }

        if (count > 0) {
            cx /= count;
            cy /= count;
            cz /= count;

            this.proxyObject.position.set(cx, cy, cz);
            this.proxyObject.rotation.set(0, 0, 0);
            this.proxyObject.scale.set(1, 1, 1);
            this.proxyObject.updateMatrixWorld();

            this.gizmo.visible = true;
            this.gizmo.enabled = true;
        }
    }

    onGizmoDragStart() {
        // Capture relative positions
        this.initialAtomStates = [];
        const inverseProxy = new THREE.Matrix4().copy(this.proxyObject.matrixWorld).invert();

        for (const id of this.system.selection) {
            const vec = new THREE.Vector3(
                this.system.pos[id * 3],
                this.system.pos[id * 3 + 1],
                this.system.pos[id * 3 + 2]
            );
            // Store local position relative to proxy
            vec.applyMatrix4(inverseProxy);
            this.initialAtomStates.push({ id: id, localPos: vec });
        }
    }

    onGizmoDragUpdate() {
        if (!this.initialAtomStates) return;

        const proxyMatrix = this.proxyObject.matrixWorld;
        const vec = new THREE.Vector3();

        for (const state of this.initialAtomStates) {
            vec.copy(state.localPos);
            vec.applyMatrix4(proxyMatrix);

            const i = state.id;
            this.system.pos[i * 3] = vec.x;
            this.system.pos[i * 3 + 1] = vec.y;
            this.system.pos[i * 3 + 2] = vec.z;
        }

        this.system.isDirty = true;
        this.molRenderer.update();
        this.molRenderer.update();
    }

    onMouseDown(e) {
        // console.log(`[Editor] MouseDown: Button=${e.button} Client=(${e.clientX}, ${e.clientY})`);
        window.lastMouseDownEvent = {
            button: e.button,
            clientX: e.clientX,
            clientY: e.clientY,
            time: Date.now(),
            target: e.target.tagName
        };

        // ONLY Allow Left Mouse Button (0) for Selection/Editor
        if (e.button !== 0) return;

        // Fix: Do not start selection if dragging gizmo
        if (this.gizmo && this.gizmo.dragging) return;

        // Fix: Check Selection Lock
        if (this.selectionLocked) return;

        this.isDragging = true;
        this.startPos.x = e.clientX;
        this.startPos.y = e.clientY;

        // Reset box
        this.selectionBox.style.display = 'none';
        this.selectionBox.style.width = '0';
        this.selectionBox.style.height = '0';
    }

    onMouseMove(e) {
        if (!this.isDragging) return;

        const dx = e.clientX - this.startPos.x;
        const dy = e.clientY - this.startPos.y;

        // If moved enough, show box (Marquee)
        if (Math.abs(dx) > 5 || Math.abs(dy) > 5) {
            if (this.selectionBox.style.display !== 'block') {
                this.selectionBox.style.display = 'block';
            }

            const x = dx < 0 ? e.clientX : this.startPos.x;
            const y = dy < 0 ? e.clientY : this.startPos.y;
            const w = Math.abs(dx);
            const h = Math.abs(dy);

            // Use transform for position to avoid layout thrashing
            // We still need to set width/height for the box size (unless we use scale, but that scales borders)
            // Setting width/height is okay if we don't read back layout properties immediately.

            this.selectionBox.style.transform = `translate(${x}px, ${y}px)`;
            this.selectionBox.style.width = `${w}px`;
            this.selectionBox.style.height = `${h}px`;
        }
    }

    onMouseUp(e) {
        // console.log(`[Editor] MouseUp: Button=${e.button} Dragging=${this.isDragging}`);
        window.lastMouseUpEvent = {
            button: e.button,
            clientX: e.clientX,
            clientY: e.clientY,
            time: Date.now(),
            wasDragging: this.isDragging
        };
        if (!this.isDragging) return;
        this.isDragging = false;
        this.selectionBox.style.display = 'none';

        const dx = e.clientX - this.startPos.x;
        const dy = e.clientY - this.startPos.y;

        // Determine Mode
        let mode = 'replace';
        if (e.shiftKey) mode = 'add';
        else if (e.ctrlKey || e.altKey) mode = 'subtract';

        if (Math.abs(dx) < 5 && Math.abs(dy) < 5) {
            // Click (Single Pick)
            // console.log("[Editor] Click detected. Calling pick()...");
            this.pick(e.clientX, e.clientY, mode);
        } else {
            // Drag (Box Select)
            // console.log("[Editor] Drag detected. Calling boxSelect()...");
            this.boxSelect(this.startPos.x, this.startPos.y, e.clientX, e.clientY, mode);
        }
    }

    pick(x, y, mode) {
        if (this.selectionLocked) return;

        // Fix: Clear gizmo state before picking
        this.initialAtomStates = null;

        // Get Ray from camera
        this.updateRay(x, y); // Updates this.rayOrigin and this.rayDir

        // DEBUG LOGGING
        let debugLog = "";
        const debug = window.VERBOSITY_LEVEL >= Logger.DEBUG;

        if (debug) {
            debugLog = "--- PICK DEBUG ---\n";
            debugLog += `Mouse: (${x}, ${y})\n`;
            debugLog += `Ray Origin: ${this.rayOrigin.x.toFixed(3)}, ${this.rayOrigin.y.toFixed(3)}, ${this.rayOrigin.z.toFixed(3)}\n`;
            debugLog += `Ray Dir:    ${this.rayDir.x.toFixed(3)}, ${this.rayDir.y.toFixed(3)}, ${this.rayDir.z.toFixed(3)}\n`;
        }

        let minT = Infinity;
        let pickedId = -1;

        // Loop through all atoms
        for (let i = 0; i < this.system.nAtoms; i++) {
            const ax = this.system.pos[i * 3];
            const ay = this.system.pos[i * 3 + 1];
            const az = this.system.pos[i * 3 + 2];
            // Reuse tempVec for center
            this.tempVec.set(ax, ay, az);

            const radius = 0.5;

            const t = this.raySphere(this.rayOrigin, this.rayDir, this.tempVec, radius);

            if (t < minT && t > 0) {
                minT = t;
                pickedId = i;
            }
        }

        if (debug) {
            debugLog += `Picked ID: ${pickedId}, minT: ${minT}\n`;
            window.logger.debug(debugLog);
        }

        if (pickedId !== -1) {
            this.system.select(pickedId, mode);
            this.molRenderer.update();
            this.updateGizmo(); // Update Gizmo Position
            if (this.onSelectionChange) this.onSelectionChange();
            window.logger.info(`Selected Atom ${pickedId} (t=${minT.toFixed(2)})`);
            return;
        }

        // If clicked nothing
        if (mode === 'replace') {
            this.system.clearSelection();
            this.molRenderer.update();
            this.updateGizmo(); // Update Gizmo Position
            if (this.onSelectionChange) this.onSelectionChange();
        }
    }

    updateRay(x, y) {
        // Normalize mouse coordinates relative to canvas
        const rect = this.renderer.domElement.getBoundingClientRect();
        const ndcX = ((x - rect.left) / rect.width) * 2 - 1;
        const ndcY = -((y - rect.top) / rect.height) * 2 + 1;

        if (this.camera.isPerspectiveCamera) {
            this.rayOrigin.setFromMatrixPosition(this.camera.matrixWorld);
            this.rayDir.set(ndcX, ndcY, 0.5).unproject(this.camera).sub(this.rayOrigin).normalize();
        } else if (this.camera.isOrthographicCamera) {
            this.rayOrigin.set(ndcX, ndcY, -1).unproject(this.camera);
            this.rayDir.set(0, 0, -1).transformDirection(this.camera.matrixWorld);
        }
    }

    // Custom Ray-Sphere Intersection (Geometric)
    // Returns distance t, or Infinity if no hit
    raySphere(rayOrigin, rayDir, center, radius) {
        // oc = center - rayOrigin
        // Use tempVec for oc? No, center is passed in.
        // We need a temp vector for calculation. Let's create one more or just use local vars if possible.
        // Vector3 operations modify in place or return new.
        // Let's use a local clone or just math.

        const cx = center.x - rayOrigin.x;
        const cy = center.y - rayOrigin.y;
        const cz = center.z - rayOrigin.z;

        const tca = cx * rayDir.x + cy * rayDir.y + cz * rayDir.z;
        if (tca < 0) return Infinity; // Behind ray

        const d2 = (cx * cx + cy * cy + cz * cz) - tca * tca;
        const r2 = radius * radius;

        if (d2 > r2) return Infinity; // Miss

        const thc = Math.sqrt(r2 - d2);
        const t0 = tca - thc;
        // const t1 = tca + thc; // Far hit

        return t0;
    }

    boxSelect(x1, y1, x2, y2, mode) {
        if (this.selectionLocked) return;

        // Fix: Clear gizmo state
        this.initialAtomStates = null;

        const minX = Math.min(x1, x2);
        const maxX = Math.max(x1, x2);
        const minY = Math.min(y1, y2);
        const maxY = Math.max(y1, y2);

        const rect = this.renderer.domElement.getBoundingClientRect();
        const width = rect.width;
        const height = rect.height;

        if (mode === 'replace') {
            this.system.clearSelection();
            // After clearing, replace behaves like add
            mode = 'add';
        }

        let count = 0;
        // Reuse tempVec

        // Optimization: Define action function outside loop or use separate loops
        // Separate loops is faster as it avoids function call overhead or condition check per atom

        if (mode === 'subtract') {
            for (let i = 0; i < this.system.nAtoms; i++) {
                this.tempVec.set(
                    this.system.pos[i * 3],
                    this.system.pos[i * 3 + 1],
                    this.system.pos[i * 3 + 2]
                );
                this.tempVec.project(this.camera);
                const sx = (this.tempVec.x * 0.5 + 0.5) * width + rect.left;
                const sy = (-(this.tempVec.y * 0.5) + 0.5) * height + rect.top;

                if (sx >= minX && sx <= maxX && sy >= minY && sy <= maxY) {
                    if (this.tempVec.z < 1.0) {
                        this.system.selection.delete(i);
                        count++;
                    }
                }
            }
        } else { // add (and replace which was converted to add)
            for (let i = 0; i < this.system.nAtoms; i++) {
                this.tempVec.set(
                    this.system.pos[i * 3],
                    this.system.pos[i * 3 + 1],
                    this.system.pos[i * 3 + 2]
                );
                this.tempVec.project(this.camera);
                const sx = (this.tempVec.x * 0.5 + 0.5) * width + rect.left;
                const sy = (-(this.tempVec.y * 0.5) + 0.5) * height + rect.top;

                if (sx >= minX && sx <= maxX && sy >= minY && sy <= maxY) {
                    if (this.tempVec.z < 1.0) {
                        this.system.selection.add(i);
                        count++;
                    }
                }
            }
        }

        this.system.isDirty = true;
        this.molRenderer.update();
        this.molRenderer.update();
        this.updateGizmo(); // Update Gizmo Position
        if (this.onSelectionChange) this.onSelectionChange();
        window.logger.info(`Box Selected ${count} atoms.`);
    }
}
