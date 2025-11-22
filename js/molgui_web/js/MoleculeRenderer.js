class MoleculeRenderer {
    constructor(scene, system, shaders, mmParams) {
        this.scene = scene;
        this.system = system;
        this.shaders = shaders;
        this.mmParams = mmParams;

        this.atomMesh = null;
        this.bondLines = null;
        this.selectionMesh = null;
        this.axesHelper = null;

        this.posTexture = null;
        this.posData = null; // Float32Array for texture
        this.texSize = 64; // Power of 2
        this.texHeight = 0;

        this.labelRenderer = null;

        this.init();
    }

    init() {
        // --- 1. Position Texture ---
        // Allocate texture large enough for max capacity
        // We use a square texture or just Nx1?
        // Max capacity is usually e.g. 10000. 1024x1024 is 1M atoms. Plenty.
        // Let's use 1024 width.
        this.texHeight = Math.ceil(this.system.capacity / this.texSize);
        const size = this.texSize * this.texHeight * 4; // RGBA float
        this.posData = new Float32Array(size);

        this.posTexture = new THREE.DataTexture(
            this.posData,
            this.texSize,
            this.texHeight,
            THREE.RGBAFormat,
            THREE.FloatType
        );
        this.posTexture.minFilter = THREE.NearestFilter;
        this.posTexture.magFilter = THREE.NearestFilter;
        this.posTexture.generateMipmaps = false;
        this.posTexture.needsUpdate = true;

        const commonUniforms = {
            uPosTex: { value: this.posTexture },
            uTexSize: { value: new THREE.Vector2(this.texSize, this.texHeight) }
        };

        // --- 2. Atoms (InstancedMesh) ---
        // Geometry: A simple plane (quad)
        const atomGeo = new THREE.PlaneBufferGeometry(1, 1);
        this.atomMesh = Draw3D.createTextureBasedInstancedMesh(
            this.system.capacity,
            atomGeo,
            this.shaders.atom,
            commonUniforms
        );

        // Initialize instance colors and scales
        const colors = new Float32Array(this.system.capacity * 3);
        const scales = new Float32Array(this.system.capacity);
        for (let i = 0; i < this.system.capacity; i++) {
            colors[i * 3] = 1; colors[i * 3 + 1] = 1; colors[i * 3 + 2] = 1;
            scales[i] = 0.0; // Hidden by default
        }
        this.atomMesh.geometry.setAttribute('instanceColor', new THREE.InstancedBufferAttribute(colors, 3));
        this.atomMesh.geometry.setAttribute('instanceScale', new THREE.InstancedBufferAttribute(scales, 1));

        this.scene.add(this.atomMesh);

        // --- 3. Bonds (LineSegments) ---
        const maxBonds = this.system.capacity * 4;
        this.bondLines = Draw3D.createTextureBasedLineSegments(
            this.system.capacity,
            maxBonds,
            this.shaders.bond,
            {
                ...commonUniforms,
                uColor: { value: new THREE.Vector4(1.0, 1.0, 1.0, 1.0) } // White
            }
        );
        this.scene.add(this.bondLines);

        // --- 4. Selection (InstancedMesh - Rings) ---
        this.selectionBaseGeo = Draw3D.createOctSphereGeometry(16, 1.0);
        this.selectionMesh = Draw3D.createTextureBasedSelectionLines(
            this.system.capacity,
            this.selectionBaseGeo,
            this.shaders.selection,
            {
                ...commonUniforms,
                uColor: { value: new THREE.Vector4(1.0, 1.0, 0.0, 0.8) } // Yellow
            }
        );
        this.scene.add(this.selectionMesh);

        // --- 5. Labels ---
        // Load shaders for labels
        // We assume shaders are passed in 'shaders' object. 
        // But we need to make sure main.js loads them.
        // For now, let's assume they are in shaders.label
        if (this.shaders.label) {
            this.labelRenderer = new LabelRenderer(this.scene, this.system, this.shaders.label);
        } else {
            console.warn("Label shaders not found!");
        }
    }

    update() {
        // Main update loop
        // If system is dirty, we assume topology or colors changed -> Full Update
        // If only positions moved (handled by Gizmo), we call updatePositions() explicitly
        // But for safety, if system.isDirty is true, we do everything.

        if (this.system.isDirty) {
            this.updateStructure();
            this.updatePositions(); // Structure change implies position update too
            this.system.isDirty = false;
        }

        // Always update labels (for camera billboard)
        if (this.labelRenderer) {
            this.labelRenderer.update(this.posTexture, new THREE.Vector2(this.texSize, this.texHeight));
        }
    }

    updatePositions() {
        const nAtoms = this.system.nAtoms;
        // Update Texture
        for (let i = 0; i < nAtoms; i++) {
            this.posData[i * 4] = this.system.pos[i * 3];
            this.posData[i * 4 + 1] = this.system.pos[i * 3 + 1];
            this.posData[i * 4 + 2] = this.system.pos[i * 3 + 2];
            this.posData[i * 4 + 3] = 1.0; // w
        }
        this.posTexture.needsUpdate = true;

        // Selection might need update if it depends on positions? 
        // No, selection depends on IDs which index into texture. 
        // So if texture updates, selection updates visually.
    }

    updateStructure() {
        const nAtoms = this.system.nAtoms;
        this.atomMesh.count = nAtoms;

        const colorAttr = this.atomMesh.geometry.getAttribute('instanceColor');
        const scaleAttr = this.atomMesh.geometry.getAttribute('instanceScale');

        for (let i = 0; i < nAtoms; i++) {
            // Color & Scale
            const type = this.system.types[i];

            // Use MMParams
            let col = [1, 0, 1]; // Default magenta
            let radius = 1.0;

            if (this.mmParams) {
                col = this.mmParams.getColor(type);
                // Scale: RvdW is usually around 1.5-2.0. 
                // For ball-and-stick, we usually want smaller, e.g. 0.25 * RvdW or fixed fraction.
                // Let's use a factor. The user didn't specify, but typical ball-and-stick is ~0.3-0.5 radius.
                // Let's try 0.4 * RvdW.
                radius = this.mmParams.getRadius(type) * 0.4;
            }

            colorAttr.setXYZ(i, col[0], col[1], col[2]);
            scaleAttr.setX(i, radius);
        }
        colorAttr.needsUpdate = true;
        scaleAttr.needsUpdate = true;

        // Update Bonds
        const bonds = this.system.bonds;
        const bondIDAttr = this.bondLines.geometry.getAttribute('aAtomID');
        let ptr = 0;
        for (let i = 0; i < bonds.length; i++) {
            bondIDAttr.setX(ptr++, bonds[i][0]);
            bondIDAttr.setX(ptr++, bonds[i][1]);
        }
        bondIDAttr.needsUpdate = true;
        this.bondLines.geometry.setDrawRange(0, bonds.length * 2);

        this.updateSelection();

        if (window.VERBOSITY_LEVEL >= Logger.INFO) {
            // window.logger.info(`Renderer structure updated: ${ nAtoms } atoms.`);
        }
    }

    updateSelection() {
        const selectedIDs = Array.from(this.system.selection);
        const count = selectedIDs.length;

        // Update instance count for InstancedBufferGeometry
        this.selectionMesh.geometry.instanceCount = count;

        const idAttr = this.selectionMesh.geometry.getAttribute('aAtomID');
        for (let i = 0; i < count; i++) {
            idAttr.setX(i, selectedIDs[i]);
        }
        idAttr.needsUpdate = true;
    }

    toggleAxes(visible) {
        if (!this.axesHelper) {
            this.axesHelper = new THREE.AxesHelper(5); // 5 units length
            this.scene.add(this.axesHelper);
        }
        this.axesHelper.visible = visible;
    }

    setLabelMode(mode) {
        if (this.labelRenderer) {
            this.labelRenderer.setMode(mode);
            this.update(); // Trigger re-render
        }
    }
}

