import { Draw3D } from './Draw3D.js';

export class MeshRenderer {
    constructor(scene, shaders, capacity) {
        this.scene = scene;
        this.shaders = shaders;
        this.capacity = capacity;

        this.posTexture = null;
        this.posData = null; // Float32Array for texture
        this.texSize = 64; // Power of 2
        this.texHeight = 0;

        this.atomMesh = null;
        this.bondLines = null;
        this.labelMesh = null;
        this.selectionMesh = null; // Highlight selected vertices
        this.fontTexture = null;

        if (this.scene && this.shaders) {
            this.init();
        }
    }

    init() {
        // --- 1. Position Texture ---
        this.texHeight = Math.ceil(this.capacity / this.texSize);
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

        this.commonUniforms = {
            uPosTex:   { value: this.posTexture },
            uTexSize:  { value: new THREE.Vector2(this.texSize, this.texHeight) },
            // Global point-size scale baseline; specific meshes may override
            // this with their own uniform objects so they can scale
            // independently (atoms vs selection).
            uPointScale: { value: 1.0 }
        };

        // --- 2. Atoms (InstancedMesh) ---
        if (this.shaders.atom) {
            const atomGeo = new THREE.PlaneBufferGeometry(1, 1);

            // Per-mesh uniform blocks so we can scale atoms and selection
            // independently while still sharing common texture uniforms.
            const atomUniforms = { ...this.commonUniforms, uPointScale: { value: 1.0 } };

            // Main atom mesh: uses atomGeo
            this.atomMesh = Draw3D.createTextureBasedInstancedMesh(
                this.capacity,
                atomGeo,
                this.shaders.atom,
                atomUniforms
            );

            // Initialize instance colors and scales
            const colors = new Float32Array(this.capacity * 3);
            const scales = new Float32Array(this.capacity);
            for (let i = 0; i < this.capacity; i++) {
                colors[i * 3] = 1; colors[i * 3 + 1] = 1; colors[i * 3 + 2] = 1;
                scales[i] = 0.0; // Hidden by default
            }
            this.atomMesh.geometry.setAttribute('instanceColor', new THREE.InstancedBufferAttribute(colors, 3));
            this.atomMesh.geometry.setAttribute('instanceScale', new THREE.InstancedBufferAttribute(scales, 1));

            this.scene.add(this.atomMesh);

            // Selection mesh: reuse same shaders but its OWN cloned quad geometry
            // and its own uniform block, so its uPointScale and instance
            // attributes are fully independent from the main atom mesh.
            const selUniforms = { ...this.commonUniforms, uPointScale: { value: 1.2 } };
            const selGeo = atomGeo.clone();
            this.selectionMesh = Draw3D.createTextureBasedInstancedMesh(
                this.capacity,
                selGeo,
                this.shaders.atom,
                selUniforms
            );

            const selColors = new Float32Array(this.capacity * 3);
            const selScales = new Float32Array(this.capacity);
            for (let i = 0; i < this.capacity; i++) {
                // Bright yellow for selection
                selColors[i * 3] = 1.0;
                selColors[i * 3 + 1] = 1.0;
                selColors[i * 3 + 2] = 0.0;
                selScales[i] = 0.0; // Hidden until selected
            }
            this.selectionMesh.geometry.setAttribute('instanceColor', new THREE.InstancedBufferAttribute(selColors, 3));
            this.selectionMesh.geometry.setAttribute('instanceScale', new THREE.InstancedBufferAttribute(selScales, 1));
            this.selectionMesh.renderOrder = 998; // Draw above atoms but below labels
            this.selectionMesh.visible = true;
            this.scene.add(this.selectionMesh);
        }

        // --- 3. Bonds (LineSegments) ---
        if (this.shaders.bond) {
            const maxBonds = this.capacity * 4;
            this.bondLines = Draw3D.createTextureBasedLineSegments(
                this.capacity,
                maxBonds,
                this.shaders.bond,
                {
                    ...this.commonUniforms,
                    // Palette of up to 8 material colors for bonds/edges.
                    // Higher-level code should fill this with desired colors.
                    uMatColors: { value: [
                        new THREE.Vector4(1.0, 1.0, 1.0, 1.0), // 0: default white
                        new THREE.Vector4(1.0, 0.0, 0.0, 1.0), // 1: example (rail)
                        new THREE.Vector4(0.0, 1.0, 0.0, 1.0), // 2: example (radial)
                        new THREE.Vector4(0.0, 0.0, 1.0, 1.0), // 3: example (diagonal)
                        new THREE.Vector4(1.0, 1.0, 0.0, 1.0), // 4
                        new THREE.Vector4(0.0, 1.0, 1.0, 1.0), // 5
                        new THREE.Vector4(1.0, 0.0, 1.0, 1.0), // 6
                        new THREE.Vector4(0.5, 0.5, 0.5, 1.0)  // 7
                    ] }
                }
            );
            this.scene.add(this.bondLines);
        }

        // --- 4. Labels ---
        if (this.shaders.label) {
            this.fontTexture = Draw3D.createFontTexture();
            this.labelMesh = Draw3D.createLabelInstancedMesh(
                this.capacity,
                this.shaders.label,
                this.fontTexture,
                this.commonUniforms
            );
            this.scene.add(this.labelMesh);
            this.labelMesh.visible = false;
        }
    }

    updatePositions(posArray, count) {
        // Update Texture
        // posArray is expected to be flat [x, y, z, ...]
        // count is number of items
        for (let i = 0; i < count; i++) {
            this.posData[i * 4] = posArray[i * 3];
            this.posData[i * 4 + 1] = posArray[i * 3 + 1];
            this.posData[i * 4 + 2] = posArray[i * 3 + 2];
            this.posData[i * 4 + 3] = 1.0; // w
        }
        this.posTexture.needsUpdate = true;
    }

    updateParticles(count, colorGetter, scaleGetter) {
        if (!this.atomMesh) return;

        this.atomMesh.count = count;
        const colorAttr = this.atomMesh.geometry.getAttribute('instanceColor');
        const scaleAttr = this.atomMesh.geometry.getAttribute('instanceScale');

        for (let i = 0; i < count; i++) {
            const col = colorGetter(i);
            const scale = scaleGetter(i);

            colorAttr.setXYZ(i, col[0], col[1], col[2]);
            scaleAttr.setX(i, scale);
        }
        colorAttr.needsUpdate = true;
        scaleAttr.needsUpdate = true;
    }

    /**
     * Update selection highlight based on an array of vertex indices.
     * Renders a small yellow marker per selected vertex on top of atoms.
     * @param {Array<number>} indices
     */
    updateSelection(indices) {
        if (!this.selectionMesh) return;

        const maxCount = Math.min(indices.length, this.capacity);
        const idAttr = this.selectionMesh.geometry.getAttribute('aAtomID');
        const scaleAttr = this.selectionMesh.geometry.getAttribute('instanceScale');
        if (!idAttr || !scaleAttr) return;

        const idArray = idAttr.array;
        const scaleArray = scaleAttr.array;

        // Reset all scales to 0 (hidden)
        for (let i = 0; i < this.capacity; i++) {
            scaleArray[i] = 0.0;
        }

        for (let i = 0; i < maxCount; i++) {
            const idx = indices[i];
            idArray[i] = idx;
            // Slightly larger than default sphere; actual visual size controlled by shader scale
            scaleArray[i] = 1.2;
        }

        idAttr.needsUpdate = true;
        scaleAttr.needsUpdate = true;
        this.selectionMesh.count = maxCount;
        this.selectionMesh.visible = maxCount > 0;
    }

    updateBonds(pairs) {
        if (!this.bondLines) return;

        const bondIDAttr = this.bondLines.geometry.getAttribute('aAtomID');
        const matIDAttr  = this.bondLines.geometry.getAttribute('aMatID');

        let ptr = 0;
        for (let i = 0; i < pairs.length; i++) {
            const entry = pairs[i];
            const i0 = entry[0];
            const i1 = entry[1];
            const matID = (entry.length >= 3) ? entry[2] : 0;
            //const matID = (entry.length >= 3) ? entry[2] : (i % 8); // If caller does not provide a material ID, cycle through palette by edge index.

            bondIDAttr.setX(ptr, i0);
            matIDAttr.setX(ptr, matID);
            ptr++;

            bondIDAttr.setX(ptr, i1);
            matIDAttr.setX(ptr, matID);
            ptr++;
        }

        bondIDAttr.needsUpdate = true;
        matIDAttr.needsUpdate = true;
        this.bondLines.geometry.setDrawRange(0, pairs.length * 2);
    }

    updateLabels(stringGetter, count) {
        if (!this.labelMesh || !this.labelMesh.visible) return;
        Draw3D.updateLabelBuffers(this.labelMesh, stringGetter, count);
    }

    setLabelStyle(colorHex, scale, screenSpace) {
        if (this.labelMesh) {
            if (colorHex) {
                this.labelMesh.material.uniforms.uColor.value.set(colorHex);
            }
            if (scale !== undefined && scale !== null) {
                this.labelMesh.material.uniforms.uScale.value = parseFloat(scale);
            }
            if (screenSpace !== undefined && screenSpace !== null) {
                this.labelMesh.material.uniforms.uScreenSpace.value = !!screenSpace;
            }
        }
    }

    updateLabelUniforms(aspect) {
        if (this.labelMesh && aspect) {
            this.labelMesh.material.uniforms.uAspect.value = aspect;
        }
    }

    // --- Visibility Toggles ---

    setAtomsVisible(flag) {
        if (this.atomMesh) this.atomMesh.visible = !!flag;
    }

    setBondsVisible(flag) {
        if (this.bondLines) this.bondLines.visible = !!flag;
    }

    setLabelsVisible(flag) {
        if (this.labelMesh) this.labelMesh.visible = !!flag;
    }

    setSelectionVisible(flag) {
        if (this.selectionMesh) this.selectionMesh.visible = !!flag;
    }

    // --- Size Controls via shader uniforms ---

    /**
     * Set global size multiplier for atom sprites (vertices).
     * Implemented purely as a shader uniform (uPointScale) so it does not
     * re-upload textures or instance buffers.
     */
    setNodeScale(scale) {
        if (!this.atomMesh || !this.atomMesh.material || !this.atomMesh.material.uniforms) return;
        if (this.atomMesh.material.uniforms.uPointScale) {
            this.atomMesh.material.uniforms.uPointScale.value = scale;
        }
    }

    /**
     * Set global size multiplier for selection markers (selected vertices).
     * Uses the same shader uniform name (uPointScale) but on the separate
     * selection mesh material, so selection and normal vertices can have
     * independent scales.
     */
    setSelectionScale(scale) {
        if (!this.selectionMesh || !this.selectionMesh.material || !this.selectionMesh.material.uniforms) return;
        if (this.selectionMesh.material.uniforms.uPointScale) {
            this.selectionMesh.material.uniforms.uPointScale.value = scale;
        }
    }
}

// Optional browser global for legacy code / debugging
if (typeof window !== 'undefined') {
    window.MeshRenderer = MeshRenderer;
}
