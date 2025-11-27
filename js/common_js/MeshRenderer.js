class MeshRenderer {
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
            uPosTex: { value: this.posTexture },
            uTexSize: { value: new THREE.Vector2(this.texSize, this.texHeight) }
        };

        // --- 2. Atoms (InstancedMesh) ---
        if (this.shaders.atom) {
            const atomGeo = new THREE.PlaneBufferGeometry(1, 1);
            this.atomMesh = Draw3D.createTextureBasedInstancedMesh(
                this.capacity,
                atomGeo,
                this.shaders.atom,
                this.commonUniforms
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
}

window.MeshRenderer = MeshRenderer;
