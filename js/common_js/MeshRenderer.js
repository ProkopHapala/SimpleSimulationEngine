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
                    uColor: { value: new THREE.Vector4(1.0, 1.0, 1.0, 1.0) } // White default
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
        let ptr = 0;
        for (let i = 0; i < pairs.length; i++) {
            bondIDAttr.setX(ptr++, pairs[i][0]);
            bondIDAttr.setX(ptr++, pairs[i][1]);
        }
        bondIDAttr.needsUpdate = true;
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
