class LabelRenderer {
    constructor(scene, system, shaders) {
        this.scene = scene;
        this.system = system;
        this.shaders = shaders;

        this.mesh = null;
        this.fontTexture = null;
        this.capacity = system.capacity;

        this.mode = 'none'; // none, id, element, type

        this.init();
    }

    init() {
        this.createFontTexture();
        this.createMesh();
    }

    createFontTexture() {
        const canvas = document.createElement('canvas');
        const size = 512;
        canvas.width = size;
        canvas.height = size;
        const ctx = canvas.getContext('2d');

        // Grid: 16x16 = 256 chars. ASCII is 128. 
        // We can use 16 cols, 8 rows for 128 chars.
        const cols = 16;
        const rows = 16;
        const charW = size / cols;
        const charH = size / rows;

        ctx.fillStyle = '#00000000'; // Transparent
        ctx.fillRect(0, 0, size, size);

        ctx.font = 'bold 24px monospace';
        ctx.fillStyle = '#FFFFFF';
        ctx.textAlign = 'center';
        ctx.textBaseline = 'middle';

        for (let i = 32; i < 127; i++) {
            const col = i % cols;
            const row = Math.floor(i / cols);

            const x = col * charW + charW / 2;
            const y = row * charH + charH / 2;

            // Adjust y for baseline
            ctx.fillText(String.fromCharCode(i), x, y + 2);
        }

        this.fontTexture = new THREE.CanvasTexture(canvas);
        this.fontTexture.minFilter = THREE.LinearFilter;
        this.fontTexture.magFilter = THREE.LinearFilter;
        this.fontTexture.needsUpdate = true;
    }

    createMesh() {
        // Geometry: 8 quads
        // We construct a BufferGeometry with 8 * 4 vertices (indexed) or just 8 quads.
        // Let's use indexed.

        const maxChars = 8;
        const baseGeo = new THREE.PlaneBufferGeometry(1, 1);
        const basePos = baseGeo.attributes.position.array;
        const baseUv = baseGeo.attributes.uv.array;
        const baseIndex = baseGeo.index.array;

        // We need to replicate this 8 times
        const vertices = [];
        const uvs = [];
        const indices = [];
        const charPosAttr = [];

        for (let c = 0; c < maxChars; c++) {
            const vOffset = (vertices.length / 3);

            for (let i = 0; i < basePos.length; i += 3) {
                vertices.push(basePos[i], basePos[i + 1], basePos[i + 2]);
                charPosAttr.push(c);
            }

            for (let i = 0; i < baseUv.length; i += 2) {
                uvs.push(baseUv[i], baseUv[i + 1]);
            }

            for (let i = 0; i < baseIndex.length; i++) {
                indices.push(baseIndex[i] + vOffset);
            }
        }

        const geometry = new THREE.BufferGeometry();
        geometry.setAttribute('position', new THREE.Float32BufferAttribute(vertices, 3));
        geometry.setAttribute('uv', new THREE.Float32BufferAttribute(uvs, 2));
        geometry.setAttribute('aCharPos', new THREE.Float32BufferAttribute(charPosAttr, 1));
        geometry.setIndex(indices);

        // Instanced Attributes
        const instID = new Float32Array(this.capacity);
        for (let i = 0; i < this.capacity; i++) instID[i] = i;

        const instLabel1 = new Float32Array(this.capacity * 4);
        const instLabel2 = new Float32Array(this.capacity * 4);
        const instStrLen = new Float32Array(this.capacity);

        const instGeo = new THREE.InstancedBufferGeometry();
        instGeo.copy(geometry);
        instGeo.setAttribute('aAtomID', new THREE.InstancedBufferAttribute(instID, 1));
        instGeo.setAttribute('aLabel1', new THREE.InstancedBufferAttribute(instLabel1, 4));
        instGeo.setAttribute('aLabel2', new THREE.InstancedBufferAttribute(instLabel2, 4));
        instGeo.setAttribute('aStrLen', new THREE.InstancedBufferAttribute(instStrLen, 1));

        // Uniforms
        const uniforms = {
            uPosTex: { value: null }, // Set by MoleculeRenderer
            uTexSize: { value: new THREE.Vector2(1, 1) },
            uFontTex: { value: this.fontTexture },
            uFontTex: { value: this.fontTexture },
            uFontGrid: { value: new THREE.Vector2(16, 16) },
            uScale: { value: 0.5 },
            uColor: { value: new THREE.Color(0xffffff) } // Default White
        };

        const material = new THREE.ShaderMaterial({
            vertexShader: this.shaders.vertex,
            fragmentShader: this.shaders.fragment,
            uniforms: uniforms,
            transparent: true,
            depthTest: false, // Always on top
            depthWrite: false,
            side: THREE.DoubleSide
        });

        this.mesh = new THREE.InstancedMesh(instGeo, material, this.capacity);
        this.mesh.frustumCulled = false;
        this.mesh.renderOrder = 999; // Render last
        this.scene.add(this.mesh);
        this.mesh.visible = false;
    }

    update(posTexture, texSize) {
        if (this.mode === 'none') {
            this.mesh.visible = false;
            return;
        }
        this.mesh.visible = true;

        // Update Uniforms
        this.mesh.material.uniforms.uPosTex.value = posTexture;
        this.mesh.material.uniforms.uTexSize.value.copy(texSize);

        // Update Labels if dirty (we assume they are dirty if update is called with mode change logic)
        // But we should only update buffers if needed.
        // Let's assume we update buffers here.

        this.updateLabels();
    }

    updateLabels() {
        const nAtoms = this.system.nAtoms;
        this.mesh.count = nAtoms;

        const attr1 = this.mesh.geometry.getAttribute('aLabel1');
        const attr2 = this.mesh.geometry.getAttribute('aLabel2');
        const attrLen = this.mesh.geometry.getAttribute('aStrLen');
        const arr1 = attr1.array;
        const arr2 = attr2.array;
        const arrLen = attrLen.array;

        // Helper to encode string
        const encode = (str, idx) => {
            const len = Math.min(str.length, 8);
            arrLen[idx] = len;

            // Clear
            for (let k = 0; k < 8; k++) {
                if (k < 4) arr1[idx * 4 + k] = 0;
                else arr2[idx * 4 + (k - 4)] = 0;
            }

            for (let k = 0; k < len; k++) {
                const code = str.charCodeAt(k);
                if (k < 4) arr1[idx * 4 + k] = code;
                else arr2[idx * 4 + (k - 4)] = code;
            }
        };

        for (let i = 0; i < nAtoms; i++) {
            let str = "";
            if (this.mode === 'id') {
                str = i.toString();
            } else if (this.mode === 'element') {
                const type = this.system.types[i];
                // Need to lookup element name from MMParams
                if (window.app.mmParams && window.app.mmParams.byAtomicNumber[type]) {
                    str = window.app.mmParams.byAtomicNumber[type].name;
                } else {
                    str = type.toString();
                }
            } else if (this.mode === 'type') {
                // Atom Type Name
                // We don't store atom type name in system.types (it stores atomic number).
                // Wait, system.types is Uint8Array. It stores atomic number?
                // Let's check MoleculeSystem.js. 
                // "this.types = new Uint8Array(this.capacity); // element type (atomic number)"
                // So we don't have specific atom type (like "C_3") stored per atom yet?
                // The user said: "atom type name (i.e. subtype of given element)"
                // If we don't store it, we can't show it.
                // For now, let's show Element for 'type' as well or "N/A".
                // Or maybe the user implies we should have it. 
                // But MoleculeSystem only has `types` which is atomic number.
                // Ah, maybe `types` IS the atom type index?
                // "this.types = new Uint8Array(this.capacity); // element type (atomic number)"
                // The comment says atomic number.

                // Let's assume for now we only have Element.
                // If the user wants Atom Type, we might need to add a field to MoleculeSystem.
                // But for this task, I'll just show Element for now or "?"

                str = "?";
            }

            encode(str, i);
        }

        attr1.needsUpdate = true;
        attr2.needsUpdate = true;
        attrLen.needsUpdate = true;
    }

    setColor(hex) {
        if (this.mesh) {
            this.mesh.material.uniforms.uColor.value.set(hex);
        }
    }

    setMode(mode) {
        if (this.mode !== mode) {
            this.mode = mode;
            // Trigger update
            if (window.app && window.app.molRenderer) {
                window.app.molRenderer.update(); // This will call our update
            }
        }
    }
}
