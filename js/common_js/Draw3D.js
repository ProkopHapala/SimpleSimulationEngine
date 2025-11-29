export class Draw3D {
    static drawSphereOctLines(n, R, pos, vertices) {


        // Axes
        const a = { x: 1, y: 0, z: 0 };
        const b = { x: 0, y: 1, z: 0 };
        const c = { x: 0, y: 0, z: 1 };

        // Draw 3 orthogonal circles
        // 1. XY plane (rot around Z (c)) -> Start at X (a)
        this.drawCircleAxis(n, pos, a, c, R, vertices);
        // 2. YZ plane (rot around X (a)) -> Start at Y (b)
        this.drawCircleAxis(n, pos, b, a, R, vertices);
        // 3. XZ plane (rot around Y (b)) -> Start at Z (c)
        this.drawCircleAxis(n, pos, c, b, R, vertices);
    }

    // Optimized circle drawing using recurrence (no sin/cos in loop)
    // v0: start vector (relative to pos)
    // uaxis: rotation axis (normalized)
    static drawCircleAxis(n, pos, v0, uaxis, R, vertices) {

        const dphi = 2 * Math.PI / n;
        const dca = Math.cos(dphi);
        const dsa = Math.sin(dphi);

        let x = v0.x;
        let y = v0.y;
        let z = v0.z;

        // Current point (start)
        let px = pos.x + x * R;
        let py = pos.y + y * R;
        let pz = pos.z + z * R;

        for (let i = 0; i < n; i++) {
            // Calculate next vector using Rodrigues' rotation (simplified for incremental)
            // v_next = v * cos + (axis x v) * sin
            // axis x v:
            const cx = uaxis.y * z - uaxis.z * y;
            const cy = uaxis.z * x - uaxis.x * z;
            const cz = uaxis.x * y - uaxis.y * x;

            const nx = x * dca + cx * dsa;
            const ny = y * dca + cy * dsa;
            const nz = z * dca + cz * dsa;

            // Next point
            const p2x = pos.x + nx * R;
            const p2y = pos.y + ny * R;
            const p2z = pos.z + nz * R;

            // Add line segment
            vertices.push(px, py, pz);
            vertices.push(p2x, p2y, p2z);

            // Update for next iteration
            x = nx;
            y = ny;
            z = nz;
            px = p2x;
            py = p2y;
            pz = p2z;
        }
    }

    static createOctSphereGeometry(n, R) {
        const vertices = [];
        this.drawSphereOctLines(n, R, { x: 0, y: 0, z: 0 }, vertices);
        const geometry = new THREE.BufferGeometry();
        geometry.setAttribute('position', new THREE.Float32BufferAttribute(vertices, 3));
        return geometry;
    }

    static updateMergedLineGeometry(geometry, baseGeometry, indices, positionArray) {
        const baseVertices = baseGeometry.attributes.position.array;
        const vertexCountPerInstance = baseVertices.length / 3;
        const instanceCount = indices.length;
        const totalVertexCount = instanceCount * vertexCountPerInstance;

        // Resize buffer if needed
        if (!geometry.attributes.position || geometry.attributes.position.array.length < totalVertexCount * 3) {
            // Allocate with some headroom or exact? Let's do exact for now or 2x?
            // Since we want to avoid reallocation, maybe allocate for capacity if passed?
            // For now, just allocate what's needed if too small.
            geometry.setAttribute('position', new THREE.BufferAttribute(new Float32Array(totalVertexCount * 3), 3));
        }

        const targetArray = geometry.attributes.position.array;

        for (let i = 0; i < instanceCount; i++) {
            const id = indices[i];
            const cx = positionArray[id * 3];
            const cy = positionArray[id * 3 + 1];
            const cz = positionArray[id * 3 + 2];

            const offset = i * vertexCountPerInstance * 3;

            for (let j = 0; j < vertexCountPerInstance; j++) {
                targetArray[offset + j * 3] = baseVertices[j * 3] + cx;
                targetArray[offset + j * 3 + 1] = baseVertices[j * 3 + 1] + cy;
                targetArray[offset + j * 3 + 2] = baseVertices[j * 3 + 2] + cz;
            }
        }

        geometry.attributes.position.needsUpdate = true;
        geometry.setDrawRange(0, totalVertexCount);
    }

    // --- Generic Texture-Based Rendering Helpers ---

    static createTextureBasedInstancedMesh(capacity, geometry, shaders, uniforms) {
        const material = new THREE.ShaderMaterial({
            vertexShader: shaders.vertex,
            fragmentShader: shaders.fragment,
            uniforms: uniforms,
            side: THREE.DoubleSide
        });

        const mesh = new THREE.InstancedMesh(geometry, material, capacity);
        mesh.instanceMatrix.setUsage(THREE.DynamicDrawUsage);

        // Default ID attribute (0..N)
        const ids = new Float32Array(capacity);
        for (let i = 0; i < capacity; i++) ids[i] = i;
        mesh.geometry.setAttribute('aAtomID', new THREE.InstancedBufferAttribute(ids, 1));
        mesh.count = 0; // Start with 0 visible

        return mesh;
    }

    static createTextureBasedLineSegments(capacity, maxSegments, shaders, uniforms) {
        const material = new THREE.ShaderMaterial({
            vertexShader: shaders.vertex,
            fragmentShader: shaders.fragment,
            uniforms: uniforms,
            depthTest: true,
            depthWrite: true
        });

        const geometry = new THREE.BufferGeometry();
        // Dummy position buffer (needed for frustum culling or just ignore)
        geometry.setAttribute('position', new THREE.BufferAttribute(new Float32Array(maxSegments * 2 * 3), 3));

        // Atom ID attribute (index into position texture)
        const ids = new Float32Array(maxSegments * 2);
        geometry.setAttribute('aAtomID', new THREE.BufferAttribute(ids, 1));

        // Material ID attribute (index into uMatColors palette)
        const matIDs = new Float32Array(maxSegments * 2);
        geometry.setAttribute('aMatID', new THREE.BufferAttribute(matIDs, 1));

        const lines = new THREE.LineSegments(geometry, material);
        lines.frustumCulled = false;
        lines.geometry.setDrawRange(0, 0); // Start with 0 visible
        return lines;
    }

    static createTextureBasedSelectionLines(capacity, baseGeometry, shaders, uniforms) {
        const material = new THREE.ShaderMaterial({
            vertexShader: shaders.vertex,
            fragmentShader: shaders.fragment,
            uniforms: uniforms,
            depthTest: false,
            depthWrite: false,
            transparent: true
        });

        const geometry = new THREE.InstancedBufferGeometry();
        geometry.setAttribute('position', baseGeometry.getAttribute('position'));

        // ID Attribute
        const ids = new Float32Array(capacity);
        geometry.setAttribute('aAtomID', new THREE.InstancedBufferAttribute(ids, 1));

        const lines = new THREE.LineSegments(geometry, material);
        lines.frustumCulled = false;
        lines.geometry.instanceCount = 0; // Start with 0 visible
        return lines;
    }

    // --- Generic Label Rendering ---

    static createFontTexture() {
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

        const texture = new THREE.CanvasTexture(canvas);
        texture.minFilter = THREE.LinearFilter;
        texture.magFilter = THREE.LinearFilter;
        texture.needsUpdate = true;
        return texture;
    }

    static createLabelInstancedMesh(capacity, shaders, fontTexture, uniforms) {
        // Geometry: 8 quads (max 8 chars per label instance).
        // We build a base geometry with 8 copies of a unit quad and use instancing
        // to draw one label per atom. The actual characters are provided via
        // instanced attributes (aLabel1/aLabel2/aStrLen). Per-vertex attribute
        // aCharPos (0..7) selects which packed character this quad should show.
        const maxChars = 8;
        const baseGeo = new THREE.PlaneBufferGeometry(1, 1);
        const basePos = baseGeo.attributes.position.array;
        const baseUv = baseGeo.attributes.uv.array;
        const baseIndex = baseGeo.index.array;

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

        // Instanced Attributes (per-label / per-atom):
        //  - aAtomID  : which atom this label is anchored to
        //  - aLabel1  : ASCII codes for characters [0..3]
        //  - aLabel2  : ASCII codes for characters [4..7]
        //  - aStrLen  : actual length of string (<= 8), used for centering
        const instID = new Float32Array(capacity);
        for (let i = 0; i < capacity; i++) instID[i] = i;

        const instLabel1 = new Float32Array(capacity * 4);
        const instLabel2 = new Float32Array(capacity * 4);
        const instStrLen = new Float32Array(capacity);

        const instGeo = new THREE.InstancedBufferGeometry();
        instGeo.copy(geometry);
        instGeo.setAttribute('aAtomID', new THREE.InstancedBufferAttribute(instID, 1));
        instGeo.setAttribute('aLabel1', new THREE.InstancedBufferAttribute(instLabel1, 4));
        instGeo.setAttribute('aLabel2', new THREE.InstancedBufferAttribute(instLabel2, 4));
        instGeo.setAttribute('aStrLen', new THREE.InstancedBufferAttribute(instStrLen, 1));

        const finalUniforms = {
            ...uniforms,
            uFontTex: { value: fontTexture },
            uFontGrid: { value: new THREE.Vector2(16, 16) },
            uScale: { value: 0.5 },
            uColor: { value: new THREE.Color(0xffffff) },
            uScreenSpace: { value: false },
            uAspect: { value: 1.0 }
        };

        const material = new THREE.ShaderMaterial({
            vertexShader: shaders.vertex,
            fragmentShader: shaders.fragment,
            uniforms: finalUniforms,
            transparent: true,
            depthTest: false, // Always on top
            depthWrite: false,
            side: THREE.DoubleSide
        });

        const mesh = new THREE.InstancedMesh(instGeo, material, capacity);
        mesh.frustumCulled = false;
        mesh.renderOrder = 999; // Render last
        mesh.count = 0; // Start with 0 visible
        return mesh;
    }

    static updateLabelBuffers(mesh, stringGetter, count) {
        mesh.count = count;

        // This function fills the per-instance label attributes used by the
        // label shaders. For each instance i we take a short string and:
        //  - write its length into aStrLen[i]
        //  - pack up to 8 ASCII codes into two vec4s:
        //      aLabel1 = chars 0..3, aLabel2 = chars 4..7
        // The vertex shader receives aCharPos=0..7 per quad and picks the
        // corresponding character code from aLabel1/aLabel2.
        const attr1 = mesh.geometry.getAttribute('aLabel1');
        const attr2 = mesh.geometry.getAttribute('aLabel2');
        const attrLen = mesh.geometry.getAttribute('aStrLen');
        const arr1 = attr1.array;
        const arr2 = attr2.array;
        const arrLen = attrLen.array;

        for (let i = 0; i < count; i++) {
            const str = stringGetter(i) || "";
            const len = Math.min(str.length, 8);
            arrLen[i] = len;

            // Clear
            for (let k = 0; k < 8; k++) {
                if (k < 4) arr1[i * 4 + k] = 0;
                else arr2[i * 4 + (k - 4)] = 0;
            }

            for (let k = 0; k < len; k++) {
                const code = str.charCodeAt(k);
                if (k < 4) arr1[i * 4 + k] = code;
                else arr2[i * 4 + (k - 4)] = code;
            }
        }

        attr1.needsUpdate = true;
        attr2.needsUpdate = true;
        attrLen.needsUpdate = true;
    }

}

// Optional browser global for legacy code / debugging
if (typeof window !== 'undefined') {
    window.Draw3D = Draw3D;
}
