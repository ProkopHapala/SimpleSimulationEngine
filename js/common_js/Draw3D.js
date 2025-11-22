class Draw3D {
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

        // ID Attribute
        const ids = new Float32Array(maxSegments * 2);
        geometry.setAttribute('aAtomID', new THREE.BufferAttribute(ids, 1));

        const lines = new THREE.LineSegments(geometry, material);
        lines.frustumCulled = false;
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
        return lines;
    }
}

window.Draw3D = Draw3D;
