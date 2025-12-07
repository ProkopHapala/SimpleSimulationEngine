// Simple Three.js-based renderer for the MHD demo (r-z view)
import * as THREE from 'https://cdn.jsdelivr.net/npm/three@0.160.0/build/three.module.js';
import { OrbitControls } from 'https://cdn.jsdelivr.net/npm/three@0.160.0/examples/jsm/controls/OrbitControls.js';
import { sampleFieldAtPoint } from './physics.js';
import { Draw3D } from '../common_js/Draw3D.js';
import { Vec3 } from '../common_js/Vec3.js';

export class DemoRenderer {
    constructor(container) {
        this.container = container;
        const w = container.clientWidth;
        const h = container.clientHeight;
        this.scene = new THREE.Scene();
        this.scene.background = new THREE.Color(0x0d0f12);

        const aspect = w / h;
        const size = 3.0;
        this.camera = new THREE.OrthographicCamera(-size * aspect, size * aspect, size, -size, -10, 10);
        this.camera.position.set(0, 0, 5);
        this.camera.lookAt(0, 0, 0);

        this.renderer = new THREE.WebGLRenderer({ antialias: true });
        this.renderer.setSize(w, h);
        container.appendChild(this.renderer.domElement);

        this.controls = new OrbitControls(this.camera, this.renderer.domElement);
        this.controls.enableRotate = true;
        this.controls.enablePan = true;
        this.controls.enableZoom = true;
        this.controls.enableDamping = true;

        // Use LineSegments for coils to allow independent loops
        this.plasmaLine = this.makeSegments(0x00c2ff);
        this.cageLine = this.makeSegments(0xffa500);
        this.scLine = this.makeSegments(0x66ff66);
        this.scene.add(this.plasmaLine);
        this.scene.add(this.cageLine);
        this.scene.add(this.scLine);

        this.fieldPoints = this.makeFieldPoints();
        this.scene.add(this.fieldPoints);

        this.tempVec = new Vec3(); // reuse

        window.addEventListener('resize', () => this.onResize());
    }

    makeSegments(color) {
        const geo = new THREE.BufferGeometry();
        // Initial buffer size
        const maxVerts = 2000;
        geo.setAttribute('position', new THREE.BufferAttribute(new Float32Array(maxVerts * 3), 3));
        const mat = new THREE.LineBasicMaterial({ color });
        const mesh = new THREE.LineSegments(geo, mat);
        mesh.frustumCulled = false;
        return mesh;
    }

    makeFieldPoints() {
        const n = 32;
        const positions = new Float32Array(n * n * 3);
        const colors = new Float32Array(n * n * 3);
        const geo = new THREE.BufferGeometry();
        geo.setAttribute('position', new THREE.BufferAttribute(positions, 3));
        geo.setAttribute('color', new THREE.BufferAttribute(colors, 3));
        const mat = new THREE.PointsMaterial({ size: 0.04, vertexColors: true, transparent: true, opacity: 0.9 });
        return new THREE.Points(geo, mat);
    }

    updateField(sim) {
        const n = 32;
        const sizeR = 1.8;
        const sizeZ = 1.8;
        const pos = this.fieldPoints.geometry.getAttribute('position').array;
        const col = this.fieldPoints.geometry.getAttribute('color').array;
        let idx = 0;
        let maxB = 1e-5; // sensible min
        const mags = [];
        for (let j = 0; j < n; j++) {
            const z = -sizeZ + (2 * sizeZ * j) / (n - 1);
            for (let i = 0; i < n; i++) {
                const r = (2 * sizeR * i) / (n - 1);
                pos[idx * 3 + 0] = z; // Map Z to X for viewing
                pos[idx * 3 + 1] = r; // Map R to Y
                pos[idx * 3 + 2] = 0;
                const B = sampleFieldAtPoint(Math.abs(r), z, sim);
                let m = Math.sqrt(B.Br * B.Br + B.Bz * B.Bz);
                if (!isFinite(m)) m = 0;
                mags[idx] = m;
                if (m > maxB) maxB = m;
                idx++;
            }
        }
        idx = 0;
        for (let k = 0; k < mags.length; k++) {
            const v = Math.min(1.0, mags[k] / (maxB + 1e-9));
            const rcol = v;
            const gcol = Math.max(0, 1.0 - Math.abs(v - 0.5) * 2.0);
            const bcol = 1.0 - v;
            col[idx * 3 + 0] = rcol;
            col[idx * 3 + 1] = gcol;
            col[idx * 3 + 2] = bcol;
            idx++;
        }
        this.fieldPoints.geometry.attributes.position.needsUpdate = true;
        this.fieldPoints.geometry.attributes.color.needsUpdate = true;
    }

    updateLines(sim) {
        const updateMesh = (mesh, coilList, isPlasma) => {
            const vertices = [];
            const pos = { x: 0, y: 0, z: 0 }; // center of coil
            const axis = { x: 1, y: 0, z: 0 }; // axis of coil (Z-axis in physics maps to X-axis in view?)
            // Wait, coordinate mapping: R->Y, Z->X. 
            // Physics: Z is symmetry axis. R is radial.
            // Screen: X is Horizontal, Y is Vertical.
            // Typical 2D plot: Z on X-axis, R on Y-axis.
            // In 3D view: 
            // We want to see the rings.
            // If Z is X-axis, then the rings are in Y-Z plane (of screen coords)?
            // No, the rings are around the Z-axis (physics).
            // So they form circles in the X-Y plane (physics).
            // But we display Z (physics) as X (screen).
            // This is getting confusing.
            // Let's stick to: Physics Z -> Screen X. Physics R -> Screen Y (radius).
            // So the ring is "standing up" perpendicular to Screen X?
            // Yes. The ring is in the plane Z_phys = constant.
            // So it lives in (X_phys, Y_phys).
            // We only see R_phys = sqrt(X^2 + Y^2).
            // If we rotate the view in 3D, we want to see the actual 3D rings.
            // Physics Z -> 3D World X.
            // Physics R -> radial distance from 3D World X.
            // Start vector can be (0, R, 0) (which is Y in 3D world).
            // Axis of rotation is (1, 0, 0) (X in 3D world).

            const nSegs = 32;

            // Check mapping from updateField:
            // pos[idx*3+0] = z; // X
            // pos[idx*3+1] = r; // Y
            // So Z_phys -> X_world.
            // R_phys -> Y_world (at phi=0).

            const uaxis = { x: 1, y: 0, z: 0 };

            if (isPlasma) {
                // Plasma is a list of nodes forming a loop in R-Z plane.
                // It is actually a "Surface of Revolution".
                // Detailed 3D view would be a mesh.
                // But for "Simple Lines", maybe we just draw the R-Z profile?
                // User asked to see "Expanding plasma coils".
                // If we draw rings for every node, it might be cluttered.
                // But let's try drawing rings for every node first.
                coilList.forEach(p => {
                    const startV = { x: 0, y: p.r, z: 0 };
                    const center = { x: p.z, y: 0, z: 0 };
                    Draw3D.drawCircleAxis(nSegs, center, startV, uaxis, p.r, vertices);
                });
                // Also draw the profile?
                // Connect centers (p.z, p.r)
                /*
                for(let i=0; i<coilList.length; i++) {
                    const p0 = coilList[i];
                    const p1 = coilList[(i+1)%coilList.length];
                    vertices.push(p0.z, p0.r, 0);
                    vertices.push(p1.z, p1.r, 0);
                    // mirror ?
                    vertices.push(p0.z, -p0.r, 0);
                    vertices.push(p1.z, -p1.r, 0);
                }
                */
            } else {
                coilList.forEach(c => {
                    const startV = { x: 0, y: c.r, z: 0 };
                    const center = { x: c.z, y: 0, z: 0 };
                    Draw3D.drawCircleAxis(nSegs, center, startV, uaxis, c.r, vertices);
                });
            }

            // Expand buffer if needed
            const needed = vertices.length;
            const currentSize = mesh.geometry.attributes.position.array.length;
            if (needed > currentSize) {
                mesh.geometry.setAttribute('position', new THREE.BufferAttribute(new Float32Array(needed * 2), 3));
            }

            const arr = mesh.geometry.attributes.position.array;
            for (let i = 0; i < needed; i++) arr[i] = vertices[i];

            mesh.geometry.setDrawRange(0, needed / 3);
            mesh.geometry.attributes.position.needsUpdate = true;
            mesh.geometry.computeBoundingSphere();
        };

        updateMesh(this.plasmaLine, sim.plasmaNodes, true);
        updateMesh(this.cageLine, sim.cageCoils, false);
        updateMesh(this.scLine, sim.scCoils, false);
    }

    onResize() {
        const w = this.container.clientWidth;
        const h = this.container.clientHeight;
        const aspect = w / h;
        const size = 3.0;
        this.camera.left = -size * aspect;
        this.camera.right = size * aspect;
        this.camera.top = size;
        this.camera.bottom = -size;
        this.camera.updateProjectionMatrix();
        this.renderer.setSize(w, h);
    }

    render(sim) {
        this.updateLines(sim);
        this.updateField(sim);
        this.renderer.render(this.scene, this.camera);
    }
}
