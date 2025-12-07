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
        // Use LineSegments for vector field
        // Each vector is a small line. 
        // 2 vertices per point.
        const positions = new Float32Array(n * n * 2 * 3);
        const colors = new Float32Array(n * n * 2 * 3);
        const geo = new THREE.BufferGeometry();
        geo.setAttribute('position', new THREE.BufferAttribute(positions, 3));
        geo.setAttribute('color', new THREE.BufferAttribute(colors, 3));
        const mat = new THREE.LineBasicMaterial({ vertexColors: true, transparent: true, opacity: 0.6 });
        return new THREE.LineSegments(geo, mat);
    }

    updateField(sim) {
        const n = 32;
        const sizeR = 2.5; // Larger view
        const sizeZ = 4.0;
        const pos = this.fieldPoints.geometry.getAttribute('position').array;
        const col = this.fieldPoints.geometry.getAttribute('color').array;
        let idx = 0;
        let maxB = 1e-4; // Adapt auto-scale

        // Pass 1: find maxB (magnitude) to normalize colors
        // But we want vectors. 
        // Let's define grid.

        // Z range: -1 to +4 (to cover nozzle)
        // R range: 0 to 2.5

        const zMin = -1.0;
        const zMax = 4.0;
        const rMax = 2.5;

        // Temp arrays to store B values to avoid re-sampling
        const samples = [];

        for (let j = 0; j < n; j++) {
            const z = zMin + (zMax - zMin) * (j / (n - 1));
            for (let i = 0; i < n; i++) {
                const r = rMax * (i / (n - 1));
                const B = sampleFieldAtPoint(Math.abs(r), z, sim);
                const m = Math.sqrt(B.Br * B.Br + B.Bz * B.Bz);
                if (m > maxB && isFinite(m)) maxB = m;
                samples.push({ r, z, Br: B.Br, Bz: B.Bz, m });
            }
        }

        // Pass 2: build lines
        // Vector scale length
        const segLen = 0.08;

        for (let k = 0; k < samples.length; k++) {
            const samp = samples[k];
            // Normalize B for direction
            let br = 0, bz = 1;
            if (samp.m > 1e-12) {
                br = samp.Br / samp.m;
                bz = samp.Bz / samp.m;
            }
            // Center of the vector
            const cz = samp.z; // X in view
            const cr = samp.r; // Y in view

            // Start and End
            // Visual mapping: z->x, r->y
            const x0 = cz - bz * segLen * 0.5;
            const y0 = cr - br * segLen * 0.5;
            const x1 = cz + bz * segLen * 0.5;
            const y1 = cr + br * segLen * 0.5;

            pos[idx * 3 + 0] = x0; pos[idx * 3 + 1] = y0; pos[idx * 3 + 2] = 0;
            pos[idx * 3 + 3] = x1; pos[idx * 3 + 4] = y1; pos[idx * 3 + 5] = 0;

            // Color based on magnitude (simple heatmap)
            // Blue (weak) -> Red (strong)
            const v = Math.min(1.0, samp.m / (maxB * 0.8 + 1e-9));
            const rcol = v;
            const gcol = Math.max(0, 0.5 - Math.abs(v - 0.5));
            const bcol = 1.0 - v;

            // Gradient along line? Or solid. Solid for now.
            col[idx * 3 + 0] = rcol; col[idx * 3 + 1] = gcol; col[idx * 3 + 2] = bcol;
            col[idx * 3 + 3] = rcol; col[idx * 3 + 4] = gcol; col[idx * 3 + 5] = bcol;

            idx += 2;
        }

        this.fieldPoints.geometry.attributes.position.needsUpdate = true;
        this.fieldPoints.geometry.attributes.color.needsUpdate = true;
    }

    makeLine(color) {
        const geo = new THREE.BufferGeometry();
        // Initial buffer size
        const maxVerts = 2000;
        geo.setAttribute('position', new THREE.BufferAttribute(new Float32Array(maxVerts * 3), 3));
        const mat = new THREE.LineBasicMaterial({ color });
        const mesh = new THREE.Line(geo, mat);
        mesh.frustumCulled = false;
        return mesh;
    }

    updateLines(sim) {
        const updateMesh = (mesh, coilList, isPlasma) => {
            // We want to visualize the "Shape".
            // 1. Rings (Draw3D.drawCircleAxis) -> LineSegments
            // 2. Profile (connecting centers) -> LineStrip (THREE.Line)
            // But 'mesh' is currently LineSegments.
            // Let's just draw the profile for the Cage and Plasma Surface to clearly show the Parabola/Sphere.
            // And maybe small rings or ticks?
            // User complained about "3D tubes" before. "Simple lines".

            // To fix "Rugby ball" and "Cone" perception, we MUST draw the profile curve.
            // But we only have a list of coils.
            // Let's populate the buffer with segments connecting coil centers.

            // If we use LineSegments, we need pairs.
            // (p0, p1), (p1, p2)...

            const vertices = [];

            // Draw Profile (Connectivity)
            if (coilList.length > 1) {
                for (let i = 0; i < coilList.length - 1; i++) {
                    const p0 = coilList[i];
                    const p1 = coilList[i + 1];
                    // Top profile (+r)
                    vertices.push(p0.z, p0.r, 0);
                    vertices.push(p1.z, p1.r, 0);

                    // Bottom profile (-r) - mirror for visualisation
                    vertices.push(p0.z, -p0.r, 0);
                    vertices.push(p1.z, -p1.r, 0);
                }
                // Close loop for Plasma?
                if (isPlasma) {
                    // Plasma surface is not simple loop in Z-ordering usually?
                    // makeSphericalPlasma orders them by angle.
                    // It is an arc from pole to pole.
                    // We don't close pLast to p0.
                }
            }

            // Draw 3D rings at each coil position
            // startV must be a UNIT vector - Draw3D multiplies it by the radius R
            const uaxis = { x: 1, y: 0, z: 0 };
            const nSegs = 32;
            coilList.forEach(p => {
                const center = { x: p.z, y: 0, z: 0 };
                const startV = { x: 0, y: 1, z: 0 };  // Unit vector in Y direction
                Draw3D.drawCircleAxis(nSegs, center, startV, uaxis, p.r, vertices);
            });

            const posAttr = mesh.geometry.attributes.position;
            if (vertices.length > posAttr.array.length) {
                mesh.geometry.setAttribute('position', new THREE.BufferAttribute(new Float32Array(vertices.length + 1000), 3));
            }
            mesh.geometry.attributes.position.array.set(vertices);
            mesh.geometry.attributes.position.needsUpdate = true;
            mesh.geometry.setDrawRange(0, vertices.length / 3);
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
