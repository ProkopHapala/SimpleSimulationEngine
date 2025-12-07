// Simple Three.js-based renderer for the MHD demo (r-z view)
import * as THREE from 'https://cdn.jsdelivr.net/npm/three@0.160.0/build/three.module.js';
import { OrbitControls } from 'https://cdn.jsdelivr.net/npm/three@0.160.0/examples/jsm/controls/OrbitControls.js';
import { CSS2DRenderer, CSS2DObject } from 'https://cdn.jsdelivr.net/npm/three@0.160.0/examples/jsm/renderers/CSS2DRenderer.js';
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

        // WebGL Renderer
        this.renderer = new THREE.WebGLRenderer({ antialias: true });
        this.renderer.setSize(w, h);
        this.renderer.domElement.style.position = 'absolute';
        this.renderer.domElement.style.top = '0';
        this.renderer.domElement.style.left = '0';
        container.appendChild(this.renderer.domElement);

        // CSS2D Renderer for labels
        this.labelRenderer = new CSS2DRenderer();
        this.labelRenderer.setSize(w, h);
        this.labelRenderer.domElement.style.position = 'absolute';
        this.labelRenderer.domElement.style.top = '0';
        this.labelRenderer.domElement.style.left = '0';
        this.labelRenderer.domElement.style.pointerEvents = 'none'; // let clicks pass through to controls
        container.appendChild(this.labelRenderer.domElement);

        this.controls = new OrbitControls(this.camera, this.renderer.domElement);
        // OrbitControls attaches to WebGL canvas, but CSS canvas is on top overlaying it.
        // We need to attach controls to the container or CSS element?
        // Actually, if we set CSS pointerEvents none, we can attach to WebGL canvas behind it IF events propagate?
        // No, canvas catches them. But labelRenderer sits on top.
        // Better to attach controls to LABEL renderer element if it consumes events, OR set labelRenderer pointerEvents=none.
        // I set labelRenderer pointerEvents = none. So OrbitControls on renderer.domElement should work.
        this.controls.enableRotate = true;
        this.controls.enablePan = true;
        this.controls.enableZoom = true;
        this.controls.enableDamping = true;

        // Render options (controlled by UI)
        this.showField = true;
        this.showControlPoints = false;
        this.showRings = true;
        this.showProfile = true;
        this.fieldScale = 1.0;

        // Use LineSegments for coils to allow independent loops
        this.plasmaLine = this.makeSegments(0x00c2ff);
        this.cageLine = this.makeSegments(0xffa500);
        this.scLine = this.makeSegments(0x66ff66);
        this.scene.add(this.plasmaLine);
        this.scene.add(this.cageLine);
        this.scene.add(this.scLine);

        // Field: Grid Points + Vectors
        this.fieldPoints = this.makeFieldPoints();
        this.scene.add(this.fieldPoints);
        this.gridPoints = this.makeGridPoints();
        this.scene.add(this.gridPoints);

        // Control points visualization
        this.controlPointsMesh = this.makeControlPointsMesh();
        this.scene.add(this.controlPointsMesh);

        this.tempVec = new Vec3(); // reuse
        this.labels = []; // Store CSS2DObjects

        window.addEventListener('resize', () => this.onResize());
    }

    makeControlPointsMesh() {
        const geo = new THREE.BufferGeometry();
        const maxPoints = 200;
        geo.setAttribute('position', new THREE.BufferAttribute(new Float32Array(maxPoints * 3), 3));
        geo.setAttribute('color', new THREE.BufferAttribute(new Float32Array(maxPoints * 3), 3));
        const mat = new THREE.PointsMaterial({ size: 8, vertexColors: true, sizeAttenuation: false });
        const mesh = new THREE.Points(geo, mat);
        mesh.frustumCulled = false;
        return mesh;
    }

    makeGridPoints() {
        const n = 32;
        const geo = new THREE.BufferGeometry();
        // Just points for grid starts
        const positions = new Float32Array(n * n * 3);
        geo.setAttribute('position', new THREE.BufferAttribute(positions, 3));
        const mat = new THREE.PointsMaterial({ color: 0x444444, size: 2, sizeAttenuation: false });
        // Populate positions once (static grid)
        const zMin = -1.0, zMax = 4.0, rMax = 2.5;
        let idx = 0;
        for (let j = 0; j < n; j++) {
            const z = zMin + (zMax - zMin) * (j / (n - 1));
            for (let i = 0; i < n; i++) {
                const t = i / (n - 1);
                const r = -rMax + 2 * rMax * t;
                positions[idx * 3 + 0] = z;
                positions[idx * 3 + 1] = r;
                positions[idx * 3 + 2] = 0;
                idx++;
            }
        }
        return new THREE.Points(geo, mat);
    }

    updateControlPoints(sim) {
        if (!sim.controlPoints) return;

        const pos = this.controlPointsMesh.geometry.attributes.position.array;
        const col = this.controlPointsMesh.geometry.attributes.color.array;

        let idx = 0;
        for (const cp of sim.controlPoints) {
            // Position: z -> x, r -> y (matching field visualization)
            pos[idx * 3 + 0] = cp.z;
            pos[idx * 3 + 1] = cp.r;
            pos[idx * 3 + 2] = 0;

            // Color by type
            let r = 1, g = 0, b = 1; // Default: magenta
            if (cp.type === 'axis') { r = 1; g = 1; b = 0; } // Yellow
            else if (cp.type === 'vertex') { r = 0; g = 1; b = 0; } // Green
            else if (cp.type === 'focus') { r = 1; g = 0.5; b = 0; } // Orange
            else if (cp.type === 'cage') { r = 1; g = 0.6; b = 0.2; } // Orange-ish
            else if (cp.type === 'plasma-interior') { r = 0; g = 0.8; b = 1; } // Cyan

            col[idx * 3 + 0] = r;
            col[idx * 3 + 1] = g;
            col[idx * 3 + 2] = b;
            idx++;
        }

        this.controlPointsMesh.geometry.attributes.position.needsUpdate = true;
        this.controlPointsMesh.geometry.attributes.color.needsUpdate = true;
        this.controlPointsMesh.geometry.setDrawRange(0, sim.controlPoints.length);
    }

    makeSegments(color) {
        const geo = new THREE.BufferGeometry();
        // Initial buffer size
        const maxVerts = 4000; // Increased buffer size
        geo.setAttribute('position', new THREE.BufferAttribute(new Float32Array(maxVerts * 3), 3));
        const mat = new THREE.LineBasicMaterial({ color });
        const mesh = new THREE.LineSegments(geo, mat);
        mesh.frustumCulled = false;
        return mesh;
    }

    makeFieldPoints() {
        const n = 32;
        // Use LineSegments for vector field - 2 vertices per vector
        const positions = new Float32Array(n * n * 2 * 3);
        const colors = new Float32Array(n * n * 2 * 3);
        const geo = new THREE.BufferGeometry();
        geo.setAttribute('position', new THREE.BufferAttribute(positions, 3));
        geo.setAttribute('color', new THREE.BufferAttribute(colors, 3));
        const mat = new THREE.LineBasicMaterial({ vertexColors: true, transparent: true, opacity: 0.8 });
        return new THREE.LineSegments(geo, mat);
    }

    updateLabels(sim) {
        // Clear existing
        this.labels.forEach(l => {
            this.scene.remove(l);
            if (l.element && l.element.parentNode) l.element.parentNode.removeChild(l.element);
        });
        this.labels = [];

        const addLabel = (r, z, I, name) => {
            if (Math.abs(I) < 1.0 && name !== 'SC') return; // Skip near-zero currents except SC

            const div = document.createElement('div');
            div.className = 'label';
            // Show kA
            div.textContent = `I=${(I / 1000).toFixed(1)}kA`;

            const label = new CSS2DObject(div);
            // Position: z -> x, r -> y
            label.position.set(z, r, 0);
            this.scene.add(label);
            this.labels.push(label);
        };

        // SC Coils
        sim.scCoils.forEach(c => addLabel(c.r, c.z, c.I, 'SC'));

        // Cage Coils
        sim.cageCoils.forEach((c, i) => {
            // Only label a few to avoid clutter? Or all? User wants "to each coil"
            // Let's do all for now, maybe decimate if too many
            const I = sim.cageCurrents[i] || 0;
            addLabel(c.r, c.z, I, 'Cage');
        });

        // Plasma Coils? Usually too many overlapping nodes.
        // Maybe just center?
        // Let's sum total plasma current and put label at center?
        // Or user insists "each coil", we might need to be careful.
        // Let's do middle node of plasma.
        if (sim.plasmaNodes.length > 0) {
            const mid = Math.floor(sim.plasmaNodes.length / 2);
            const p = sim.plasmaNodes[mid];
            const I = sim.plasmaCurrents[mid] || 0;
            // Total plasma current?
            const totalI = sim.plasmaCurrents.reduce((a, b) => a + b, 0);
            addLabel(p.r, p.z, totalI, 'PlasmaTotal');
        }
    }

    updateField(sim) {
        const n = 32;
        const zMin = -1.0;
        const zMax = 4.0;
        const rMax = 2.5;

        const pos = this.fieldPoints.geometry.getAttribute('position').array;
        const col = this.fieldPoints.geometry.getAttribute('color').array;

        let idx = 0;
        const scale = this.fieldScale * 0.05; // Adjustable scale factor

        // Symmetric Grid Visualization
        // We draw vectors starting at grid points
        for (let j = 0; j < n; j++) {
            const z = zMin + (zMax - zMin) * (j / (n - 1));
            for (let i = 0; i < n; i++) {
                const t = i / (n - 1);
                const r = -rMax + 2 * rMax * t;

                const B = sampleFieldAtPoint(Math.abs(r), z, sim);

                // Flip radial component for visualization in lower half
                // Br from physics is always "outward" (positive r direction)
                // If r < 0 (lower half), "outward" is -y direction.
                // If Br is positive (outward), then dy should be negative.
                // Wait.
                // Physics: Br is component in \hat{r}.
                // Position vector R = r * \hat{r} + z * \hat{z}.
                // B = Br * \hat{r} + Bz * \hat{z}.
                // In Cartesian (y, x):
                // If point is at y_coord (which maps to r), x_coord (maps to z).
                // If y_coord > 0: \hat{r} is +y. B_y = Br.
                // If y_coord < 0: \hat{r} is -y. B_y = -Br?
                // Actually, rotational symmetry.
                // If Br > 0 (outward), and we are at y=-5, force is outward (down), so B_y should be -5.
                // Yes. B vector should point away from axis.

                const sign = r >= 0 ? 1 : -1;
                const B_y = B.Br * sign;
                const B_x = B.Bz; // Bz is along Z axis (+x in view)

                const mag = Math.sqrt(B_y * B_y + B_x * B_x);

                // Start Point (at grid)
                const sx = z;
                const sy = r;

                // End Point = Start + Vector * Scale
                // Asymmetric: Vector starts at grid point
                const ex = sx + B_x * scale;
                const ey = sy + B_y * scale;

                pos[idx * 3 + 0] = sx; pos[idx * 3 + 1] = sy; pos[idx * 3 + 2] = 0;
                pos[idx * 3 + 3] = ex; pos[idx * 3 + 4] = ey; pos[idx * 3 + 5] = 0;

                // Color based on magnitude
                // Log scale for better dynamic range?
                // Use a clamped linear scale for now
                const maxVal = 0.5;
                const v = Math.min(1.0, mag / maxVal);

                // Heatmap: Blue -> Cyan -> Green -> Yellow -> Red
                // Simple: Blue -> Red
                const rcol = v;
                const gcol = v < 0.5 ? v * 2 : (1 - v) * 2;
                const bcol = 1.0 - v;

                col[idx * 3 + 0] = rcol; col[idx * 3 + 1] = gcol; col[idx * 3 + 2] = bcol;
                col[idx * 3 + 3] = rcol; col[idx * 3 + 4] = gcol; col[idx * 3 + 5] = bcol;

                idx += 2;
            }
        }

        this.fieldPoints.geometry.attributes.position.needsUpdate = true;
        this.fieldPoints.geometry.attributes.color.needsUpdate = true;
    }

    makeLine(color) {
        const geo = new THREE.BufferGeometry();
        // Initial buffer size
        const maxVerts = 4000;
        geo.setAttribute('position', new THREE.BufferAttribute(new Float32Array(maxVerts * 3), 3));
        const mat = new THREE.LineBasicMaterial({ color });
        const mesh = new THREE.Line(geo, mat);
        mesh.frustumCulled = false;
        return mesh;
    }

    updateLines(sim) {
        const updateMesh = (mesh, coilList, isPlasma) => {
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
            // Check resize if needed
            if (vertices.length > posAttr.array.length) {
                // Warning: buffer resize is expensive, better start large
                // Re-create geometry would be better but let's assume maxVerts=4000 is enough
                // If not, we should reallocate.
                const newArr = new Float32Array(vertices.length + 1000);
                mesh.geometry.setAttribute('position', new THREE.BufferAttribute(newArr, 3));
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
        if (this.labelRenderer) this.labelRenderer.setSize(w, h);
    }

    render(sim) {
        // Update visibility based on flags
        this.plasmaLine.visible = this.showRings || this.showProfile;
        this.cageLine.visible = this.showRings || this.showProfile;
        this.scLine.visible = this.showRings || this.showProfile;
        this.fieldPoints.visible = this.showField;
        this.gridPoints.visible = this.showField;
        this.controlPointsMesh.visible = this.showControlPoints;

        this.updateLines(sim);
        if (this.showField) {
            this.updateField(sim);
        }
        if (this.showControlPoints) {
            this.updateControlPoints(sim);
        }

        // CSS Labels update
        // We call updateLabels explicitly from UI on change, but maybe also here if currents change continuously?
        // Sim currents change during solve step only.
        // But let's leave valid logic in updateLabels and call it if needed.

        this.renderer.render(this.scene, this.camera);
        if (this.labelRenderer) this.labelRenderer.render(this.scene, this.camera);
    }
}
