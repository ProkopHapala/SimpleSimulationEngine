import { Renderer } from './OrientedParticles.js';
import { OrientedParticles_tex } from './OrientedParticles_tex.js';
import { UI } from './ui.js';
import { Physics } from './physics.js';
import { RigidBody } from './RigidBody.js';
import { geomFromString, paramsFromString, buildPoints, computePointForce } from './AeroSurface.js';
import * as THREE from 'three';
import { OrbitControls } from 'three/examples/jsm/controls/OrbitControls.js';

class App {
    constructor() {
        this.container = document.getElementById('container');
        
        // 1. Initialize Renderer (Scene/Camera/WebGL)
        this.renderer = new Renderer(this.container);
        this.gl = this.renderer.renderer.getContext();
        
        // Add OrbitControls
        this.controls = new OrbitControls(this.renderer.camera, this.renderer.renderer.domElement);
        this.controls.enableDamping = true;

        // 2. Initialize UI
        this.ui = new UI({
            onReset: () => this.reset(),
            onRenderModeChange: (mode) => this.setRenderMode(mode),
            onRunToggle: (running) => { this.running = running; },
            onToggleShowGPU: (show) => this.setShowGPU(show),
            onToggleShowCPU: (show) => this.setShowCPU(show),
            onToggleCPUDebug: () => this.updateDebugVisuals(),
            onHeadSizeChange: (size) => this.updateArrowHeadSize(size)
        });

        // Debug Visualizers
        this.debugHelpers = new THREE.Group();
        this.renderer.scene.add(this.debugHelpers);
        
        // CPU Debug Mesh
        this.cpuMeshGroup = new THREE.Group();
        this.debugHelpers.add(this.cpuMeshGroup);

        // CPU Debug Lines (connecting points to CoM)
        this.cpuLinesGroup = new THREE.Group();
        this.debugHelpers.add(this.cpuLinesGroup);

        this.forceArrows = [];
        this.headLen = 0.02; this.headWid = 0.015;
        this.torqueArrow = new THREE.ArrowHelper(new THREE.Vector3(0,1,0), new THREE.Vector3(0,0,0), 0, 0xff00ff, this.headLen, this.headWid);
        this.totalForceArrow = new THREE.ArrowHelper(new THREE.Vector3(0,1,0), new THREE.Vector3(0,0,0), 0, 0xffff00, this.headLen, this.headWid);
        this.debugHelpers.add(this.torqueArrow, this.totalForceArrow);

        // Resizer for UI panel
        this.panel = document.getElementById('ui-panel');
        this.resizer = document.getElementById('ui-resizer');
        this.initResizer();

        // Kick off async setup
        this.setup();
    }

    async setup() {
        // 3. Initialize Physics and Renderers
        const nParticles = this.ui.params.nParticles;
        const texSize = Math.ceil(Math.sqrt(nParticles));
        
        this.updateBodyDef();

        this.physics = new Physics(this.gl, nParticles, this.bodyDef);
        
        await this.renderer.initInstancedMesh(nParticles, this.bodyDef);
        
        this.rendererTex = new OrientedParticles_tex(this.renderer.scene, this.renderer.camera);
        await this.rendererTex.init(nParticles, texSize, this.bodyDef);
        
        this.setRenderMode(this.ui.params.renderMode || 'texture');
        this.running = true;
        this.lastTime = performance.now();

        // 4. Initial Data
        await this.physics.initShaders();
        this.reset();
        this.animate();
    }

    updateArrowHeadSize(size) {
        this.headLen = size;
        this.headWid = size * 0.75;
        this.torqueArrow.setLength(this.torqueArrow.line.scale.y, this.headLen, this.headWid);
        this.totalForceArrow.setLength(this.totalForceArrow.line.scale.y, this.headLen, this.headWid);
        this.forceArrows.forEach(a => {
            a.setLength(a.line.scale.y, this.headLen, this.headWid);
        });
    }

    initResizer() {
        if (!this.resizer || !this.panel) return;
        let isResizing = false;

        const onMouseMove = (e) => {
            if (!isResizing) return;
            const newWidth = Math.min(Math.max(window.innerWidth - e.clientX, 260), window.innerWidth * 0.9);
            this.panel.style.width = `${newWidth}px`;
            this.resizer.style.right = `${newWidth}px`;
        };

        this.resizer.addEventListener('mousedown', (e) => {
            isResizing = true;
            document.body.style.cursor = 'col-resize';
            e.preventDefault();
        });

        window.addEventListener('mousemove', onMouseMove);
        window.addEventListener('mouseup', () => {
            if (!isResizing) return;
            isResizing = false;
            document.body.style.cursor = 'default';
        });
    }

    updateBodyDef() {
        const offset = this.ui.params.comOffset; // 0 (back) to 1 (nose)
        const z_shift = -offset;

        // Parse Wing Geometry and Params via helper
        const geomLines = this.ui.elements.wingGeometry.value;
        const aeroLines = this.ui.elements.aeroParams.value;

        const geom = geomFromString(geomLines, z_shift);
        const aero = paramsFromString(aeroLines);
        const points = buildPoints(geom, aero);

        // If parsing fails, fall back to a default wing
        if (points.length === 0) {
            points.push({
                pos: new THREE.Vector3(0, 0, z_shift),
                t: new THREE.Vector3(0, 0, 1),
                n: new THREE.Vector3(0, 1, 0),
                params: { area: 1.0, CD0: 0.02, dCD: 0.9, dCDS: 0.9, dCL: 6.28, dCLS: 2.82, sStall: 0.16, wStall: 0.08 }
            });
        }

        this.bodyDef = {
            points: points,
            inertiaInv: [2.0, 2.0, 2.0],
            z_shift: z_shift
        };

        if (this.ui && typeof this.ui.setWingOptions === 'function') {
            this.ui.setWingOptions(points.length);
        }
    }

    initCPUDebugGeometry() {
        this.cpuMeshGroup.clear();
        this.cpuLinesGroup.clear();

        const z_shift = this.bodyDef.z_shift;

        // Placeholder for COG marker (to avoid undefined refs)
        this.cpuDartMesh = new THREE.Group();
        this.cpuMeshGroup.add(this.cpuDartMesh);

        // Create a mesh for each wing segment
        const lineMaterial = new THREE.LineBasicMaterial({ color: 0x888888 });
        
        this.bodyDef.points.forEach(p => {
            // Visualize the wing chord as a small triangle/line
            const geom = new THREE.BufferGeometry();
            const wingWidth = 0.5;
            const wingLength = 0.5;
            
            // local coords of the segment
            const side = new THREE.Vector3().crossVectors(p.t, p.n).normalize();
            
            const v0 = p.pos.clone().addScaledVector(side, -wingWidth);
            const v1 = p.pos.clone().addScaledVector(side,  wingWidth);
            const v2 = p.pos.clone().addScaledVector(p.t,  wingLength); // tip forward along t
            
            const positions = new Float32Array([
                v0.x, v0.y, v0.z,
                v1.x, v1.y, v1.z,
                v2.x, v2.y, v2.z
            ]);
            
            geom.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
            const mesh = new THREE.Mesh(geom, new THREE.MeshBasicMaterial({ color: 0x4444ff, wireframe: true, side: THREE.DoubleSide }));
            this.cpuMeshGroup.add(mesh);

            // CoM to wing connection
            const lineGeom = new THREE.BufferGeometry().setFromPoints([
                new THREE.Vector3(0, 0, 0),
                p.pos
            ]);
            const line = new THREE.Line(lineGeom, lineMaterial);
            this.cpuLinesGroup.add(line);
        });
    }

    async reset() {
        const nParticles = this.ui.params.nParticles;
        const texSize = Math.ceil(Math.sqrt(nParticles));
        const scatterRad = (this.ui.params.rotScatterDeg || 0) * Math.PI / 180;

        const { speedBase, speedSpread } = this.ui.params;
        this.updateBodyDef();
        this.physics.nParticles = nParticles;
        this.physics.size = texSize;
        this.physics.initTextures();
        this.physics.initFramebuffers();
        this.physics.reset(scatterRad, speedBase, speedSpread);

        const posArray = new Float32Array(nParticles * 3);
        const quatArray = new Float32Array(nParticles * 4);
        
        // Reuse renderers; just update geometry for new CoM offset (avoid duplicate meshes)
        if (this.renderer.mesh) {
            this.renderer.updateGeometry(this.bodyDef);
        } else {
            await this.renderer.initInstancedMesh(nParticles, this.bodyDef);
        }

        if (this.rendererTex.mesh) {
            this.rendererTex.updateGeometry(this.bodyDef);
        } else {
            await this.rendererTex.init(nParticles, texSize, this.bodyDef);
        }

        // SYNC RANDOM SEED/VALUES FOR FIRST PARTICLE
        const seed_pos = [ (Math.random() - 0.5) * 5.0, 5.0 + (Math.random() - 0.5) * 5.0, (Math.random() - 0.5) * 5.0 ];
        const axis = [Math.random()-0.5, Math.random()-0.5, Math.random()-0.5];
        let len = Math.sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);
        if(len > 0) { axis[0]/=len; axis[1]/=len; axis[2]/=len; } else { axis[0]=1; axis[1]=0; axis[2]=0; }
        const angle = (Math.random() - 0.5) * scatterRad;
        const c = Math.cos(angle/2);
        const s = Math.sin(angle/2);
        const seed_quat = [axis[0]*s, axis[1]*s, axis[2]*s, c];
        const seed_speed = (Math.random() - 0.5) * speedSpread + speedBase;

        for (let i = 0; i < nParticles; i++) {
            if (i === 0) {
                posArray[0] = seed_pos[0]; posArray[1] = seed_pos[1]; posArray[2] = seed_pos[2];
                quatArray[0] = seed_quat[0]; quatArray[1] = seed_quat[1]; quatArray[2] = seed_quat[2]; quatArray[3] = seed_quat[3];
            } else {
                posArray[i*3+0] = (Math.random() - 0.5) * 5.0;
                posArray[i*3+1] = 5.0 + (Math.random() - 0.5) * 5.0;
                posArray[i*3+2] = (Math.random() - 0.5) * 5.0;

                const axis_i = [Math.random()-0.5, Math.random()-0.5, Math.random()-0.5];
                const len_i = Math.sqrt(axis_i[0]*axis_i[0] + axis_i[1]*axis_i[1] + axis_i[2]*axis_i[2]);
                if(len_i > 0) { axis_i[0]/=len_i; axis_i[1]/=len_i; axis_i[2]/=len_i; } else { axis_i[0]=1; axis_i[1]=0; axis_i[2]=0; }
                const angle_i = (Math.random() - 0.5) * scatterRad;
                const c_i = Math.cos(angle_i/2);
                const s_i = Math.sin(angle_i/2);
                quatArray[i*4+0] = axis_i[0]*s_i; quatArray[i*4+1] = axis_i[1]*s_i; quatArray[i*4+2] = axis_i[2]*s_i; quatArray[i*4+3] = c_i;
            }
        }
        this.renderer.updateInstances(posArray, quatArray);
        this.updateRendererTextures();

        // Initialize CPU RigidBody for particle 0 (EXACT MATCH with physics.js:reset)
        const initialPos = [0, 5, 0];
        const initialQuat = [0, 0, 0, 1];
        const initialVel = [0, 0, speedBase];
        const initialAngVel = [0, 0, 0];
        
        this.cpuBody = new RigidBody(
            initialPos, initialVel, initialQuat, initialAngVel, 
            1.0, this.bodyDef.inertiaInv, this.bodyDef.points
        );

        this.initCPUDebugGeometry();

        // Setup debug force arrows
        this.debugHelpers.children.filter(c => c.name === 'forcePoint').forEach(c => this.debugHelpers.remove(c));
        this.forceArrows = [];
        for (let i = 0; i < this.bodyDef.points.length / 3; i++) {
            const arrow = new THREE.ArrowHelper(new THREE.Vector3(0,1,0), new THREE.Vector3(0,0,0), 0, 0x00ff00, this.headLen, this.headWid);
            arrow.name = 'forcePoint';
            this.forceArrows.push(arrow);
            this.debugHelpers.add(arrow);
        }

        this.setShowGPU(this.ui.params.showGPU);
        this.setShowCPU(this.ui.params.showCPU);
    }

    updateRendererTextures() {
        if (!this.ui.params.showGPU) return;
        const out = this.physics.getOutput();
        const texSize = this.physics.size;
        
        const wrapGLTexture = (glTex) => {
            const tex = new THREE.Texture();
            const properties = this.renderer.renderer.properties.get(tex);
            properties.__webglTexture = glTex;
            properties.__webglInit = true;
            tex.image = { width: texSize, height: texSize };
            tex.format = THREE.RGBAFormat;
            tex.type = THREE.FloatType;
            tex.minFilter = THREE.NearestFilter;
            tex.magFilter = THREE.NearestFilter;
            tex.needsUpdate = false;
            return tex;
        };

        const posTex = wrapGLTexture(out.pos);
        const quatTex = wrapGLTexture(out.quat);
        this.rendererTex.updateTextures(posTex, quatTex);
    }

    setRenderMode(mode) {
        this.renderMode = mode;
        this.setShowGPU(this.ui.params.showGPU);
        this.setShowCPU(this.ui.params.showCPU);
    }

    setShowGPU(show) {
        console.log("setShowGPU CALLED - value:", show, "mode:", this.renderMode);
        this.ui.params.showGPU = show;
        
        // No longer removing from scene to avoid state mess. 
        // We will use if-guards in the render loop or .visible flags.
        const showBuf = show && (this.renderMode === 'buffer');
        const showTex = show && (this.renderMode === 'texture');
        
        if (this.renderer && this.renderer.mesh) {
            this.renderer.mesh.visible = showBuf;
        }
        
        if (this.rendererTex && this.rendererTex.mesh) {
            this.rendererTex.mesh.visible = showTex;
        }
    }

    setShowCPU(show) {
        console.log("setShowCPU called with:", show);
        this.ui.params.showCPU = show;
        this.debugHelpers.visible = show;
    }

    updateDebugVisuals() {
        if (!this.cpuBody) return;
        
        // Update visibility based on UI toggles
        this.cpuMeshGroup.visible = this.ui.params.showCPUMesh;
        this.cpuLinesGroup.visible = this.ui.params.showCPULines;

        const Rm4 = new THREE.Matrix4().makeRotationFromQuaternion(this.cpuBody.quat);
        const R = new THREE.Matrix3().setFromMatrix4(Rm4);
        
        // Update mesh position/rotation
        if (this.cpuDartMesh) {
            this.cpuDartMesh.position.copy(this.cpuBody.pos);
            this.cpuDartMesh.setRotationFromQuaternion(this.cpuBody.quat);
        }

        // Update line group position/rotation
        this.cpuLinesGroup.position.copy(this.cpuBody.pos);
        this.cpuLinesGroup.setRotationFromQuaternion(this.cpuBody.quat);

        const v_wind = this.ui.params.windTunnelMode ? new THREE.Vector3(0, 0, -5.0) : new THREE.Vector3(0, 0, 0);
        for (let i = 0; i < this.forceArrows.length; i++) {
            const arrow = this.forceArrows[i];
            const p = this.bodyDef.points[i];
            if (!p) continue;

            const { force, position } = computePointForce(p, this.cpuBody, v_wind, R);
            arrow.position.copy(position);
            if (force.length() > 0.001) {
                arrow.setDirection(force.clone().normalize());
                arrow.setLength(force.length() * 0.1, this.headLen, this.headWid); // Scaled for visibility
                arrow.visible = true;
            } else {
                arrow.visible = false;
            }
        }

        // Total Force
        this.totalForceArrow.position.copy(this.cpuBody.pos);
        if (this.cpuBody.force.length() > 0.001) {
            this.totalForceArrow.setDirection(this.cpuBody.force.clone().normalize());
            this.totalForceArrow.setLength(this.cpuBody.force.length() * 5);
            this.totalForceArrow.visible = true;
        } else {
            this.totalForceArrow.visible = false;
        }

        // Total Torque
        this.torqueArrow.position.copy(this.cpuBody.pos);
        if (this.cpuBody.torque.length() > 0.001) {
            this.torqueArrow.setDirection(this.cpuBody.torque.clone().normalize());
            this.torqueArrow.setLength(this.cpuBody.torque.length() * 5);
            this.torqueArrow.visible = true;
        } else {
            this.torqueArrow.visible = false;
        }
    }

    animate() {
        requestAnimationFrame(() => this.animate());
        
        const now = performance.now();
        const dt_real = (now - this.lastTime) / 1000;
        this.lastTime = now;
        this.ui.updateFPS(1 / dt_real);

        if (this.running && !this.ui.params.paused) {
            const dt = this.ui.params.dt;
            const gravity = this.ui.params.gravity;
            const substeps = Math.max(1, this.ui.params.substeps || 1);
            const subDt = dt / substeps;

            for (let s = 0; s < substeps; s++) {
                // CPU Step
                this.cpuBody.step(subDt, gravity, this.ui.params.windTunnelMode);
                // GPU Step
                this.renderer.renderer.state.reset(); 
                this.physics.step(subDt, gravity, this.bodyDef);
            }

            this.updateDebugVisuals();
            this.ui.drawPlots(this.bodyDef);
            
            this.gl.bindFramebuffer(this.gl.FRAMEBUFFER, null);
            this.gl.bindVertexArray(null);
            if (this.gl.drawBuffers) this.gl.drawBuffers([this.gl.BACK]);
            this.gl.useProgram(null);
            this.gl.bindTexture(this.gl.TEXTURE_2D, null);
            this.gl.activeTexture(this.gl.TEXTURE0);
            
            if (this.renderer.renderer.resetState) {
                this.renderer.renderer.resetState();
            } else {
                this.renderer.renderer.state.reset();
            }

            if (this.ui.params.verbosity > 0) {
                this.physics.readback(this.ui.params.verbosity);
            }
        }

        // --- VISIBILITY GUARDS FOR RENDERING ---
        const showGPU = this.ui.params.showGPU;
        const showCPU = this.ui.params.showCPU;

        // Ensure visibility is forced every frame to the correct state
        const shouldShowBuf = showGPU && (this.renderMode === 'buffer');
        const shouldShowTex = showGPU && (this.renderMode === 'texture');

        if (this.renderer && this.renderer.mesh) {
            this.renderer.mesh.visible = shouldShowBuf;
            if (this.renderer.mesh.material) this.renderer.mesh.material.visible = shouldShowBuf;
        }

        if (this.rendererTex && this.rendererTex.mesh) {
            this.rendererTex.mesh.visible = shouldShowTex;
            if (this.rendererTex.mesh.material) this.rendererTex.mesh.material.visible = shouldShowTex;
        }

        if (this.debugHelpers) {
            this.debugHelpers.visible = showCPU;
        }

        if (showGPU && this.renderMode === 'texture') {
            this.updateRendererTextures();
        }

        this.controls.update();

        // 5. FINAL RENDER CALL
        this.gl.bindFramebuffer(this.gl.FRAMEBUFFER, null);
        if (this.gl.drawBuffers) this.gl.drawBuffers([this.gl.BACK]);
        
        // Debug every ~100 frames
        if (Math.random() < 0.01) {
            console.log("LOOP - showGPU:", showGPU, "mode:", this.renderMode, "visibleBuf:", this.renderer?.mesh?.visible, "visibleTex:", this.rendererTex?.mesh?.visible);
        }

        this.renderer.renderer.clear(); // Manual clear to ensure no ghosting
        this.renderer.renderer.render(this.renderer.scene, this.renderer.camera);
    }
}

new App();
