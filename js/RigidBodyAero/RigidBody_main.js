import { Renderer } from './OrientedParticles.js';
import { OrientedParticles_tex } from './OrientedParticles_tex.js';
import { UI } from './ui.js';
import { Physics } from './physics.js';
import { RigidBody } from './RigidBody.js';
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

        // Kick off async setup
        this.setup();
    }

    async setup() {
        // 3. Initialize Physics and Renderers
        const nParticles = this.ui.params.nParticles;
        const texSize = Math.ceil(Math.sqrt(nParticles));
        
        this.updateBodyDef();

        this.physics = new Physics(this.gl, nParticles, this.bodyDef);
        
        await this.renderer.initInstancedMesh(nParticles, this.ui.params.comOffset);
        
        this.rendererTex = new OrientedParticles_tex(this.renderer.scene, this.renderer.camera);
        await this.rendererTex.init(nParticles, texSize, this.ui.params.comOffset);
        
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

    updateBodyDef() {
        const offset = this.ui.params.comOffset; // 0 (back) to 1 (nose)
        // Nose is at z=1, Back is at z=0 (relative to local origin before shift)
        // We want to shift the COG such that it is at 'offset' between back and nose.
        // If offset=1, COG is at nose. If offset=0, COG is at back.
        // New local coordinates r' = r - r_com
        // r_com = [0, 0, offset]
        const z_shift = -offset;

        this.bodyDef = {
            points: new Float32Array([
                -0.5, 0.0, 0.0 + z_shift,  0.5, 0.0, 0.0 + z_shift,  0.0, 0.0, 1.0 + z_shift, // Wing
                 0.0, 0.0, 0.0 + z_shift,  0.0, 0.5, 0.0 + z_shift,  0.0, 0.0, 1.0 + z_shift  // Rudder
            ]),
            inertiaInv: [2.0, 2.0, 2.0]
        };
    }

    initCPUDebugGeometry() {
        this.cpuMeshGroup.clear();
        this.cpuLinesGroup.clear();

        // 1. Create a simple dart mesh for the CPU body
        const positions = [
            -0.5, 0, 0,   0, 0, 0,   0, 0, 1.0,  // Left Wing
             0.5, 0, 0,   0, 0, 0,   0, 0, 1.0,  // Right Wing
             0, 0.5, 0,   0, 0, 0,   0, 0, 1.0   // Rudder
        ];
        const colors = [
            1, 0, 0,  1, 0, 0,  1, 0, 0,
            0, 1, 0,  0, 1, 0,  0, 1, 0,
            0, 0, 1,  0, 0, 1,  0, 0, 1
        ];

        const geometry = new THREE.BufferGeometry();
        geometry.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
        geometry.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));
        const material = new THREE.MeshBasicMaterial({ vertexColors: true, side: THREE.DoubleSide, wireframe: false });
        this.cpuDartMesh = new THREE.Mesh(geometry, material);
        this.cpuMeshGroup.add(this.cpuDartMesh);

        // 2. Create lines from COG to force points
        const lineMaterial = new THREE.LineBasicMaterial({ color: 0x888888 });
        for (let i = 0; i < this.bodyDef.points.length; i += 3) {
            const lineGeom = new THREE.BufferGeometry().setFromPoints([
                new THREE.Vector3(0, 0, 0),
                new THREE.Vector3(this.bodyDef.points[i], this.bodyDef.points[i+1], this.bodyDef.points[i+2])
            ]);
            const line = new THREE.Line(lineGeom, lineMaterial);
            this.cpuLinesGroup.add(line);
        }
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
            this.renderer.updateGeometry(this.ui.params.comOffset);
        } else {
            await this.renderer.initInstancedMesh(nParticles, this.ui.params.comOffset);
        }

        if (this.rendererTex.mesh) {
            this.rendererTex.updateGeometry(this.ui.params.comOffset);
        } else {
            await this.rendererTex.init(nParticles, texSize, this.ui.params.comOffset);
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
        this.cpuDartMesh.position.copy(this.cpuBody.pos);
        this.cpuDartMesh.setRotationFromQuaternion(this.cpuBody.quat);

        // Update line group position/rotation
        this.cpuLinesGroup.position.copy(this.cpuBody.pos);
        this.cpuLinesGroup.setRotationFromQuaternion(this.cpuBody.quat);

        // Update force point arrows
        for (let i = 0; i < this.forceArrows.length; i++) {
            const arrow = this.forceArrows[i];
            const r_local = new THREE.Vector3(this.bodyDef.points[i*3], this.bodyDef.points[i*3+1], this.bodyDef.points[i*3+2]);
            const r_world = r_local.clone().applyMatrix3(R);
            const p_world = this.cpuBody.pos.clone().add(r_world);
            
            // Re-calculate local force for visualization (simplified drag only for now)
            const v_point = this.cpuBody.vel.clone().add(new THREE.Vector3().crossVectors(this.cpuBody.angVel, r_world));
            const vel_sq = v_point.lengthSq();
            const f_point = new THREE.Vector3(0,0,0);
            if (vel_sq > 0.0001) {
                f_point.copy(v_point).normalize().multiplyScalar(-0.01 * vel_sq);
            }
            
            arrow.position.copy(p_world);
            if (f_point.length() > 0.001) {
                arrow.setDirection(f_point.clone().normalize());
                arrow.setLength(f_point.length() * 10, this.headLen, this.headWid); // Scale for visibility
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

            // CPU Step
            this.cpuBody.step(dt, gravity);
            this.updateDebugVisuals();

            // GPU Step
            this.renderer.renderer.state.reset(); 
            this.physics.step(dt, gravity, this.bodyDef);
            
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
