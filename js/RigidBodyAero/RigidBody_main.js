import { Renderer } from './OrientedParticles.js';
import { OrientedParticles_tex } from './OrientedParticles_tex.js';
import { UI } from './ui.js';
import { Physics } from './physics.js';
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
            onRunToggle: (running) => { this.running = running; }
        });

        // 3. Initialize Physics and Renderers
        const nParticles = 16;
        const texSize = 4; // 4x4 for 16 particles
        
        const bodyDef = {
            points: new Float32Array([
                -0.5, 0, 0.5,  0.5, 0, 0.5,  0, 0, -0.5, // Wing
                0, 0, 0.5,  0, 0.5, 0.5,  0, 0, -0.5  // Rudder
            ]),
            inertiaInv: [2.0, 2.0, 2.0]
        };

        this.physics = new Physics(this.gl, nParticles, bodyDef);
        
        this.renderer.initInstancedMesh(nParticles);
        
        this.rendererTex = new OrientedParticles_tex(this.renderer.scene, this.renderer.camera);
        this.rendererTex.init(nParticles, texSize);
        this.rendererTex.setVisibility(false);

        this.renderMode = 'buffer';
        this.running = false;
        this.lastTime = performance.now();

        // 4. Initial Data
        this.physics.initShaders().then(() => {
            this.reset();
            this.animate();
        });
    }

    reset() {
        const nParticles = 16;
        const texSize = 4;
        const scatterRad = (this.ui.params.rotScatterDeg || 0) * Math.PI / 180;

        this.physics.reset(scatterRad);

        const posArray = new Float32Array(nParticles * 3);
        const quatArray = new Float32Array(nParticles * 4);
        
        // We can use random data for initial buffer-based static test
        for (let i = 0; i < nParticles; i++) {
            posArray[i*3+0] = (Math.random() - 0.5) * 5.0;
            posArray[i*3+1] = (Math.random() - 0.5) * 5.0;
            posArray[i*3+2] = (Math.random() - 0.5) * 5.0;

            // random axis
            const axis = [Math.random()-0.5, Math.random()-0.5, Math.random()-0.5];
            const len = Math.sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);
            if(len > 0) { axis[0]/=len; axis[1]/=len; axis[2]/=len; } else { axis[0]=1; axis[1]=0; axis[2]=0; }
            const angle = (Math.random() - 0.5) * scatterRad;
            const c = Math.cos(angle/2);
            const s = Math.sin(angle/2);
            const q = [axis[0]*s, axis[1]*s, axis[2]*s, c];
            quatArray[i*4+0] = q[0];
            quatArray[i*4+1] = q[1];
            quatArray[i*4+2] = q[2];
            quatArray[i*4+3] = q[3];
        }
        this.renderer.updateInstances(posArray, quatArray);

        this.updateRendererTextures();
    }

    updateRendererTextures() {
        const out = this.physics.getOutput();
        // Here we convert WebGLTexture to THREE.Texture
        // Three.js doesn't easily wrap external WebGLTextures for general use in ShaderMaterial
        // but we can create a Raw Shader if needed. 
        // For now, we will try to pass the raw textures by wrapping them.
        
        const texSize = this.physics.size;
        
        if (!this.threeTexPos) {
            this.threeTexPos = new THREE.FramebufferTexture(texSize, texSize, THREE.RGBAFormat);
            this.threeTexQuat = new THREE.FramebufferTexture(texSize, texSize, THREE.RGBAFormat);
        }

        // NOTE: FramebufferTexture is usually for copying from screen.
        // To use raw GL textures from Physics class in THREE.ShaderMaterial, 
        // we might need to use a custom property or a DataTexture with a fake __webglTexture.
        
        const wrapGLTexture = (glTex) => {
            const tex = new THREE.Texture();
            const properties = this.renderer.renderer.properties.get(tex);
            properties.__webglTexture = glTex;
            properties.__webglInit = true; // Mark as initialized
            tex.format = THREE.RGBAFormat;
            tex.type = THREE.FloatType;
            tex.minFilter = THREE.NearestFilter;
            tex.magFilter = THREE.NearestFilter;
            tex.needsUpdate = true; // Trigger internal update logic
            return tex;
        };

        const posTex = wrapGLTexture(out.pos);
        const quatTex = wrapGLTexture(out.quat);
        
        this.rendererTex.updateTextures(posTex, quatTex);
    }

    setRenderMode(mode) {
        this.renderMode = mode;
        if (mode === 'buffer') {
            this.renderer.mesh.visible = true;
            this.rendererTex.setVisibility(false);
        } else {
            this.renderer.mesh.visible = false;
            this.rendererTex.setVisibility(true);
        }
    }

    animate() {
        requestAnimationFrame(() => this.animate());
        
        const now = performance.now();
        const dt_real = (now - this.lastTime) / 1000;
        this.lastTime = now;
        this.ui.updateFPS(1 / dt_real);

        if (this.running && !this.ui.params.paused) {
            // Rebind THREE's state to prevent conflicts when we use raw GL
            this.renderer.renderer.state.reset(); 
            
            this.physics.step(this.ui.params.dt, this.ui.params.gravity);
            
            // After raw GL calls, tell Three.js to reset its state tracking
            this.gl.bindFramebuffer(this.gl.FRAMEBUFFER, null);
            this.gl.bindVertexArray(null);
            // Reset Three.js state tracking and GL draw buffers
            if (this.gl.drawBuffers) this.gl.drawBuffers([this.gl.BACK]);
            if (this.renderer.renderer.resetState) {
                this.renderer.renderer.resetState();
            } else {
                this.renderer.renderer.state.reset();
            }
            
            if (this.renderMode === 'texture') {
                this.updateRendererTextures();
            }
        }

        this.controls.update();
        // Ensure default framebuffer/draw buffers before Three renders
        this.gl.bindFramebuffer(this.gl.FRAMEBUFFER, null);
        if (this.gl.drawBuffers) this.gl.drawBuffers([this.gl.BACK]);
        this.renderer.render();
    }
}

new App();
