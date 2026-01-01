import { UI } from './ui.js';
import { Physics } from './physics.js';
import { Renderer } from './render.js';

class App {
    constructor() {
        this.container = document.getElementById('container');
        this.ui = new UI({
            onReset: () => this.reset()
        });

        const bodyDef = {
            points: new Float32Array([
                -0.5, 0, 0.5,  0.5, 0, 0.5,  0, 0, -0.5, // Wing
                0, 0, 0.5,  0, 0.5, 0.5,  0, 0, -0.5  // Rudder
            ]),
            inertiaInv: [2.0, 2.0, 2.0] // Simple uniform inertia for now
        };

        this.renderer = new Renderer(this.container, this.ui.params.nParticles, Math.ceil(Math.sqrt(this.ui.params.nParticles)));
        const gl = this.renderer.renderer.getContext();
        if (!gl.getExtension('EXT_color_buffer_float')) {
            alert('EXT_color_buffer_float not supported');
        }

        this.physics = new Physics(gl, this.ui.params.nParticles, bodyDef);
        this.lastTime = performance.now();
        
        this.renderer.material.uniforms.u_static_test.value = 1.0; // Start with static test
        this.frameCounter = 0;
        this.physics.initShaders().then(() => {
            this.animate();
        });
    }

    reset() {
        this.physics.nParticles = this.ui.params.nParticles;
        this.physics.size = Math.ceil(Math.sqrt(this.physics.nParticles));
        // Re-init textures if size changed
        this.physics.initTextures();
        this.physics.initFramebuffers();
        this.physics.reset();
        
        // Update renderer
        this.renderer.scene.remove(this.renderer.mesh);
        this.renderer.nParticles = this.physics.nParticles;
        this.renderer.size = this.physics.size;
        this.renderer.initInstancedMesh();
        this.renderer.material.uniforms.u_static_test.value = 1.0;
    }

    animate() {
        requestAnimationFrame(() => this.animate());

        const now = performance.now();
        const dt_real = (now - this.lastTime) / 1000;
        this.lastTime = now;
        this.ui.updateFPS(1 / dt_real);

        if (!this.ui.params.paused) {
            // Take multiple steps if needed or just one with fixed dt
            this.physics.step(this.ui.params.dt, this.ui.params.gravity);
            
            this.frameCounter++;
            if(this.frameCounter % 60 === 0){
                console.log(`--- Frame ${this.frameCounter} ---`);
                this.physics.debugReadback();
            }
        }

        const out = this.physics.getOutput();
        this.renderer.render(out.pos, out.quat);
    }
}

new App();
