import { Renderer } from './render.js';

class App {
    constructor() {
        this.container = document.getElementById('container');
        
        // 1. Initialize Renderer
        this.renderer = new Renderer(this.container);
        
        // 2. Initialize Instances (Step 1.8)
        const nParticles = 16;
        this.renderer.initInstancedMesh(nParticles);
        
        // 3. Generate Random Data
        const posArray = new Float32Array(nParticles * 3);
        const quatArray = new Float32Array(nParticles * 4);
        
        console.log(`--- Initializing ${nParticles} Particles (Step 1.8) ---`);
        
        for (let i = 0; i < nParticles; i++) {
            // Random Position in [-2.5, 2.5] box
            posArray[i*3+0] = (Math.random() - 0.5) * 5.0;
            posArray[i*3+1] = (Math.random() - 0.5) * 5.0;
            posArray[i*3+2] = (Math.random() - 0.5) * 5.0;
            
            // Random Rotation (Random axis + Random angle)
            const axis = new Float32Array([Math.random()-0.5, Math.random()-0.5, Math.random()-0.5]);
            const len = Math.sqrt(axis[0]*axis[0] + axis[1]*axis[1] + axis[2]*axis[2]);
            if(len > 0) { axis[0]/=len; axis[1]/=len; axis[2]/=len; }
            const angle = Math.random() * Math.PI * 2;
            const c = Math.cos(angle/2);
            const s = Math.sin(angle/2);
            
            quatArray[i*4+0] = axis[0] * s;
            quatArray[i*4+1] = axis[1] * s;
            quatArray[i*4+2] = axis[2] * s;
            quatArray[i*4+3] = c;
        }

        console.log("Positions:", posArray);
        console.log("Quaternions:", quatArray);

        // 4. Update Renderer
        this.renderer.updateInstances(posArray, quatArray);

        // 5. Start Animation Loop
        this.animate();
    }

    animate() {
        requestAnimationFrame(() => this.animate());
        this.renderer.render();
    }
}

new App();
