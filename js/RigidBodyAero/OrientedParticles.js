import * as THREE from 'three';

export class Renderer {
    constructor(container) {
        this.container = container;
        
        // 1. Scene Setup
        this.scene = new THREE.Scene();
        this.scene.background = new THREE.Color(0x222222); // Dark grey background

        // 2. Camera Setup
        this.camera = new THREE.PerspectiveCamera(60, window.innerWidth / window.innerHeight, 0.1, 1000);
        this.camera.position.set(20, 20, 20);
        this.camera.lookAt(0, 0, 0);

        // 3. WebGL2 Renderer Setup
        const canvas = document.createElement('canvas');
        const context = canvas.getContext('webgl2', { alpha: false, antialias: true });
        
        this.renderer = new THREE.WebGLRenderer({ canvas: canvas, context: context });
        this.renderer.setSize(window.innerWidth, window.innerHeight);
        this.container.appendChild(this.renderer.domElement);

        // 4. Static Helpers (Axes & Grid)
        this.axesHelper = new THREE.AxesHelper(10);
        this.scene.add(this.axesHelper);

        this.gridHelper = new THREE.GridHelper(20, 20);
        this.scene.add(this.gridHelper);

        // 5. Lighting (Good practice to have, even if only for helpers initially)
        const ambientLight = new THREE.AmbientLight(0xffffff, 0.6);
        this.scene.add(ambientLight);
    }

    async initInstancedMesh(nParticles) {
        this.nParticles = nParticles;
        const vertSrc = await fetch('./shaders/orientedParticle_buffer.glslv').then(r => r.text());
        const fragSrc = await fetch('./shaders/orientedParticle.glslf').then(r => r.text());

        // 1. Geometry (Reuse Step 1.5 Geometry)
        const positions = [
            // Left Wing (Left Tip, Back Center, Nose)
            -0.5, 0, 0,   0, 0, 0,   0, 0, 1.0,
            // Right Wing (Right Tip, Back Center, Nose)
            0.5, 0, 0,    0, 0, 0,   0, 0, 1.0,
            // Rudder (Top Tip, Back Center, Nose)
            0, 0.5, 0,    0, 0, 0,   0, 0, 1.0
        ];

        const colors = [
            // Left Wing Colors (Red, Black, White)
            1, 0, 0,      0, 0, 0,   1, 1, 1,
            // Right Wing Colors (Green, Black, White)
            0, 1, 0,      0, 0, 0,   1, 1, 1,
            // Rudder Colors (Blue, Black, White)
            0, 0, 1,      0, 0, 0,   1, 1, 1
        ];

        const baseGeometry = new THREE.BufferGeometry();
        baseGeometry.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
        baseGeometry.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));

        // 2. Instanced Geometry
        this.instancedGeometry = new THREE.InstancedBufferGeometry();
        this.instancedGeometry.copy(baseGeometry);
        this.instancedGeometry.instanceCount = nParticles;

        // 3. Instanced Attributes
        this.instancedGeometry.setAttribute('a_instance_pos', new THREE.InstancedBufferAttribute(new Float32Array(nParticles * 3), 3));
        this.instancedGeometry.setAttribute('a_instance_quat', new THREE.InstancedBufferAttribute(new Float32Array(nParticles * 4), 4));

        // 4. Shader Material
        this.material = new THREE.ShaderMaterial({
            vertexColors: true,
            side: THREE.DoubleSide,
            uniforms: {},
            vertexShader: vertSrc,
            fragmentShader: fragSrc
        });

        // 5. Mesh
        this.mesh = new THREE.Mesh(this.instancedGeometry, this.material);
        this.mesh.frustumCulled = false; // Disable culling for safety
        this.scene.add(this.mesh);
    }

    updateInstances(posArray, quatArray) {
        if (this.instancedGeometry) {
            this.instancedGeometry.attributes.a_instance_pos.array.set(posArray);
            this.instancedGeometry.attributes.a_instance_pos.needsUpdate = true;
            
            this.instancedGeometry.attributes.a_instance_quat.array.set(quatArray);
            this.instancedGeometry.attributes.a_instance_quat.needsUpdate = true;
        }
    }

    render() {
        this.renderer.render(this.scene, this.camera);
    }

    onWindowResize() {
        this.camera.aspect = window.innerWidth / window.innerHeight;
        this.camera.updateProjectionMatrix();
        this.renderer.setSize(window.innerWidth, window.innerHeight);
    }
}
