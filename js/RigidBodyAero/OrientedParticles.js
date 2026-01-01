import * as THREE from 'three';

export class Renderer {
    constructor(container) {
        this.container = container;
        
        // 1. Scene Setup
        this.scene = new THREE.Scene();
        this.scene.background = new THREE.Color(0x222222); // Dark grey background

        // 2. Camera Setup
        const aspect = window.innerWidth / window.innerHeight;
        const d = 20; // Orthographic frustum vertical half-size
        this.camera = new THREE.OrthographicCamera(-d * aspect, d * aspect, d, -d, 0.1, 1000);
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

    async initInstancedMesh(nParticles, bodyDef) {
        this.nParticles = nParticles;
        const vertSrc = await fetch('./shaders/orientedParticle_buffer.glslv').then(r => r.text());
        const fragSrc = await fetch('./shaders/orientedParticle.glslf').then(r => r.text());

        const geometry = this.createGeometry(bodyDef);
        geometry.instanceCount = nParticles;

        // 3. Instanced Attributes
        geometry.setAttribute('a_instance_pos', new THREE.InstancedBufferAttribute(new Float32Array(nParticles * 3), 3));
        geometry.setAttribute('a_instance_quat', new THREE.InstancedBufferAttribute(new Float32Array(nParticles * 4), 4));

        // 4. Shader Material
        this.material = new THREE.ShaderMaterial({
            vertexColors: true,
            side: THREE.DoubleSide,
            uniforms: {},
            vertexShader: vertSrc,
            fragmentShader: fragSrc
        });

        // 5. Mesh
        this.mesh = new THREE.Mesh(geometry, this.material);
        this.instancedGeometry = geometry;
        this.mesh.frustumCulled = false; // Disable culling for safety
        this.scene.add(this.mesh);
    }

    createGeometry(bodyDef) {
        const positions = [];
        const colors = [];
        const wingWidth = 0.5;
        const wingLength = 0.5;

        bodyDef.points.forEach((p, idx) => {
            const side = new THREE.Vector3().crossVectors(p.t, p.n).normalize();
            const v0 = p.pos.clone().addScaledVector(side, -wingWidth);
            const v1 = p.pos.clone().addScaledVector(side, wingWidth);
            const v2 = p.pos.clone().addScaledVector(p.t, -wingLength);

            positions.push(v0.x, v0.y, v0.z, v1.x, v1.y, v1.z, v2.x, v2.y, v2.z);
            
            const color = new THREE.Color().setHSL((idx * 0.3) % 1.0, 0.8, 0.5);
            colors.push(color.r, color.g, color.b, color.r, color.g, color.b, color.r, color.g, color.b);
        });

        const baseGeometry = new THREE.BufferGeometry();
        baseGeometry.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
        baseGeometry.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));

        const geometry = new THREE.InstancedBufferGeometry();
        geometry.copy(baseGeometry);
        return geometry;
    }

    updateGeometry(bodyDef) {
        if (!this.mesh) return;
        
        // Preserve instanced attributes
        const oldPosAttr = this.instancedGeometry.getAttribute('a_instance_pos');
        const oldQuatAttr = this.instancedGeometry.getAttribute('a_instance_quat');
        
        const newGeom = this.createGeometry(bodyDef);
        newGeom.instanceCount = this.nParticles;
        newGeom.setAttribute('a_instance_pos', oldPosAttr);
        newGeom.setAttribute('a_instance_quat', oldQuatAttr);
        
        this.mesh.geometry.dispose();
        this.mesh.geometry = newGeom;
        this.instancedGeometry = newGeom;
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
        const aspect = window.innerWidth / window.innerHeight;
        const d = 20;
        this.camera.left = -d * aspect;
        this.camera.right = d * aspect;
        this.camera.top = d;
        this.camera.bottom = -d;
        this.camera.updateProjectionMatrix();
        this.renderer.setSize(window.innerWidth, window.innerHeight);
    }
}
