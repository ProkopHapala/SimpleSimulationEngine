import * as THREE from 'three';

export class OrientedParticles_tex {
    constructor(scene, camera) {
        this.scene = scene;
        this.camera = camera;
        this.mesh = null;
        this.nParticles = 0;
        this.size = 0;
    }

    async init(nParticles, size, comOffset = 0.5) {
        this.nParticles = nParticles;
        this.size = size;

        // Correct shader path (texture-driven vertex shader)
        const vertSrc = await fetch('./shaders/orientedParticle_texture.glslv').then(r => r.text());
        const fragSrc = await fetch('./shaders/orientedParticle.glslf').then(r => r.text());

        this.basePositions = this.computePositions(comOffset);
        const colors = [
            1, 0, 0,  0, 0, 0,  1, 1, 1, // Red
            0, 1, 0,  0, 0, 0,  1, 1, 1, // Green
            0, 0, 1,  0, 0, 0,  1, 1, 1  // Blue
        ];

        const baseGeometry = new THREE.BufferGeometry();
        baseGeometry.setAttribute('position', new THREE.Float32BufferAttribute(this.basePositions, 3));
        baseGeometry.setAttribute('color', new THREE.Float32BufferAttribute(colors, 3));

        const geometry = new THREE.InstancedBufferGeometry();
        geometry.copy(baseGeometry);
        geometry.instanceCount = nParticles;

        this.material = new THREE.ShaderMaterial({
            vertexColors: true,
            side: THREE.DoubleSide,
            uniforms: {
                u_tex_pos: { value: null },
                u_tex_quat: { value: null },
                u_tex_size: { value: size }
            },
            vertexShader: vertSrc,
            fragmentShader: fragSrc
        });

        this.mesh = new THREE.Mesh(geometry, this.material);
        this.mesh.frustumCulled = false;
        this.scene.add(this.mesh);
    }

    computePositions(comOffset = 0.5) {
        const z_shift = -comOffset;
        return [
            -0.5, 0, 0 + z_shift,   0, 0, 0 + z_shift,   0, 0, 1.0 + z_shift,
             0.5, 0, 0 + z_shift,   0, 0, 0 + z_shift,   0, 0, 1.0 + z_shift,
             0, 0.5, 0 + z_shift,   0, 0, 0 + z_shift,   0, 0, 1.0 + z_shift
        ];
    }

    updateGeometry(comOffset = 0.5) {
        if (!this.mesh || !this.mesh.geometry) return;
        this.basePositions = this.computePositions(comOffset);
        this.mesh.geometry.setAttribute('position', new THREE.Float32BufferAttribute(this.basePositions, 3));
        this.mesh.geometry.attributes.position.needsUpdate = true;
    }

    updateTextures(texPos, texQuat) {
        if (this.material) {
            this.material.uniforms.u_tex_pos.value = texPos;
            this.material.uniforms.u_tex_quat.value = texQuat;
        }
    }

    setVisibility(visible) {
        if (this.mesh) this.mesh.visible = visible;
    }
}
