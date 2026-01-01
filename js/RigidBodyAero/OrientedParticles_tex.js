import * as THREE from 'three';

export class OrientedParticles_tex {
    constructor(scene, camera) {
        this.scene = scene;
        this.camera = camera;
        this.mesh = null;
        this.nParticles = 0;
        this.size = 0;
    }

    async init(nParticles, size, bodyDef) {
        this.nParticles = nParticles;
        this.size = size;

        const vertSrc = await fetch('./shaders/orientedParticle_texture.glslv').then(r => r.text());
        const fragSrc = await fetch('./shaders/orientedParticle.glslf').then(r => r.text());

        const geometry = this.createGeometry(bodyDef);
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
        const newGeom = this.createGeometry(bodyDef);
        newGeom.instanceCount = this.nParticles;
        this.mesh.geometry.dispose();
        this.mesh.geometry = newGeom;
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
