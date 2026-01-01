import * as THREE from 'three';

export class OrientedParticles_tex {
    constructor(scene, camera) {
        this.scene = scene;
        this.camera = camera;
        this.mesh = null;
        this.nParticles = 0;
        this.texSize = 0;
    }

    init(nParticles, texSize) {
        this.nParticles = nParticles;
        this.texSize = texSize;

        // Dart Geometry
        const positions = [
            -0.5, 0, 0,   0, 0, 0,   0, 0, 1.0,  // Left Wing
             0.5, 0, 0,   0, 0, 0,   0, 0, 1.0,  // Right Wing
             0, 0.5, 0,   0, 0, 0,   0, 0, 1.0   // Rudder
        ];
        const colors = [
            1, 0, 0,  0, 0, 0,  1, 1, 1, // Red
            0, 1, 0,  0, 0, 0,  1, 1, 1, // Green
            0, 0, 1,  0, 0, 0,  1, 1, 1  // Blue
        ];

        const baseGeometry = new THREE.BufferGeometry();
        baseGeometry.setAttribute('position', new THREE.Float32BufferAttribute(positions, 3));
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
                u_tex_size: { value: texSize }
            },
            vertexShader: `
                uniform sampler2D u_tex_pos;
                uniform sampler2D u_tex_quat;
                uniform float u_tex_size;
                varying vec3 vColor;

                mat3 quat_to_mat3(vec4 q) {
                    float x2 = q.x + q.x; float y2 = q.y + q.y; float z2 = q.z + q.z;
                    float xx = q.x * x2;  float xy = q.x * y2;  float xz = q.x * z2;
                    float yy = q.y * y2;  float yz = q.y * z2;  float zz = q.z * z2;
                    float wx = q.w * x2;  float wy = q.w * y2;  float wz = q.w * z2;
                    return mat3(
                        1.0 - (yy + zz), xy + wz,          xz - wy,
                        xy - wz,         1.0 - (xx + zz),  yz + wx,
                        xz + wy,         yz - wx,          1.0 - (xx + yy)
                    );
                }

                void main() {
                    vColor = color;
                    float idx = float(gl_InstanceID);
                    vec2 uv = vec2(
                        mod(idx, u_tex_size) + 0.5,
                        floor(idx / u_tex_size) + 0.5
                    ) / u_tex_size;

                    vec4 p_data = texture2D(u_tex_pos, uv);
                    vec4 q_data = texture2D(u_tex_quat, uv);

                    vec3 pos = p_data.xyz;
                    mat3 R = quat_to_mat3(q_data);
                    vec3 transformed = R * position + pos;
                    
                    gl_Position = projectionMatrix * modelViewMatrix * vec4(transformed, 1.0);
                }
            `,
            fragmentShader: `
                varying vec3 vColor;
                void main() {
                    gl_FragColor = vec4(vColor, 1.0);
                }
            `
        });

        this.mesh = new THREE.Mesh(geometry, this.material);
        this.mesh.frustumCulled = false;
        this.scene.add(this.mesh);
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
