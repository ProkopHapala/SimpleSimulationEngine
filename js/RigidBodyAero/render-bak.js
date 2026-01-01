import * as THREE from 'three';

export class Renderer {
    constructor(container, nParticles, size) {
        this.container = container;
        this.nParticles = nParticles;
        this.size = size;

        this.scene = new THREE.Scene();
        this.camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);
        this.camera.position.set(0, 30, 60);
        this.camera.lookAt(0, 15, 0);

        this.renderer = new THREE.WebGLRenderer({ antialias: true });
        this.renderer.setSize(window.innerWidth, window.innerHeight);
        this.renderer.setClearColor(0x111111);
        this.container.appendChild(this.renderer.domElement);

        this.initLights();
        this.initInstancedMesh();
        
        window.addEventListener('resize', () => this.onWindowResize());
    }

    initLights() {
        const ambientLight = new THREE.AmbientLight(0xffffff, 0.5);
        this.scene.add(ambientLight);
        const directionalLight = new THREE.DirectionalLight(0xffffff, 0.8);
        directionalLight.position.set(1, 1, 1).normalize();
        this.scene.add(directionalLight);
    }

    initInstancedMesh() {
        // Create a dart shape: 3 triangles (Inverted T)
        // 1. Wing-platform (horizontal triangle)
        // 2. Rudder (vertical triangle)
        const geometry = new THREE.BufferGeometry();
        const vertices = new Float32Array([
            // Wing (horizontal)
            -0.5, 0, 0.5,  0.5, 0, 0.5,  0, 0, -0.5,
            // Rudder (vertical)
            0, 0, 0.5,  0, 0.5, 0.5,  0, 0, -0.5
        ]);
        geometry.setAttribute('position', new THREE.BufferAttribute(vertices, 3));
        
        this.material = new THREE.ShaderMaterial({
            uniforms: {
                u_tex_pos: { value: null },
                u_tex_quat: { value: null },
                u_size: { value: this.size },
                u_static_test: { value: 0.0 }
            },
            glslVersion: THREE.GLSL3,
            vertexShader: `
                uniform sampler2D u_tex_pos;
                uniform sampler2D u_tex_quat;
                uniform float u_size;
                uniform float u_static_test;
                
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
                    float id = float(gl_InstanceID);
                    vec2 uv = vec2(mod(id, u_size) + 0.5, floor(id / u_size) + 0.5) / u_size;
                    
                    vec4 pos_data = texture(u_tex_pos, uv);
                    vec4 quat = texture(u_tex_quat, uv);
                    
                    if( u_static_test > 0.5 ){
                        pos_data = vec4(0.0, 15.0, 0.0, 1.0);
                        quat     = vec4(0.0, 0.0, 0.0, 1.0);
                    }

                    mat3 R = quat_to_mat3(quat);
                    vec3 transformed = R * position + pos_data.xyz;
                    
                    gl_Position = projectionMatrix * modelViewMatrix * vec4(transformed, 1.0);
                }
            `,
            fragmentShader: `
                out vec4 pc_fragColor;
                void main() {
                    pc_fragColor = vec4(1.0, 0.5, 0.2, 1.0);
                }
            `,
            side: THREE.DoubleSide
        });

        this.mesh = new THREE.InstancedMesh(geometry, this.material, this.nParticles);
        this.mesh.frustumCulled = false;
        this.scene.add(this.mesh);
    }

    onWindowResize() {
        this.camera.aspect = window.innerWidth / window.innerHeight;
        this.camera.updateProjectionMatrix();
        this.renderer.setSize(window.innerWidth, window.innerHeight);
    }

    render(posTex, quatTex) {
        if (!this.threePosTex) {
            this.threePosTex = new THREE.Texture();
            this.threeQuatTex = new THREE.Texture();
            
            // Hack to wrap raw WebGLTexture
            const propsP = this.renderer.properties.get(this.threePosTex);
            propsP.__webglTexture = posTex;
            propsP.__webglInit = true;
            
            const propsQ = this.renderer.properties.get(this.threeQuatTex);
            propsQ.__webglTexture = quatTex;
            propsQ.__webglInit = true;

            this.material.uniforms.u_tex_pos.value = this.threePosTex;
            this.material.uniforms.u_tex_quat.value = this.threeQuatTex;
        } else {
            // Update the underlying WebGLTexture in case it changed (ping-pong)
            this.renderer.properties.get(this.threePosTex).__webglTexture = posTex;
            this.renderer.properties.get(this.threeQuatTex).__webglTexture = quatTex;
        }
        
        this.renderer.setViewport(0, 0, window.innerWidth, window.innerHeight);
        this.renderer.render(this.scene, this.camera);
    }
}
