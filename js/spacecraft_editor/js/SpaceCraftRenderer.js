
class SpaceCraftRenderer {
    constructor(engine) {
        this.engine = engine;
        this.scene = null;
        this.camera = null;
        this.renderer = null;
        this.controls = null;

        // Meshes
        this.nodeMesh = null;
        this.girderMesh = null;
    }

    init(container) {
        // 1. Setup Three.js
        this.scene = new THREE.Scene();
        this.scene.background = new THREE.Color(0x222222);

        this.camera = new THREE.PerspectiveCamera(60, container.clientWidth / container.clientHeight, 0.1, 1000);
        this.camera.position.set(5, 5, 10);

        this.renderer = new THREE.WebGLRenderer({ antialias: true });
        this.renderer.setSize(container.clientWidth, container.clientHeight);
        container.appendChild(this.renderer.domElement);

        // 2. Controls
        this.controls = new THREE.OrbitControls(this.camera, this.renderer.domElement);
        this.controls.enableDamping = false; // Disabled as per user request for snappier control

        // 3. Grid & Axes
        this.scene.add(new THREE.GridHelper(20, 20));
        this.scene.add(new THREE.AxesHelper(2));

        // 4. Init GPU Data
        this.initDataTexture();

        // 5. Init Meshes (Empty initially)
        this.initNodeMesh();
        this.initGirderMesh();

        // Handle Resize
        window.addEventListener('resize', () => this.resize());
    }

    resize() {
        const container = this.renderer.domElement.parentElement;
        if (!container) return;

        this.camera.aspect = container.clientWidth / container.clientHeight;
        this.camera.updateProjectionMatrix();
        this.renderer.setSize(container.clientWidth, container.clientHeight);
    }

    initDataTexture() {
        // Create DataTexture from Engine's MeshBuilder Verts
        // Format: RGBA (4 floats per pixel)
        // Type: FloatType
        // Max size: 4096 verts for now
        this.maxVerts = 4096;
        this.posData = new Float32Array(this.maxVerts * 4);

        this.dataTexture = new THREE.DataTexture(
            this.posData,
            this.maxVerts,
            1,
            THREE.RGBAFormat,
            THREE.FloatType
        );
        this.dataTexture.minFilter = THREE.NearestFilter;
        this.dataTexture.magFilter = THREE.NearestFilter;
        this.dataTexture.needsUpdate = true;
    }

    initNodeMesh() {
        // Geometry: Simple Plane (Billboard)
        const geometry = new THREE.PlaneBufferGeometry(1, 1);

        // Instanced Attributes
        const maxInstances = this.maxVerts;
        const instancedGeometry = new THREE.InstancedBufferGeometry();
        instancedGeometry.index = geometry.index;
        instancedGeometry.attributes = geometry.attributes;

        // aNodeID attribute (0, 1, 2...)
        const ids = new Float32Array(maxInstances);
        for (let i = 0; i < maxInstances; i++) ids[i] = i;
        instancedGeometry.setAttribute('aNodeID', new THREE.InstancedBufferAttribute(ids, 1));

        // Material
        const material = new THREE.RawShaderMaterial({
            vertexShader: Shaders.nodeVertex,
            fragmentShader: Shaders.nodeFragment,
            uniforms: {
                uPosTex: { value: this.dataTexture },
                uTexSize: { value: new THREE.Vector2(this.maxVerts, 1) },
                uScale: { value: 0.2 },
                // Standard uniforms
                modelViewMatrix: { value: new THREE.Matrix4() },
                projectionMatrix: { value: new THREE.Matrix4() }
            },
            side: THREE.DoubleSide,
            transparent: true
        });

        // Hook for standard uniforms update
        material.onBeforeCompile = (shader) => {
            // This is RawShaderMaterial, so we must manage uniforms manually or use ShaderMaterial
            // But we are using RawShaderMaterial, so we need to pass standard uniforms manually in render loop?
            // Actually, Three.js does NOT auto-update uniforms for RawShaderMaterial unless we use `uniformsLib` or similar.
            // For simplicity, let's assume we update them in `render()` or switch to `ShaderMaterial`.
            // Let's stick to Raw but ensure we pass matrices.
        };

        this.nodeMesh = new THREE.Mesh(instancedGeometry, material);
        this.nodeMesh.frustumCulled = false;
        this.scene.add(this.nodeMesh);
    }

    initGirderMesh() {
        // Geometry: Cylinder (Y-up)
        const geometry = new THREE.CylinderBufferGeometry(0.5, 0.5, 1.0, 8, 1, true);
        geometry.translate(0, 0.5, 0); // Pivot at base? No, shader handles it.

        const maxInstances = 4096; // Max edges
        const instancedGeometry = new THREE.InstancedBufferGeometry();
        instancedGeometry.index = geometry.index;
        instancedGeometry.attributes = geometry.attributes;

        // Attributes to be updated
        this.girderIDs = new Float32Array(maxInstances * 2); // [idA, idB]
        instancedGeometry.setAttribute('aNodeIDs', new THREE.InstancedBufferAttribute(this.girderIDs, 2));

        const material = new THREE.RawShaderMaterial({
            vertexShader: Shaders.girderVertex,
            fragmentShader: Shaders.girderFragment,
            uniforms: {
                uPosTex: { value: this.dataTexture },
                uTexSize: { value: new THREE.Vector2(this.maxVerts, 1) },
                uThickness: { value: 0.05 },
                modelViewMatrix: { value: new THREE.Matrix4() },
                projectionMatrix: { value: new THREE.Matrix4() },
                viewMatrix: { value: new THREE.Matrix4() }
            },
            side: THREE.DoubleSide
        });

        this.girderMesh = new THREE.Mesh(instancedGeometry, material);
        this.girderMesh.frustumCulled = false;
        this.scene.add(this.girderMesh);
    }

    updateGeometry(meshBuilder) {
        if (!meshBuilder) return;

        // 1. Update Texture (Vertices)
        const verts = meshBuilder.verts; // Flat array [x,y,z, x,y,z...]
        const numVerts = verts.length / 3;

        for (let i = 0; i < numVerts; i++) {
            this.posData[i * 4] = verts[i * 3];
            this.posData[i * 4 + 1] = verts[i * 3 + 1];
            this.posData[i * 4 + 2] = verts[i * 3 + 2];
            this.posData[i * 4 + 3] = 0; // Type
        }
        this.dataTexture.needsUpdate = true;

        // 2. Update Node Count
        this.nodeMesh.geometry.instanceCount = numVerts;

        // 3. Update Girders (Edges)
        const edges = meshBuilder.edges; // Flat array [a,b,type,0...]
        const numEdges = edges.length / 4;
        const ids = this.girderMesh.geometry.attributes.aNodeIDs.array;

        for (let i = 0; i < numEdges; i++) {
            ids[i * 2] = edges[i * 4];     // Node A index
            ids[i * 2 + 1] = edges[i * 4 + 1]; // Node B index
        }

        this.girderMesh.geometry.attributes.aNodeIDs.needsUpdate = true;
        this.girderMesh.geometry.instanceCount = numEdges;
    }

    update() {
        this.controls.update();

        // Manually update uniforms for RawShaderMaterial if needed
        // (Three.js usually handles modelViewMatrix and projectionMatrix for Mesh,
        // but for RawShaderMaterial we might need to be careful.
        // However, since we passed them as uniforms, Three.js renderer logic *should* populate them
        // if we don't override them with static values.
        // Actually, for RawShaderMaterial, Three.js DOES NOT populate standard uniforms automatically
        // unless we use `uniforms: THREE.UniformsUtils.merge([THREE.UniformsLib.common...])` AND set `lights: true` etc.
        // OR we manually update them here.)

        if (this.nodeMesh) {
            this.nodeMesh.material.uniforms.modelViewMatrix.value.copy(this.nodeMesh.modelViewMatrix);
            this.nodeMesh.material.uniforms.projectionMatrix.value.copy(this.camera.projectionMatrix);
        }
        if (this.girderMesh) {
            this.girderMesh.material.uniforms.modelViewMatrix.value.copy(this.girderMesh.modelViewMatrix);
            this.girderMesh.material.uniforms.projectionMatrix.value.copy(this.camera.projectionMatrix);
            this.girderMesh.material.uniforms.viewMatrix.value.copy(this.camera.matrixWorldInverse);
        }
    }

    render() {
        this.renderer.render(this.scene, this.camera);
    }
}
