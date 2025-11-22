
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
        window.addEventListener('resize', () => {
            this.camera.aspect = container.clientWidth / container.clientHeight;
            this.camera.updateProjectionMatrix();
            this.renderer.setSize(container.clientWidth, container.clientHeight);
        });
    }

    initDataTexture() {
        // Create DataTexture from Engine's Float32Array
        // Format: RGBA (4 floats per pixel)
        // Type: FloatType
        this.dataTexture = new THREE.DataTexture(
            this.engine.nodes,
            this.engine.maxNodes,
            1,
            THREE.RGBAFormat,
            THREE.FloatType
        );
        this.dataTexture.minFilter = THREE.NearestFilter;
        this.dataTexture.magFilter = THREE.NearestFilter;
        this.dataTexture.needsUpdate = true;

        this.engine.dataTexture = this.dataTexture; // Link back
    }

    initNodeMesh() {
        // Geometry: Simple Plane (Billboard)
        const geometry = new THREE.PlaneBufferGeometry(1, 1);

        // Instanced Attributes
        // We need to manually allocate max instances
        const maxInstances = this.engine.maxNodes;
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
                uTexSize: { value: new THREE.Vector2(this.engine.maxNodes, 1) },
                uScale: { value: 0.2 }
            },
            side: THREE.DoubleSide,
            transparent: true
        });

        this.nodeMesh = new THREE.Mesh(instancedGeometry, material);
        this.nodeMesh.frustumCulled = false; // Important: Bounds are dynamic
        this.scene.add(this.nodeMesh);
    }

    initGirderMesh() {
        // Geometry: Cylinder (Y-up)
        // Radius 0.5, Height 1.0
        const geometry = new THREE.CylinderBufferGeometry(0.5, 0.5, 1.0, 8, 1, true);
        geometry.translate(0, 0.5, 0); // Pivot at bottom? No, shader handles center. Let's keep center at 0.
        // Actually shader expects local Y to be length.
        // Let's use a unit cylinder centered at 0.

        const maxInstances = 4096; // Max girders
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
                uTexSize: { value: new THREE.Vector2(this.engine.maxNodes, 1) },
                uThickness: { value: 0.05 }
            },
            side: THREE.DoubleSide
        });

        this.girderMesh = new THREE.Mesh(instancedGeometry, material);
        this.girderMesh.frustumCulled = false;
        this.scene.add(this.girderMesh);
    }

    updateGeometry() {
        // 1. Update Texture
        this.dataTexture.needsUpdate = true;

        // 2. Update Node Count
        this.nodeMesh.geometry.instanceCount = this.engine.nodeCount;

        // 3. Update Girders
        const girderCount = this.engine.girders.length;
        const ids = this.girderMesh.geometry.attributes.aNodeIDs.array;

        for (let i = 0; i < girderCount; i++) {
            const g = this.engine.girders[i];
            ids[i * 2] = g.n1;
            ids[i * 2 + 1] = g.n2;
        }

        this.girderMesh.geometry.attributes.aNodeIDs.needsUpdate = true;
        this.girderMesh.geometry.instanceCount = girderCount;
    }

    update() {
        this.controls.update();
    }

    render() {
        this.renderer.render(this.scene, this.camera);
    }
}
