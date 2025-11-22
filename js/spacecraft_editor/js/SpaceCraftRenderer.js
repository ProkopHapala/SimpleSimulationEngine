
class SpaceCraftRenderer extends MeshRenderer {
    constructor(engine) {
        super(null, null, 0);
        this.engine = engine;
        this.labelMode = 'none';
    }

    init(container, shaders) {
        this.scene = new THREE.Scene();
        this.scene.background = new THREE.Color(0x111111);

        // Camera
        const aspect = container.clientWidth / container.clientHeight;
        this.camera = new THREE.PerspectiveCamera(60, aspect, 0.1, 1000);
        this.camera.position.set(0, 0, 20);

        // Renderer
        this.renderer = new THREE.WebGLRenderer({ antialias: true });
        this.renderer.setSize(container.clientWidth, container.clientHeight);
        container.appendChild(this.renderer.domElement);

        // Controls
        this.controls = new THREE.OrbitControls(this.camera, this.renderer.domElement);

        // Lights
        const ambientLight = new THREE.AmbientLight(0x404040);
        this.scene.add(ambientLight);
        const dirLight = new THREE.DirectionalLight(0xffffff, 0.8);
        dirLight.position.set(1, 1, 1);
        this.scene.add(dirLight);

        // MeshRenderer Init
        this.shaders = shaders;
        this.capacity = 10000; // Max verts

        MeshRenderer.prototype.init.call(this);

        // 3. Grid & Axes (from original SpaceCraftRenderer)
        this.scene.add(new THREE.GridHelper(20, 20));
        this.scene.add(new THREE.AxesHelper(2));

        // Handle Resize
        window.addEventListener('resize', () => this.onWindowResize());
    }

    updateGeometry(meshBuilder) {
        if (!meshBuilder) return;

        const verts = meshBuilder.verts;
        const numVerts = verts.length;

        if (numVerts > this.capacity) {
            console.warn("Mesh exceeds max vertices capacity!");
            return;
        }

        // 1. Update Positions
        for (let i = 0; i < numVerts; i++) {
            const v = verts[i];
            this.posData[i * 4] = v.pos.x;
            this.posData[i * 4 + 1] = v.pos.y;
            this.posData[i * 4 + 2] = v.pos.z;
            this.posData[i * 4 + 3] = 1.0;
        }
        this.posTexture.needsUpdate = true;

        // 2. Update Particles (Nodes)
        this.updateParticles(numVerts, (i) => {
            return [0.5, 0.5, 1.0]; // Blue-ish
        }, (i) => {
            return 0.2; // Fixed size for nodes
        });

        // 3. Update Bonds (Girders)
        const edges = meshBuilder.edges;
        const pairs = [];
        for (let i = 0; i < edges.length; i++) {
            pairs.push([edges[i].x, edges[i].y]);
        }
        this.updateBonds(pairs);

        // 4. Update Labels
        this.updateLabelsContent();
    }

    updateLabelsContent() {
        if (this.labelMode === 'none' || !this.labelMesh) return;

        const numVerts = this.atomMesh.count; // Or track separately

        this.updateLabels((i) => {
            if (this.labelMode === 'id') {
                return i.toString();
            }
            return "";
        }, numVerts);
    }

    setLabelMode(mode) {
        if (this.labelMode !== mode) {
            this.labelMode = mode;
            if (this.labelMesh) {
                this.labelMesh.visible = (mode !== 'none');
                if (this.labelMesh.visible) {
                    this.updateLabelsContent();
                }
            }
        }
    }

    update() {
        if (this.controls) this.controls.update();
    }

    render() {
        this.renderer.render(this.scene, this.camera);
    }

    onWindowResize() {
        if (this.camera && this.renderer && this.renderer.domElement.parentElement) {
            const container = this.renderer.domElement.parentElement;
            this.camera.aspect = container.clientWidth / container.clientHeight;
            this.camera.updateProjectionMatrix();
            this.renderer.setSize(container.clientWidth, container.clientHeight);
        }
    }
}
