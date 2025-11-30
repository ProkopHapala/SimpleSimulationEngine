import { MeshRenderer } from '../../common_js/MeshRenderer.js';

export class SpaceCraftRenderer extends MeshRenderer {
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
        const height = 20;
        const width = height * aspect;
        this.camera = new THREE.OrthographicCamera(-width / 2, width / 2, height / 2, -height / 2, -1000, 1000);
        this.camera.position.set(0, 0, 20);
        this.cameraMode = 'ortho';

        // Renderer
        this.renderer = new THREE.WebGLRenderer({ antialias: true });
        this.renderer.setSize(container.clientWidth, container.clientHeight);
        container.appendChild(this.renderer.domElement);

        // Controls (basic OrbitControls; detailed mouse/key bindings are configured by GUI)
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
        this.gridHelper = new THREE.GridHelper(20, 20);
        this.scene.add(this.gridHelper);
        
        this.axesHelper = new THREE.AxesHelper(2);
        this.scene.add(this.axesHelper);

        // Handle Resize
        window.addEventListener('resize', () => this.onWindowResize());
    }

    setGridVisible(visible) {
        if (this.gridHelper) this.gridHelper.visible = visible;
    }

    setAxisVisible(visible) {
        if (this.axesHelper) this.axesHelper.visible = visible;
    }

    updateGeometry(meshBuilder) {
        if (!meshBuilder) return;

        const verts = meshBuilder.verts;
        const numVerts = verts.length;

        if (numVerts > this.capacity) {
            if (typeof logger !== 'undefined') logger.warn("Mesh exceeds max vertices capacity!");
            else console.warn("Mesh exceeds max vertices capacity!");
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
            const e = edges[i];
            // Pass edge.z as material ID so post-processed rails (e.g. from addSkipEdges)
            // can be colored differently by MeshRenderer.
            pairs.push([e.x, e.y, e.z]);
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

    setLabelScreenSpace(enabled) {
        this.setLabelStyle(null, null, enabled);
        // If switching to screen space, we might want to adjust default scale?
        // User controls scale via GUI, so let them adjust.
    }

    setCameraMode(mode) {
        if (this.cameraMode === mode) return;
        this.cameraMode = mode;

        const oldCam = this.camera;
        const aspect = oldCam.aspect;
        const fov = 60;
        const near = 0.1;
        const far = 1000;

        // Save state
        const pos = oldCam.position.clone();
        const target = this.controls.target.clone();
        const zoom = this.controls.object.zoom; // For ortho

        if (mode === 'ortho') {
            // Estimate frustum size based on distance to target
            const dist = pos.distanceTo(target);
            const height = 2 * Math.tan((fov * Math.PI / 180) / 2) * dist;
            const width = height * aspect;

            this.camera = new THREE.OrthographicCamera(-width / 2, width / 2, height / 2, -height / 2, near, far);
            this.camera.zoom = 1; // Reset zoom for ortho, as width/height already capture view
        } else {
            this.camera = new THREE.PerspectiveCamera(fov, aspect, near, far);
        }

        this.camera.position.copy(pos);
        this.camera.up.copy(oldCam.up);
        this.camera.lookAt(target);

        // Update controls
        this.controls.object = this.camera;
        this.controls.update();

        // Update aspect uniform
        this.updateLabelUniforms(aspect);
    }

    update() {
        if (this.controls) this.controls.update();
        // Ensure aspect is up to date if needed? No, onWindowResize handles it.
    }

    render() {
        this.renderer.render(this.scene, this.camera);
    }

    onWindowResize() {
        if (this.camera && this.renderer && this.renderer.domElement.parentElement) {
            const container = this.renderer.domElement.parentElement;
            const width = container.clientWidth;
            const height = container.clientHeight;
            const aspect = width / height;

            if (this.camera.isPerspectiveCamera) {
                this.camera.aspect = aspect;
            } else {
                // Ortho: Maintain vertical size, adjust horizontal
                // Or maintain scale?
                // Let's keep the height constant in world units, adjust width
                const camH = (this.camera.top - this.camera.bottom) / this.camera.zoom;
                const camW = camH * aspect;
                this.camera.left = -camW / 2;
                this.camera.right = camW / 2;
                this.camera.top = camH / 2;
                this.camera.bottom = -camH / 2;
            }

            this.camera.updateProjectionMatrix();
            this.renderer.setSize(width, height);

            this.updateLabelUniforms(aspect);
        }
    }
}

export default SpaceCraftRenderer;
