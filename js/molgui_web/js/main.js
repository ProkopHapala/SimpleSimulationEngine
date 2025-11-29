import { logger } from '../../common_js/Logger.js';
import { MoleculeSystem } from './MoleculeSystem.js';
import { MMParams } from './MMParams.js';
import { MoleculeRenderer } from './MoleculeRenderer.js';
import { IO } from './IO.js';
import { GUI } from './GUI.js';
import { Editor } from './Editor.js';
import { ShortcutManager } from './ShortcutManager.js';

class MolGUIApp {
    constructor() {
        this.container = document.getElementById('canvas-container');
        this.scene = null;
        this.camera = null;
        this.renderer = null;
        this.controls = null;
    }

    async init() {
        logger.info("Initializing MolGUI...");

        // Load Shaders first
        try {
            const vPromise  = fetch('../common_resources/shaders/atom.glslv').then(r => r.text());
            const fPromise  = fetch('../common_resources/shaders/atom.glslf').then(r => r.text());
            const bvPromise = fetch('../common_resources/shaders/bond.glslv').then(r => r.text());
            const svPromise = fetch('../common_resources/shaders/selection.glslv').then(r => r.text());
            const bfPromise = fetch('../common_resources/shaders/bond_color.glslf').then(r => r.text());
            const cfPromise = fetch('../common_resources/shaders/color.glslf').then(r => r.text());
            const lvPromise = fetch('../common_resources/shaders/label.glslv').then(r => r.text());
            const lfPromise = fetch('../common_resources/shaders/label.glslf').then(r => r.text());

            const [vertex, fragment, bVertex, sVertex, bondFrag, colorFrag, lVertex, lFragment] = await Promise.all([
                vPromise, fPromise, bvPromise, svPromise, bfPromise, cfPromise, lvPromise, lfPromise
            ]);

            this.shaders = {
                atom:      { vertex,      fragment },
                bond:      { vertex: bVertex, fragment: bondFrag },   // bond_color.glslf (uses vColor)
                selection: { vertex: sVertex, fragment: colorFrag },  // color.glslf (uses uColor)
                label:     { vertex: lVertex, fragment: lFragment }
            };
            window.logger.info("Shaders loaded.");
        } catch (e) {
            window.logger.error("Failed to load shaders: " + e);
            return;
        }

        // 1. Scene
        this.scene = new THREE.Scene();
        this.scene.background = new THREE.Color(0x222222);

        // 2. Camera
        const width = this.container.clientWidth;
        const height = this.container.clientHeight;
        const aspect = width / height;
        this.frustumSize = 20;
        this.camera = new THREE.OrthographicCamera(
            this.frustumSize * aspect / -2,
            this.frustumSize * aspect / 2,
            this.frustumSize / 2,
            this.frustumSize / -2,
            0.1,
            1000
        );
        this.camera.position.set(10, 10, 10);
        this.camera.lookAt(0, 0, 0);

        // 3. Renderer
        this.renderer = new THREE.WebGLRenderer({
            antialias: true,
            powerPreference: "high-performance"
        });
        this.renderer.setSize(width, height);
        this.renderer.setPixelRatio(window.devicePixelRatio);
        this.container.appendChild(this.renderer.domElement);

        // 4. Controls
        if (THREE.OrbitControls) {
            this.controls = new THREE.OrbitControls(this.camera, this.renderer.domElement);
            this.controls.enableDamping = false; // Disable inertia for raw performance feel
            this.controls.dampingFactor = 0.05;

            // Unity Style Controls:
            // LMB: Selection (Handled by Editor) -> Set to NULL here
            // MMB: Zoom/Dolly (Standard)
            // RMB: Rotate (Standard)
            // Shift + RMB: Pan (Custom)

            this.controls.mouseButtons = {
                LEFT: null, // Disable Left Click for Camera
                MIDDLE: THREE.MOUSE.DOLLY,
                RIGHT: THREE.MOUSE.ROTATE
            };

            // Listen for Shift Key to toggle Pan Mode on RMB
            window.addEventListener('keydown', (e) => {
                if (e.key === 'Shift') {
                    this.controls.mouseButtons.RIGHT = THREE.MOUSE.PAN;
                }
            });

            window.addEventListener('keyup', (e) => {
                if (e.key === 'Shift') {
                    this.controls.mouseButtons.RIGHT = THREE.MOUSE.ROTATE;
                }
            });
        } else {
            window.logger.error("OrbitControls not loaded!");
        }

        // 5. Molecule System
        this.system = new MoleculeSystem();

        // 5. MMParams (Load resources)
        this.mmParams = new MMParams();
        await this.mmParams.loadResources('../common_resources/ElementTypes.dat', '../common_resources/AtomTypes.dat');

        // 6. Molecule Renderer (Pass loaded shaders and mmParams)
        this.molRenderer = new MoleculeRenderer(this.scene, this.system, this.shaders, this.mmParams);

        // 7. IO
        this.io = new IO(this.system, this.molRenderer);
        this.gui = new GUI(this.io);

        // 7. Editor (Selection, Gizmo)
        this.editor = new Editor(this.scene, this.camera, this.renderer, this.system, this.molRenderer);

        // 8. Selection Rendering (Centralized in MoleculeRenderer)
        // No extra code needed here, MoleculeRenderer handles it.

        this.editor.onSelectionChange = () => {
            this.gui.updateSelectionCount();
            this.molRenderer.updateSelection();
        };

        // 9. Shortcut Manager
        this.shortcuts = new ShortcutManager(this.editor);

        // Hook into Editor or System updates?
        // Ideally System should dispatch events, but we can just patch the update method or call it manually.
        // For now, let's monkey-patch the renderer update or just call it in animate loop?
        // Better: Pass it to Editor? Or make Editor call a global update?
        // Let's make Editor call it.
        // this.editor.onSelectionChange is already set above.

        // Also hook GUI input back to renderer
        // Also hook GUI input back to renderer
        this.gui.onSelectionChanged = () => {
            this.molRenderer.updateSelection();
        };

        // Create Test Molecule (Water: H-O-H)
        // O at 0,0,0
        const o = this.system.addAtom(0, 0, 0, 8);
        // H at 0.8, 0.6, 0
        const h1 = this.system.addAtom(0.8, 0.6, 0, 1);
        // H at -0.8, 0.6, 0
        const h2 = this.system.addAtom(-0.8, 0.6, 0, 1);

        this.system.addBond(o, h1);
        this.system.addBond(o, h2);

        // Add some more atoms to test performance/visuals
        // Methane-like structure nearby
        const c = this.system.addAtom(3, 0, 0, 6);
        const h3 = this.system.addAtom(3, 1, 0, 1);
        const h4 = this.system.addAtom(3, -0.5, 0.8, 1);
        const h5 = this.system.addAtom(3, -0.5, -0.8, 1);

        this.system.addBond(c, h3);
        this.system.addBond(c, h4);
        this.system.addBond(c, h5);

        this.molRenderer.update();

        // Lights
        const ambientLight = new THREE.AmbientLight(0xffffff, 0.5);
        this.scene.add(ambientLight);
        const dirLight = new THREE.DirectionalLight(0xffffff, 1);
        dirLight.position.set(5, 10, 7);
        this.scene.add(dirLight);

        // Events
        window.addEventListener('resize', this.onWindowResize.bind(this));

        // Start Loop
        this.animate();

        window.logger.info("Initialization Complete.");
    }

    onWindowResize() {
        const width = this.container.clientWidth;
        const height = this.container.clientHeight;
        const aspect = width / height;

        this.camera.left = -this.frustumSize * aspect / 2;
        this.camera.right = this.frustumSize * aspect / 2;
        this.camera.top = this.frustumSize / 2;
        this.camera.bottom = -this.frustumSize / 2;

        this.camera.updateProjectionMatrix();
        this.renderer.setSize(width, height);
    }

    animate() {
        requestAnimationFrame(this.animate.bind(this));

        if (this.controls) this.controls.update();

        this.renderer.render(this.scene, this.camera);
    }
}

// Start
window.onload = () => {
    const app = new MolGUIApp();
    window.app = app; // Expose for debugging
    app.init();
};
