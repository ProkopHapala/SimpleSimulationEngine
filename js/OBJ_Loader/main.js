
// Global variables
let scene, renderer, controls;
let camera, cameraPersp, cameraOrtho;
let meshRenderer = null;
let loadedObject = null;
let shaders = null;

// Configuration
const MAX_CAPACITY = 65536; // Max vertices/edges

document.addEventListener('DOMContentLoaded', () => {
    init();
    loadShaders();
    setupUI();
    // Don't load default OBJ here - wait for shaders
});

function loadDefaultOBJ() {
    fetch('turret.obj')
        .then(response => response.text())
        .then(text => parseOBJ(text))
        .catch(err => console.error("Failed to load default OBJ:", err));
}

function init() {
    // Scene Setup
    scene = new THREE.Scene();
    scene.background = new THREE.Color(0x222222);

    // Cameras
    const aspect = window.innerWidth / window.innerHeight;
    cameraPersp = new THREE.PerspectiveCamera(60, aspect, 0.1, 1000);
    cameraPersp.position.set(5, 5, 5);

    const frustumSize = 10;
    cameraOrtho = new THREE.OrthographicCamera(
        frustumSize * aspect / -2, frustumSize * aspect / 2,
        frustumSize / 2, frustumSize / -2,
        0.1, 1000
    );
    cameraOrtho.position.set(5, 5, 5);

    camera = cameraPersp; // Default

    // Renderer
    renderer = new THREE.WebGLRenderer({ antialias: true });
    renderer.setSize(window.innerWidth, window.innerHeight);
    document.getElementById('canvas-container').appendChild(renderer.domElement);

    // Controls
    controls = new THREE.OrbitControls(camera, renderer.domElement);
    controls.enableDamping = false; // Default off
    controls.dampingFactor = 0.05;

    // Lights
    const ambientLight = new THREE.AmbientLight(0x404040);
    scene.add(ambientLight);
    const dirLight = new THREE.DirectionalLight(0xffffff, 0.8);
    dirLight.position.set(1, 1, 1);
    scene.add(dirLight);
    const dirLight2 = new THREE.DirectionalLight(0xffffff, 0.5);
    dirLight2.position.set(-1, -1, -1);
    scene.add(dirLight2);

    // Resize handler
    window.addEventListener('resize', onWindowResize, false);

    // Error Capture
    const originalError = console.error;
    console.error = function (...args) {
        originalError.apply(console, args);
        showError(args.join(' '));
    };

    // Animation Loop
    animate();
}

function showError(msg) {
    const log = document.getElementById('error-log');
    const content = document.getElementById('error-content');
    log.style.display = 'block';
    content.textContent += msg + '\n';
    log.scrollTop = log.scrollHeight;
}

function loadShaders() {
    const loadShader = async (url) => {
        const response = await fetch(url);
        return await response.text();
    };

    Promise.all([
        loadShader('../common_resources/shaders/atom.glslv'),
        loadShader('../common_resources/shaders/atom.glslf'),
        loadShader('../common_resources/shaders/bond.glslv'),
        loadShader('../common_resources/shaders/color.glslf'),
        loadShader('../common_resources/shaders/label.glslv'),
        loadShader('../common_resources/shaders/label.glslf')
    ]).then(([atomVert, atomFrag, bondVert, bondFrag, labelVert, labelFrag]) => {
        shaders = {
            atom: { vertex: atomVert, fragment: atomFrag },
            bond: { vertex: bondVert, fragment: bondFrag },
            label: { vertex: labelVert, fragment: labelFrag }
        };
        console.log("Shaders loaded.");

        initMeshRenderer(MAX_CAPACITY);

        // Now it's safe to load default OBJ
        loadDefaultOBJ();
    }).catch(err => {
        console.error("Failed to load shaders:", err);
    });
}

function initMeshRenderer(capacity) {
    if (meshRenderer) {
        if (meshRenderer.atomMesh) scene.remove(meshRenderer.atomMesh);
        if (meshRenderer.bondLines) scene.remove(meshRenderer.bondLines);
        if (meshRenderer.labelMesh) scene.remove(meshRenderer.labelMesh);
    }

    meshRenderer = new MeshRenderer(scene, shaders, capacity);

    // Apply initial styles from UI
    updateLabelStyle();

    updateVisibility();
}

function setupUI() {
    document.getElementById('fileInput').addEventListener('change', handleFileSelect, false);

    // Visibility
    document.getElementById('chkNodes').addEventListener('change', updateVisibility);
    document.getElementById('chkEdges').addEventListener('change', updateVisibility);
    document.getElementById('chkFaces').addEventListener('change', updateVisibility);
    document.getElementById('chkLabels').addEventListener('change', updateVisibility);

    // Camera
    document.getElementById('selCameraMode').addEventListener('change', (e) => {
        setCameraMode(e.target.value);
    });

    document.getElementById('btnViewXY').addEventListener('click', () => setView([0, 0, 20], [0, 1, 0]));
    document.getElementById('btnViewXZ').addEventListener('click', () => setView([0, 20, 0], [0, 0, -1]));
    document.getElementById('btnViewYZ').addEventListener('click', () => setView([20, 0, 0], [0, 1, 0]));

    document.getElementById('zoomInput').addEventListener('input', (e) => {
        const val = parseFloat(e.target.value);
        if (!isNaN(val)) {
            camera.zoom = val;
            camera.updateProjectionMatrix();
        }
    });

    document.getElementById('chkDamping').addEventListener('change', (e) => {
        controls.enableDamping = e.target.checked;
    });

    // Labels
    document.getElementById('colLabel').addEventListener('input', updateLabelStyle);
    document.getElementById('numLabelSize').addEventListener('input', updateLabelStyle);
    document.getElementById('chkLabelFixed').addEventListener('change', updateLabelStyle);
}

function updateLabelStyle() {
    if (!meshRenderer) return;

    const color = document.getElementById('colLabel').value;
    const size = parseFloat(document.getElementById('numLabelSize').value);
    const fixed = document.getElementById('chkLabelFixed').checked;

    meshRenderer.setLabelStyle(color, size, fixed);
}

function setCameraMode(mode) {
    const oldCam = camera;

    if (mode === 'ortho') {
        camera = cameraOrtho;
    } else {
        camera = cameraPersp;
    }

    // Copy basic transform
    camera.position.copy(oldCam.position);
    camera.rotation.copy(oldCam.rotation);
    camera.zoom = oldCam.zoom;

    controls.object = camera;
    camera.updateProjectionMatrix();
}

function setView(pos, up) {
    camera.position.set(pos[0], pos[1], pos[2]);
    camera.up.set(up[0], up[1], up[2]);
    camera.lookAt(0, 0, 0);
    camera.updateProjectionMatrix();
    controls.update();
}

function updateVisibility() {
    if (!meshRenderer) return;

    const showNodes = document.getElementById('chkNodes').checked;
    const showEdges = document.getElementById('chkEdges').checked;
    const showFaces = document.getElementById('chkFaces').checked;
    const showLabels = document.getElementById('chkLabels').checked;

    if (meshRenderer.atomMesh) meshRenderer.atomMesh.visible = showNodes;
    if (meshRenderer.bondLines) meshRenderer.bondLines.visible = showEdges;
    if (meshRenderer.labelMesh) {
        meshRenderer.labelMesh.visible = showLabels;
    }

    if (loadedObject) {
        loadedObject.visible = showFaces;
    }
}

function handleFileSelect(evt) {
    const file = evt.target.files[0];
    if (!file) return;

    const reader = new FileReader();
    reader.onload = function (e) {
        const contents = e.target.result;
        parseOBJ(contents);
    };
    reader.readAsText(file);
}

function parseOBJ(text) {
    if (!meshRenderer) {
        console.error("MeshRenderer not initialized yet!");
        return;
    }

    if (loadedObject) {
        scene.remove(loadedObject);
        loadedObject = null;
    }

    const loader = new THREE.OBJLoader();
    const object = loader.parse(text);

    loadedObject = object;
    scene.add(loadedObject);

    processGeometry(object);

    // Center camera
    const box = new THREE.Box3().setFromObject(object);
    const center = box.getCenter(new THREE.Vector3());
    const size = box.getSize(new THREE.Vector3());
    const maxDim = Math.max(size.x, size.y, size.z);

    controls.target.copy(center);
    camera.position.copy(center);
    camera.position.z += maxDim * 2;
    camera.updateProjectionMatrix();

    updateVisibility();
}

function processGeometry(object) {
    if (!meshRenderer) {
        console.error("MeshRenderer not initialized!");
        return;
    }

    let vertices = [];
    let edges = [];

    const uniqueVerts = [];
    const vertMap = new Map();

    function getVertIndex(x, y, z) {
        const key = `${x.toFixed(6)},${y.toFixed(6)},${z.toFixed(6)}`;
        if (vertMap.has(key)) return vertMap.get(key);

        const idx = uniqueVerts.length / 3;
        uniqueVerts.push(x, y, z);
        vertMap.set(key, idx);
        return idx;
    }

    object.traverse((child) => {
        if (child.isMesh) {
            child.material = new THREE.MeshPhongMaterial({
                color: 0x888888,
                side: THREE.DoubleSide,
                polygonOffset: true,
                polygonOffsetFactor: 1,
                polygonOffsetUnits: 1
            });

            const pos = child.geometry.attributes.position;
            const count = pos.count;

            for (let i = 0; i < count; i += 3) {
                const ax = pos.getX(i); const ay = pos.getY(i); const az = pos.getZ(i);
                const bx = pos.getX(i + 1); const by = pos.getY(i + 1); const bz = pos.getZ(i + 1);
                const cx = pos.getX(i + 2); const cy = pos.getY(i + 2); const cz = pos.getZ(i + 2);

                const ia = getVertIndex(ax, ay, az);
                const ib = getVertIndex(bx, by, bz);
                const ic = getVertIndex(cx, cy, cz);

                edges.push([ia, ib]);
                edges.push([ib, ic]);
                edges.push([ic, ia]);
            }
        }
    });

    const nodeCount = uniqueVerts.length / 3;

    if (nodeCount > meshRenderer.capacity) {
        console.warn("Capacity exceeded, resizing...");
        initMeshRenderer(nodeCount * 2);
    }

    // 1. Update Positions
    meshRenderer.updatePositions(uniqueVerts, nodeCount);

    // 2. Update Nodes
    meshRenderer.updateParticles(nodeCount,
        (i) => [1, 0, 0],
        (i) => 0.05
    );

    // 3. Update Edges
    const uniqueEdges = new Map();
    const finalEdges = [];
    for (const e of edges) {
        const k = e[0] < e[1] ? `${e[0]}_${e[1]}` : `${e[1]}_${e[0]}`;
        if (!uniqueEdges.has(k)) {
            uniqueEdges.set(k, true);
            finalEdges.push(e);
        }
    }

    meshRenderer.updateBonds(finalEdges);

    // 4. Update Labels
    meshRenderer.updateLabels((i) => i.toString(), nodeCount);

    // Ensure styles are applied
    updateLabelStyle();

    // Make sure visibility is synced
    updateVisibility();
}

function onWindowResize() {
    const aspect = window.innerWidth / window.innerHeight;

    cameraPersp.aspect = aspect;
    cameraPersp.updateProjectionMatrix();

    const frustumSize = 10;
    cameraOrtho.left = -frustumSize * aspect / 2;
    cameraOrtho.right = frustumSize * aspect / 2;
    cameraOrtho.top = frustumSize / 2;
    cameraOrtho.bottom = -frustumSize / 2;
    cameraOrtho.updateProjectionMatrix();

    renderer.setSize(window.innerWidth, window.innerHeight);

    if (meshRenderer) meshRenderer.updateLabelUniforms(aspect);
}

function animate() {
    requestAnimationFrame(animate);
    controls.update();

    // Update label uniforms every frame for proper billboard effect
    if (meshRenderer && meshRenderer.labelMesh) {
        meshRenderer.updateLabelUniforms(window.innerWidth / window.innerHeight);
    }

    renderer.render(scene, camera);
}
