import { MeshRenderer } from '../common_js/MeshRenderer.js';

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
    cameraPersp = new THREE.PerspectiveCamera(60, aspect, 0.1, 1000000.0);
    cameraPersp.position.set(5, 5, 5);

    const frustumSize = 10;
    cameraOrtho = new THREE.OrthographicCamera(
        frustumSize * aspect / -2, frustumSize * aspect / 2,
        frustumSize / 2, frustumSize / -2,
        0.1, 1000000.0
    );
    cameraOrtho.position.set(5, 5, 5);

    camera = cameraPersp; // Default

    // Renderer is always created to match full browser viewport size.
    // We expose antialiasing as a runtime toggle, see recreateRenderer().
    recreateRenderer(true);

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

// (Re)create the WebGL renderer with the requested antialias setting.
// The canvas drawing buffer is sized to match the visible viewport (canvas-container)
// with a 1:1 mapping between framebuffer pixels and CSS pixels (no extra resampling).
function recreateRenderer(useAntialias) {
    const container = document.getElementById('canvas-container');
    if (!container) return;

    if (renderer && renderer.domElement && renderer.domElement.parentNode === container) {
        container.removeChild(renderer.domElement);
    }

    const width = container.clientWidth || window.innerWidth;
    const height = container.clientHeight || window.innerHeight;

    renderer = new THREE.WebGLRenderer({ antialias: !!useAntialias });
    // Use physical display resolution for crisp, thin lines on high-DPI screens.
    const dpr = Math.min(window.devicePixelRatio || 1, 2.0);
    renderer.setPixelRatio(dpr);
    renderer.setSize(width, height);
    container.appendChild(renderer.domElement);

    // Re-bind controls if they already exist
    if (controls) {
        controls.dispose();
        controls = new THREE.OrbitControls(camera, renderer.domElement);
        controls.enableDamping = document.getElementById('chkDamping')?.checked || false;
        controls.dampingFactor = 0.05;
    }
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
        loadShader('../common_resources/shaders/bond_color.glslf'),
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
        if (controls) controls.enableDamping = e.target.checked;
    });

    // Antialiasing toggle: recreates renderer with or without MSAA.
    const chkAA = document.getElementById('chkAA');
    if (chkAA) {
        chkAA.addEventListener('change', (e) => {
            recreateRenderer(e.target.checked);
        });
    }

    // Centering actions
    document.getElementById('btnCenterCamera').addEventListener('click', () => {
        centerCameraOnObject();
    });

    document.getElementById('btnCenterCOG').addEventListener('click', () => {
        centerObjectCOGToOrigin();
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

    // Debug information about display and scaling for this OBJ load.
    console.log('DPR:', window.devicePixelRatio);
    console.log('CSS viewport:', window.innerWidth + 'x' + window.innerHeight);
    console.log('Physical screen approx:', screen.width + 'x' + screen.height);
    console.log('Scaling factor:', screen.width / window.innerWidth);

    if (loadedObject) {
        scene.remove(loadedObject);
        loadedObject = null;
    }

    const loader = new THREE.OBJLoader();
    const object = loader.parse(text);

    loadedObject = object;
    scene.add(loadedObject);

    processGeometry(object);
    updateVisibility();
}

// --- Camera and object centering helpers (invoked from UI buttons) ---

function computeObjectBounds() {
    if (!loadedObject) {
        console.warn("No OBJ loaded yet.");
        return null;
    }

    const box = new THREE.Box3().setFromObject(loadedObject);
    if (box.isEmpty()) {
        console.warn("Loaded OBJ has empty bounds.");
        return null;
    }
    return box;
}

function fitCameraToBounds(box) {
    if (!box) return;

    const center = box.getCenter(new THREE.Vector3());
    const size = box.getSize(new THREE.Vector3());
    const maxDim = Math.max(size.x, size.y, size.z) || 1.0;

    // Always orbit around the object center
    controls.target.copy(center);

    if (camera.isPerspectiveCamera) {
        const direction = new THREE.Vector3()
            .subVectors(camera.position, controls.target)
            .normalize();

        if (!isFinite(direction.lengthSq()) || direction.lengthSq() === 0) {
            direction.set(0, 0, 1);
        }

        const fov = camera.fov * Math.PI / 180.0;
        let distance = maxDim / (2 * Math.tan(fov / 2));
        distance *= 2.0; // add a safety factor so the whole object is visible

        camera.position.copy(center).addScaledVector(direction, distance);

        // Set clipping planes relative to object size to avoid far/near clipping.
        camera.near = Math.max(maxDim * 1e-4, 0.01);
        camera.far = Math.max(distance + maxDim * 4.0, maxDim * 10.0);

    } else if (camera.isOrthographicCamera) {
        const aspect = window.innerWidth / window.innerHeight;
        const frustumSize = 10;

        camera.left = -frustumSize * aspect / 2;
        camera.right = frustumSize * aspect / 2;
        camera.top = frustumSize / 2;
        camera.bottom = -frustumSize / 2;

        const zoomX = frustumSize * aspect / (size.x || 1.0);
        const zoomY = frustumSize / (size.y || 1.0);
        const zoom = Math.min(zoomX, zoomY) * 0.9; // margin

        if (zoom > 0 && isFinite(zoom)) {
            camera.zoom = zoom;
            const zoomInput = document.getElementById('zoomInput');
            if (zoomInput) zoomInput.value = camera.zoom.toFixed(2);
        }

        // For ortho camera, keep generous clipping planes based on object span.
        camera.near = Math.max(maxDim * 1e-3, 0.01);
        camera.far = Math.max(maxDim * 10.0, 100.0);
    }

    camera.updateProjectionMatrix();
    controls.update();
}

// Center only the camera on the current object without modifying geometry
function centerCameraOnObject() {
    const box = computeObjectBounds();
    if (!box) return;
    fitCameraToBounds(box);
}

// Move the object so that its COG (approximated by bounding-box center)
// is at the world origin and then fit the camera.
function centerObjectCOGToOrigin() {
    const box = computeObjectBounds();
    if (!box || !loadedObject) return;

    const center = box.getCenter(new THREE.Vector3());

    // Shift the whole loaded object so that COG ~ (0,0,0)
    loadedObject.position.sub(center);

    // After moving, recompute bounds (size is unchanged, but center moved)
    const newBox = new THREE.Box3().setFromObject(loadedObject);

    // Target origin and fit camera to new bounds
    controls.target.set(0, 0, 0);
    fitCameraToBounds(newBox);
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

    const container = document.getElementById('canvas-container');
    const width = (container && container.clientWidth) ? container.clientWidth : window.innerWidth;
    const height = (container && container.clientHeight) ? container.clientHeight : window.innerHeight;

    if (renderer) {
        const dpr = Math.min(window.devicePixelRatio || 1, 2.0);
        renderer.setPixelRatio(dpr);
        renderer.setSize(width, height);
    }

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
