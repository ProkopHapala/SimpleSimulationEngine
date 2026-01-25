export const SHADER_HEADER = `
struct GlobalParams {
    time: f32,
    dt: f32,
    frame: u32,
    resX: f32,
    resY: f32,
    mouseX: f32,
    mouseY: f32,
    camX: f32,
    camY: f32,
    zoom: f32,
    minHeight: f32,
    maxHeight: f32
};
@group(0) @binding(0) var<uniform> globals : GlobalParams;
`;

export class GPUContext {
    constructor(canvasId) {
        this.canvas = document.getElementById(canvasId);
        if (!this.canvas) throw new Error(`Canvas ${canvasId} not found`);
        
        this.device = null;
        this.context = null;
        this.format = navigator.gpu.getPreferredCanvasFormat();
        
        this.mapSize = { width: 1024, height: 1024 };
        this.camera = { x: 512, y: 512, zoom: 1.0 }; // Start centered
        this.mouse = { screenX: 0, screenY: 0, worldX: 0, worldY: 0, down: false };
        
        this.uniformValues = new ArrayBuffer(64);
        this.views = {
            time: new Float32Array(this.uniformValues, 0, 1),
            dt: new Float32Array(this.uniformValues, 4, 1),
            frame: new Uint32Array(this.uniformValues, 8, 1),
            res: new Float32Array(this.uniformValues, 12, 2),
            mouse: new Float32Array(this.uniformValues, 20, 2),
            cam: new Float32Array(this.uniformValues, 28, 2),
            zoom: new Float32Array(this.uniformValues, 36, 1),
            heightRange: new Float32Array(this.uniformValues, 40, 2)
        };
        
        this.frameCount = 0;
        this.startTime = performance.now();
        this.lastTime = this.startTime;
    }

    async init() {
        if (!navigator.gpu) throw new Error("WebGPU not supported");
        const adapter = await navigator.gpu.requestAdapter({ powerPreference: "high-performance" });
        this.device = await adapter.requestDevice();
        this.context = this.canvas.getContext('webgpu');
        this._resizeCanvas(true);
        
        this.uniformBuffer = this.device.createBuffer({
            size: 64, usage: GPUBufferUsage.UNIFORM | GPUBufferUsage.COPY_DST
        });

        this.globalLayout = this.device.createBindGroupLayout({
            entries: [{ binding: 0, visibility: GPUShaderStage.COMPUTE | GPUShaderStage.VERTEX | GPUShaderStage.FRAGMENT, buffer: { type: 'uniform' } }]
        });

        this.globalBindGroup = this.device.createBindGroup({
            layout: this.globalLayout,
            entries: [{ binding: 0, resource: { buffer: this.uniformBuffer } }]
        });

        window.addEventListener('resize', () => this._resizeCanvas());
        this._setupInput();
    }

    update() {
        const now = performance.now();
        const dt = (now - this.lastTime) / 1000;
        this.lastTime = now;
        this.frameCount++;

        this._resizeCanvas();

        // Mouse to World Calculation
        const rect = this.canvas.getBoundingClientRect();
        const ndcX = (this.mouse.screenX / rect.width) * 2 - 1;
        const ndcY = -((this.mouse.screenY / rect.height) * 2 - 1); // WebGPU Clip Y is up
        
        // Aspect Ratio Correction
        const aspect = rect.width / rect.height;
        
        // Calculate visible world height
        // We define zoom=1.0 as "1000 world units height fits screen"
        const visibleHeight = 1000.0 / this.camera.zoom;
        const visibleWidth = visibleHeight * aspect;

        this.mouse.worldX = this.camera.x + (ndcX * visibleWidth * 0.5);
        this.mouse.worldY = this.camera.y + (ndcY * visibleHeight * 0.5);

        this.views.time[0] = (now - this.startTime) / 1000;
        this.views.dt[0] = dt;
        this.views.frame[0] = this.frameCount;
        this.views.res[0] = this.mapSize.width;
        this.views.res[1] = this.mapSize.height;
        this.views.mouse[0] = this.mouse.worldX;
        this.views.mouse[1] = this.mouse.worldY;
        this.views.cam[0] = this.camera.x;
        this.views.cam[1] = this.camera.y;
        this.views.zoom[0] = this.camera.zoom;

        this.device.queue.writeBuffer(this.uniformBuffer, 0, this.uniformValues);
    }

    getRenderPass(commandEncoder) {
        return commandEncoder.beginRenderPass({
            colorAttachments: [{
                view: this.context.getCurrentTexture().createView(),
                clearValue: [0.02, 0.02, 0.05, 1], // Deep dark blue background
                loadOp: 'clear', storeOp: 'store'
            }]
        });
    }

    _setupInput() {
        this.canvas.addEventListener('wheel', e => {
            e.preventDefault();
            this.camera.zoom *= (1 - e.deltaY * 0.001);
            this.camera.zoom = Math.max(0.05, Math.min(this.camera.zoom, 20.0));
        }, { passive: false });
        
        this.canvas.addEventListener('mousedown', () => this.mouse.down = true);
        window.addEventListener('mouseup', () => this.mouse.down = false);
        this.canvas.addEventListener('mousemove', e => {
            const r = this.canvas.getBoundingClientRect();
            this.mouse.screenX = e.clientX - r.left;
            this.mouse.screenY = e.clientY - r.top;
            if (this.mouse.down) {
                // Drag speed depends on zoom
                const visibleHeight = 1000.0 / this.camera.zoom;
                const scale = visibleHeight / r.height;
                this.camera.x -= e.movementX * scale;
                this.camera.y += e.movementY * scale; // Inverted to match screen drag feel
            }
        });
    }

    _resizeCanvas(force = false) {
        if (!this.context) return;
        const dpr = window.devicePixelRatio || 1;
        const displayWidth = Math.max(1, Math.floor(this.canvas.clientWidth * dpr));
        const displayHeight = Math.max(1, Math.floor(this.canvas.clientHeight * dpr));
        if (!force && this.canvas.width === displayWidth && this.canvas.height === displayHeight) return;
        this.canvas.width = displayWidth;
        this.canvas.height = displayHeight;
        this.context.configure({
            device: this.device,
            format: this.format,
            size: [displayWidth, displayHeight],
            alphaMode: 'opaque'
        });
    }
}

// Data Container for Simulation (Ping-Pong)
export class Field {
    constructor(gpu, label = 'Field', format = 'r32float') {
        this.gpu = gpu;
        this.width = gpu.mapSize.width;
        this.height = gpu.mapSize.height;
        this.format = format;
        
        const desc = {
            label: label,
            size: [this.width, this.height],
            format: format,
            usage: GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.STORAGE_BINDING | GPUTextureUsage.COPY_SRC
        };
        
        this.texA = gpu.device.createTexture(desc);
        this.texB = gpu.device.createTexture(desc);
        this.readTex = this.texA;
        this.writeTex = this.texB;
        this.readBuffer = null;
        this.readBufferSize = 0;
        this.pendingReadback = false;
    }

    swap() {
        const t = this.readTex;
        this.readTex = this.writeTex;
        this.writeTex = t;
    }

    getReadView() { return this.readTex.createView(); }
    getWriteView() { return this.writeTex.createView(); }

    _ensureReadBuffer() {
        const size = this.width * this.height * 4;
        if (!this.readBuffer || this.readBufferSize !== size) {
            this.readBuffer?.destroy?.();
            this.readBuffer = this.gpu.device.createBuffer({
                size,
                usage: GPUBufferUsage.COPY_DST | GPUBufferUsage.MAP_READ
            });
            this.readBufferSize = size;
        }
        return this.readBuffer;
    }

    encodeReadback(commandEncoder) {
        if (this.pendingReadback) { return false; }
        const buffer = this._ensureReadBuffer();
        commandEncoder.copyTextureToBuffer(
            { texture: this.readTex },
            { buffer, bytesPerRow: this.width * 4, rowsPerImage: this.height },
            { width: this.width, height: this.height, depthOrArrayLayers: 1 }
        );
        this.pendingReadback = true;
        return true;
    }

    async collectReadback() {
        if (!this.pendingReadback || !this.readBuffer) { return null; }
        await this.readBuffer.mapAsync(GPUMapMode.READ);
        const mapped = this.readBuffer.getMappedRange();
        const copy = mapped.slice(0);
        this.readBuffer.unmap();
        this.pendingReadback = false;
        return new Float32Array(copy);
    }
}

// Visual Container (Single texture, filterable)
export class VisualField {
    constructor(gpu, label = 'Visual', format = 'rgba8unorm') {
        this.gpu = gpu;
        this.tex = gpu.device.createTexture({
            label: label,
            size: [gpu.mapSize.width, gpu.mapSize.height],
            format: format,
            usage: GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.STORAGE_BINDING
        });
        this.view = this.tex.createView();
        this.format = format;
    }
}