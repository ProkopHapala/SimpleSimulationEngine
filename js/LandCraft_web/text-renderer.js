import { SHADER_HEADER } from './gpu-core.js';
import { COMMON_VS } from './renderers.js';

export class TextRenderer {
    constructor(gpu, maxChars = 2048) {
        this.gpu = gpu;
        this.maxChars = maxChars;
        this.charCount = 0;

        // 1. Generate Font Atlas (Canvas CPU -> Texture GPU)
        this.atlasSize = 512;
        this.gridSize = 32; // 32x32 pixels per char
        this.cols = 16;
        this.rows = 16;
        this.charMap = {}; // 'A' -> index
        this.atlasTexture = this._createFontAtlas();

        // 2. Instance Buffer
        // Struct: [posX, posY, uvIdx, scale]
        this.buffer = gpu.device.createBuffer({
            size: maxChars * 4 * 4,
            usage: GPUBufferUsage.VERTEX | GPUBufferUsage.COPY_DST
        });

        // 3. Pipeline
        const code = `
        ${SHADER_HEADER}
        ${COMMON_VS}

        struct InstanceInput {
            @location(0) pos: vec2f,
            @location(1) uvIdx: f32,
            @location(2) scale: f32
        };

        struct VSOut {
            @builtin(position) pos: vec4f,
            @location(0) uv: vec2f
        };

        @vertex fn vs(@builtin(vertex_index) vIdx: u32, inst: InstanceInput) -> VSOut {
            // Quad vertices
            var quad = array<vec2f, 6>(
                vec2f(0,0), vec2f(1,0), vec2f(0,1),
                vec2f(0,1), vec2f(1,0), vec2f(1,1)
            );
            
            let charSizeWorld = 20.0 * inst.scale; // Base size
            let localPos = quad[vIdx] * charSizeWorld;
            
            // Text is top-left anchored usually, let's center vertically
            let worldPos = inst.pos + localPos - vec2f(0.0, charSizeWorld * 0.5);
            
            var out: VSOut;
            out.pos = worldToScreen(worldPos);

            // Calculate UVs in Atlas (flip vertically because canvas Y is down)
            let idx = u32(inst.uvIdx);
            let cols = 16u;
            let col = idx % cols;
            let row = idx / cols;
            let uvSize = 1.0 / f32(cols);
            let baseUV = vec2f(f32(col), f32(row + 1u)) * uvSize;
            out.uv = baseUV + vec2f(quad[vIdx].x, -quad[vIdx].y) * uvSize;
            
            return out;
        }

        @group(1) @binding(0) var tex: texture_2d<f32>;
        @group(1) @binding(1) var samp: sampler;

        @fragment fn fs(in: VSOut) -> @location(0) vec4f {
            let alpha = textureSample(tex, samp, in.uv).a;
            if(alpha < 0.1) { discard; }
            return vec4f(1.0, 1.0, 1.0, alpha);
        }
        `;

        this.module = gpu.device.createShaderModule({ code });
        this.layout = gpu.device.createBindGroupLayout({
            entries: [
                { binding: 0, visibility: GPUShaderStage.FRAGMENT, texture: { sampleType: 'float', viewDimension: '2d' } },
                { binding: 1, visibility: GPUShaderStage.FRAGMENT, sampler: {} }
            ]
        });

        this.pipeline = gpu.device.createRenderPipeline({
            layout: gpu.device.createPipelineLayout({ bindGroupLayouts: [gpu.globalLayout, this.layout] }),
            vertex: {
                module: this.module, entryPoint: 'vs',
                buffers: [{ 
                    arrayStride: 16, stepMode: 'instance',
                    attributes: [
                        { shaderLocation: 0, offset: 0, format: 'float32x2' }, // pos
                        { shaderLocation: 1, offset: 8, format: 'float32' },   // uvIdx
                        { shaderLocation: 2, offset: 12, format: 'float32' }   // scale
                    ] 
                }]
            },
            fragment: { module: this.module, entryPoint: 'fs', targets: [{ 
                format: gpu.format,
                blend: { color: { srcFactor: 'src-alpha', dstFactor: 'one-minus-src-alpha', operation: 'add' }, alpha: { srcFactor: 'one', dstFactor: 'one', operation: 'add' } }
            }] },
            primitive: { topology: 'triangle-list' }
        });

        this.sampler = gpu.device.createSampler({ magFilter: 'linear', minFilter: 'linear' });
        this.bindGroup = gpu.device.createBindGroup({
            layout: this.layout,
            entries: [{ binding: 0, resource: this.atlasTexture.createView() }, { binding: 1, resource: this.sampler }]
        });
        
        // Staging for CPU text construction
        this.cpuBuffer = new Float32Array(maxChars * 4);
    }

    _createFontAtlas() {
        const c = document.createElement('canvas');
        c.width = this.atlasSize;
        c.height = this.atlasSize;
        const ctx = c.getContext('2d');
        ctx.font = 'bold 24px monospace';
        ctx.textAlign = 'center';
        ctx.textBaseline = 'middle';
        ctx.fillStyle = 'white';
        
        // ASCII printable range
        let idx = 0;
        for(let i=32; i<127; i++) {
            const char = String.fromCharCode(i);
            this.charMap[char] = idx;
            const x = (idx % this.cols) * this.gridSize + this.gridSize/2;
            const y = Math.floor(idx / this.cols) * this.gridSize + this.gridSize/2;
            ctx.fillText(char, x, y);
            idx++;
        }
        
        // Create GPU Texture
        const tex = this.gpu.device.createTexture({
            size: [this.atlasSize, this.atlasSize],
            format: 'rgba8unorm',
            usage: GPUTextureUsage.TEXTURE_BINDING | GPUTextureUsage.COPY_DST | GPUTextureUsage.RENDER_ATTACHMENT
        });
        this.gpu.device.queue.copyExternalImageToTexture({ source: c }, { texture: tex }, [this.atlasSize, this.atlasSize]);
        return tex;
    }

    // Call this every frame to rebuild text buffer
    begin() { this.charCount = 0; }

    addText(str, x, y, scale = 1.0) {
        let cursorX = x;
        for (let i = 0; i < str.length; i++) {
            if(this.charCount >= this.maxChars) break;
            const char = str[i];
            const idx = this.charMap[char] || 0;
            
            this.cpuBuffer[this.charCount*4 + 0] = cursorX;
            this.cpuBuffer[this.charCount*4 + 1] = y;
            this.cpuBuffer[this.charCount*4 + 2] = idx;
            this.cpuBuffer[this.charCount*4 + 3] = scale;
            
            cursorX += 12 * scale; // Simple monospace advance
            this.charCount++;
        }
    }

    draw(pass) {
        if(this.charCount === 0) return;
        this.gpu.device.queue.writeBuffer(this.buffer, 0, this.cpuBuffer, 0, this.charCount * 4);
        
        pass.setPipeline(this.pipeline);
        pass.setBindGroup(0, this.gpu.globalBindGroup);
        pass.setBindGroup(1, this.bindGroup);
        pass.setVertexBuffer(0, this.buffer);
        pass.draw(6, this.charCount);
    }
}