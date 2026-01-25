import { SHADER_HEADER } from './gpu-core.js';

export const COMMON_VS = `
struct VertexOutput {
    @builtin(position) pos: vec4f,
    @location(0) uv: vec2f,
    @location(1) worldPos: vec2f
};

fn worldToScreen(p: vec2f) -> vec4f {
    let camPos = vec2f(globals.camX, globals.camY);
    let viewSize = vec2f(globals.resX, globals.resY);
    let aspect = viewSize.x / viewSize.y;
    
    let rel = p - camPos;
    let visibleHeight = 1000.0 / globals.zoom; 
    let scaleY = 2.0 / visibleHeight; 
    let scaleX = scaleY / aspect;

    let ndcX = rel.x * scaleX;
    let ndcY = rel.y * scaleY;
    
    return vec4f(ndcX, ndcY, 0.0, 1.0);
}
`;

export class ViewportRenderer {
    constructor(gpu) {
        this.gpu = gpu;
        const code = `
        ${SHADER_HEADER}
        struct VSOut { @builtin(position) pos: vec4f, @location(0) uv: vec2f };

        @vertex fn vs(@builtin(vertex_index) idx: u32) -> VSOut {
            var pos = array<vec2f, 4>(vec2f(-1,-1), vec2f(1,-1), vec2f(-1,1), vec2f(1,1));
            var out: VSOut;
            out.pos = vec4f(pos[idx], 0.0, 1.0);
            
            // Reconstruct World Pos for texture sampling
            let viewSize = vec2f(globals.resX, globals.resY);
            let aspect = viewSize.x / viewSize.y;
            let visibleHeight = 1000.0 / globals.zoom;
            let visibleWidth = visibleHeight * aspect;
            
            // Screen UV 0..1 (Flip Y for calculation)
            var sUV = pos[idx] * 0.5 + 0.5;
            sUV.y = 1.0 - sUV.y; 
            
            let worldPos = vec2f(globals.camX, globals.camY) + (sUV - 0.5) * vec2f(visibleWidth, visibleHeight);
            
            // Map World 0..mapSize to UV 0..1
            out.uv = worldPos / vec2f(globals.resX, globals.resY);
            
            return out;
        }

        @group(1) @binding(0) var tex: texture_2d<f32>;
        @group(1) @binding(1) var samp: sampler;

        @fragment fn fs(in: VSOut) -> @location(0) vec4f {
            // Check bounds
            if(in.uv.x < 0.0 || in.uv.x > 1.0 || in.uv.y < 0.0 || in.uv.y > 1.0) {
                return vec4f(0.01, 0.01, 0.02, 1.0); // Background color
            }
            
            // FIX: Use textureSampleLevel to avoid non-uniform control flow errors
            // Level 0.0 is the base texture, which is exactly what we want for a map view
            return textureSampleLevel(tex, samp, in.uv, 0.0);
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
            vertex: { module: this.module, entryPoint: 'vs' },
            fragment: { module: this.module, entryPoint: 'fs', targets: [{ format: gpu.format }] },
            primitive: { topology: 'triangle-strip' }
        });
        
        this.sampler = gpu.device.createSampler({ magFilter: 'linear', minFilter: 'linear' });
    }

    draw(pass, visualField) {
        const bg = this.gpu.device.createBindGroup({
            layout: this.layout,
            entries: [
                { binding: 0, resource: visualField.view },
                { binding: 1, resource: this.sampler }
            ]
        });
        pass.setPipeline(this.pipeline);
        pass.setBindGroup(0, this.gpu.globalBindGroup);
        pass.setBindGroup(1, bg);
        pass.draw(4);
    }
}

export class LineRenderer {
    constructor(gpu, maxLines = 10000) {
        this.gpu = gpu;
        this.count = 0;
        this.buffer = gpu.device.createBuffer({ size: maxLines * 16, usage: GPUBufferUsage.VERTEX | GPUBufferUsage.COPY_DST });

        const code = `
        ${SHADER_HEADER}
        ${COMMON_VS}
        @vertex fn vs(@location(0) posIn: vec2f) -> VertexOutput {
            var out: VertexOutput;
            out.pos = worldToScreen(posIn);
            return out;
        }
        @fragment fn fs() -> @location(0) vec4f { return vec4f(0.0, 0.5, 1.0, 1.0); }
        `;

        this.pipeline = gpu.device.createRenderPipeline({
            layout: gpu.device.createPipelineLayout({ bindGroupLayouts: [gpu.globalLayout] }),
            vertex: {
                module: gpu.device.createShaderModule({ code }), entryPoint: 'vs',
                buffers: [{ arrayStride: 8, attributes: [{ shaderLocation: 0, offset: 0, format: 'float32x2' }] }]
            },
            fragment: { module: gpu.device.createShaderModule({ code }), entryPoint: 'fs', targets: [{ format: gpu.format }] },
            primitive: { topology: 'line-list' }
        });
    }

    updateData(points) {
        this.gpu.device.queue.writeBuffer(this.buffer, 0, points);
        this.count = points.length / 2;
    }

    draw(pass) {
        if(this.count === 0) return;
        pass.setPipeline(this.pipeline);
        pass.setBindGroup(0, this.gpu.globalBindGroup);
        pass.setVertexBuffer(0, this.buffer);
        pass.draw(this.count);
    }
}