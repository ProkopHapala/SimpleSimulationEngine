import { SHADER_HEADER } from './gpu-core.js';

export class MapAlgorithm {
    constructor(gpu, wgslBody, inputs = [], outputs = []) {
        this.gpu = gpu;
        this.inputs = inputs;
        this.outputs = outputs;

        let header = `\n// --- Auto-Generated Bindings ---\n`;
        let bIdx = 0;
        inputs.forEach((f, i) => header += `@group(1) @binding(${bIdx++}) var in_${i} : texture_2d<f32>;\n`);
        outputs.forEach((f, i) => header += `@group(1) @binding(${bIdx++}) var out_${i} : texture_storage_2d<${f.format}, write>;\n`);

        const fullCode = `${SHADER_HEADER}\n${header}\n${wgslBody}`;
        
        const entries = [];
        bIdx = 0;
        inputs.forEach(() => entries.push({ 
            binding: bIdx++, visibility: GPUShaderStage.COMPUTE, 
            texture: { sampleType: 'unfilterable-float', viewDimension: '2d' } 
        }));
        outputs.forEach(f => entries.push({ 
            binding: bIdx++, visibility: GPUShaderStage.COMPUTE, 
            storageTexture: { access: 'write-only', format: f.format, viewDimension: '2d' } 
        }));

        this.group1Layout = gpu.device.createBindGroupLayout({ entries });
        this.pipeline = gpu.device.createComputePipeline({
            layout: gpu.device.createPipelineLayout({ bindGroupLayouts: [gpu.globalLayout, this.group1Layout] }),
            compute: { module: gpu.device.createShaderModule({ code: fullCode }), entryPoint: 'main' }
        });
    }

    run(commandEncoder) {
        const entries = [];
        let bIdx = 0;
        this.inputs.forEach(f => entries.push({ binding: bIdx++, resource: f.getReadView ? f.getReadView() : f.view })); // Support Field or VisualField
        this.outputs.forEach(f => entries.push({ binding: bIdx++, resource: f.getWriteView ? f.getWriteView() : f.view }));

        const bg1 = this.gpu.device.createBindGroup({ layout: this.group1Layout, entries });

        const pass = commandEncoder.beginComputePass();
        pass.setPipeline(this.pipeline);
        pass.setBindGroup(0, this.gpu.globalBindGroup);
        pass.setBindGroup(1, bg1);
        pass.dispatchWorkgroups(Math.ceil(this.gpu.mapSize.width / 8), Math.ceil(this.gpu.mapSize.height / 8));
        pass.end();

        this.outputs.forEach(f => { if(f.swap) f.swap(); });
    }
}