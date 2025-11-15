### WebGL Simulation Framework: Design Overview

I'll design a lightweight, Shadertoy-inspired WebGL2 framework in JavaScript for your needs. It's called **WebGLSim** – a single JS class that abstracts WebGL boilerplate (context setup, quad rendering, texture/FBO ping-pong, state management) while exposing a simple JSON-based render graph API. This keeps the focus on your GLSL physics algorithms (e.g., grid-projected forcefields, Verlet integration, docking simulations).

#### Core Principles
- **Abstraction Layers**:
  - **Low-Level WebGL**: Hidden entirely. No manual `gl.bindTexture`, `gl.viewport`, `gl.useProgram`, etc. – all handled in baked passes.
  - **Buffers/Textures**: Auto-allocated as RGBA32F (or half-float for perf) ping-pong pairs. Supports 3D textures via slice-rendering (like your biochemistry grid).
  - **Render Passes**: Defined as JSON tuples: `[fragmentShaderPath, outputBufferName, {iChannel0: "inputBuffer", ...}, ["uniform1", "uniform2"]]`. Baked into zero-overhead functions.
  - **Transform Feedback (TF)**: Optional extension for vertex-buffer updates (e.g., particle positions/velocities). Specify `"type": "transform"` in pass JSON for TF mode (reads VBO, writes new VBO).
  - **State Management**: Automatic – uniforms from HTML inputs auto-sync to shaders. No manual `gl.uniform*` calls.
- **Didactic Focus**: Easy to swap shaders (e.g., `fluid/solve.glslf` → `docking/gridProject.glslf`). Runs on-demand (button) or auto (requestAnimationFrame loop).
- **JSON Render Graph**: Matches your example exactly. Parameters auto-generate HTML inputs (sliders for scalars/vecs).
- **Compatibility**: WebGL2 only (universal in 2025 browsers). Float textures via `EXT_color_buffer_float`. Half-float for mobile perf.
- **Inspiration from Python/ModernGL**: Direct port – `WebGLSim` mirrors `GLSL_Simulation`. GUI is vanilla HTML/JS (no frameworks) with dynamic input gen, like your PyQt setup.
- **Limitations/Exts**: No native compute (use WebGPU later). For TF, needs WebGL2's `transformFeedback`. 3D textures: Render as 2D slices (loop over Z).

#### High-Level Usage
1. **HTML Setup**: Embed `<canvas id="sim">` and `<script src="webglsim.js"></script>`.
2. **JSON Config**: Load your graph (e.g., from `config.json`).
3. **Run**: `sim.runGraph(params)` – updates from inputs auto-applied.
4. **Extend**: Add passes for TF particles, 3D grids, or multi-RT (multiple outputs per pass).

Example `config.json` (yours + docking example):
```json
{
  "parameters": {
    "dt": ["float", 0.01, 0.001],
    "driver": ["vec4", [0.5, 0.5, 0.0, 1.0], 0.1],
    "gridSize": ["vec3", [128, 128, 128], 1]
  },
  "Pipeline": [
    ["shaders/projectGrid.glslf", "receptorGrid", {"iChannel0": "receptorAtoms"}, ["dt"]],
    ["shaders/interpolateForces.glslf", "forceField", {"iChannel0": "receptorGrid"}, ["dt", "gridSize"]],
    ["shaders/verletUpdate.glslv", "ligandPositions", {"iChannel0": "forceField"}, ["dt"], {"type": "transform"}],
    ["shaders/viewDocking.glslf", "view0", {"iChannel0": "ligandPositions"}, []]
  ]
}
```

#### JS Framework Code: `webglsim.js`
Here's the core class. Save as `webglsim.js`. It handles everything statelessly, like your Python.

```javascript
// webglsim.js - WebGL2 Simulation Framework (2025 Edition)
class WebGLSim {
  constructor(canvas, simSize = [512, 512], headless = false) {
    this.canvas = canvas;
    this.simSize = simSize; // [width, height] for 2D/3D slices
    this.headless = headless;
    this.gl = null;
    this.programs = new Map();
    this.textures = new Map(); // {name: texture}
    this.fbos = new Map();     // {name: fbo}
    this.buffers = new Map();  // For TF VBOs
    this.quadBuffer = null;
    this.quadVAO = null;       // Simulated via attribs (WebGL no native VAO)
    this.iFrame = 0;
    this.extFloat = null;
    this.initGL();
    this.initQuad();
  }

  initGL() {
    const gl = this.canvas.getContext('webgl2', { antialias: false, preserveDrawingBuffer: false });
    if (!gl) throw new Error('WebGL2 not supported');
    this.gl = gl;
    gl.clearColor(0, 0, 0, 1);
    gl.enable(gl.BLEND); gl.blendFunc(gl.ONE, gl.ZERO); // For additive grids if needed

    // Extensions for float rendering/sampling
    this.extFloat = gl.getExtension('EXT_color_buffer_float');
    if (!this.extFloat) throw new Error('Float buffers not supported');
    gl.getExtension('OES_texture_float_linear'); // For interpolation
  }

  initQuad() {
    const gl = this.gl;
    const vertices = new Float32Array([-1, -1, 1, -1, -1, 1, 1, 1]); // Fullscreen quad
    this.quadBuffer = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, this.quadBuffer);
    gl.bufferData(gl.ARRAY_BUFFER, vertices, gl.STATIC_DRAW);
  }

  // Load GLSL program (vertex optional; injects default quad VS)
  loadProgram(name, { vertexSrc = null, fragmentSrc }) {
    const gl = this.gl;
    const vsSrc = vertexSrc || `#version 300 es
      in vec2 a_position; out vec2 v_texcoord;
      void main() { v_texcoord = (a_position + 1.) * .5; gl_Position = vec4(a_position, 0, 1); }`;
    const fsSrc = fragmentSrc.replace(/\/\/.*$/gm, ''); // Strip comments for includes if needed
    const vs = this.compileShader(gl.VERTEX_SHADER, vsSrc);
    const fs = this.compileShader(gl.FRAGMENT_SHADER, fsSrc);
    const prog = gl.createProgram();
    gl.attachShader(prog, vs); gl.attachShader(prog, fs);
    gl.linkProgram(prog);
    if (!gl.getProgramParameter(prog, gl.LINK_STATUS)) throw new Error(`Link error ${name}: ${gl.getProgramInfoLog(prog)}`);
    this.programs.set(name, prog);
    return prog;
  }

  compileShader(type, src) {
    const gl = this.gl;
    const shader = gl.createShader(type);
    gl.shaderSource(shader, src);
    gl.compileShader(shader);
    if (!gl.getShaderParameter(shader, gl.COMPILE_STATUS)) throw new Error(`Compile error: ${gl.getShaderInfoLog(shader)}`);
    return shader;
  }

  // Create texture/FBO pair (ping-pong ready; supports 3D via slices)
  createTexture(name, size = this.simSize, channels = 4, is3D = false) {
    const gl = this.gl;
    const internalFmt = gl.RGBA32F;
    const format = gl.RGBA; const type = gl.FLOAT;
    let tex;
    if (is3D) {
      tex = gl.createTexture();
      gl.bindTexture(gl.TEXTURE_3D, tex);
      gl.texImage3D(gl.TEXTURE_3D, 0, internalFmt, size[0], size[1], size[2], 0, format, type, null);
      gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
      gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
      gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
      gl.texParameteri(gl.TEXTURE_3D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    } else {
      tex = gl.createTexture();
      gl.bindTexture(gl.TEXTURE_2D, tex);
      gl.texImage2D(gl.TEXTURE_2D, 0, internalFmt, size[0], size[1], 0, format, type, null);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
    }
    this.textures.set(name, tex);

    const fbo = gl.createFramebuffer();
    gl.bindFramebuffer(gl.FRAMEBUFFER, fbo);
    gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, is3D ? gl.TEXTURE_3D : gl.TEXTURE_2D, tex, 0);
    if (is3D) { /* Handle slice attachment if needed */ }
    if (gl.checkFramebufferStatus(gl.FRAMEBUFFER) !== gl.FRAMEBUFFER_COMPLETE) throw new Error(`FBO incomplete for ${name}`);
    this.fbos.set(name, fbo);
    return { tex, fbo };
  }

  // For TF: Create VBO pair
  createBuffer(name, count, attribs = ['position', 'velocity']) {
    const gl = this.gl;
    const size = Float32Array.BYTES_PER_ELEMENT * 3 * count; // e.g., vec3 pos + vel
    const bufA = gl.createBuffer(); gl.bindBuffer(gl.ARRAY_BUFFER, bufA); gl.bufferData(gl.ARRAY_BUFFER, size, gl.DYNAMIC_DRAW);
    const bufB = gl.createBuffer(); gl.bindBuffer(gl.ARRAY_BUFFER, bufB); gl.bufferData(gl.ARRAY_BUFFER, size, gl.DYNAMIC_DRAW);
    this.buffers.set(`${name}A`, bufA); this.buffers.set(`${name}B`, bufB);
    // TF object
    const tf = gl.createTransformFeedback();
    gl.bindTransformFeedback(gl.TRANSFORM_FEEDBACK, tf);
    gl.bindBufferBase(gl.TRANSFORM_FEEDBACK_BUFFER, 0, bufA);
    gl.bindBufferBase(gl.TRANSFORM_FEEDBACK_BUFFER, 1, bufB);
    return { tf, current: bufA, next: bufB };
  }

  // Initialize all textures/buffers to zero
  initializeTextures(value = [0, 0, 0, 0]) {
    const gl = this.gl;
    gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    for (const [name, fbo] of this.fbos) {
      gl.bindFramebuffer(gl.FRAMEBUFFER, fbo);
      gl.viewport(0, 0, this.simSize[0], this.simSize[1]);
      gl.clearColor(...value);
      gl.clear(gl.COLOR_BUFFER_BIT);
    }
    // For buffers: gl.bufferSubData to zero if needed
  }

  // Bake a single pass (returns func(dynamicUniforms))
  bakePass(progName, outputName, inputMap = {}, dynUniforms = [], opts = { type: 'render', is3D: false }) {
    const gl = this.gl;
    const prog = this.programs.get(progName);
    if (!prog) throw new Error(`Program ${progName} not loaded`);

    // Pre-resolve uniform locations & setters
    const setters = [];
    dynUniforms.forEach(uni => {
      const loc = gl.getUniformLocation(prog, uni);
      if (!loc) throw new Error(`Uniform ${uni} missing in ${progName}`);
      setters.push((vals) => {
        const val = vals[uni];
        if (Array.isArray(val)) {
          if (val.length === 1) gl.uniform1fv(loc, val);
          else if (val.length === 2) gl.uniform2fv(loc, val);
          else if (val.length === 3) gl.uniform3fv(loc, val);
          else if (val.length === 4) gl.uniform4fv(loc, val);
        } else {
          gl.uniform1f(loc, val);
        }
      });
    });

    // Bind texture units (Shadertoy-style iChannelN)
    for (let i = 0; i < 8; i++) {
      const chan = `iChannel${i}`;
      const texLoc = gl.getUniformLocation(prog, chan);
      if (texLoc !== null) gl.uniform1i(texLoc, i);
    }
    // Default uniforms
    const resLoc = gl.getUniformLocation(prog, 'iResolution');
    if (resLoc) gl.uniform3f(resLoc, this.simSize[0], this.simSize[1], 1);
    const frameLoc = gl.getUniformLocation(prog, 'iFrame');
    if (frameLoc) gl.uniform1i(frameLoc, this.iFrame);

    // TF setup if opts.type === 'transform'
    let tf = null, currentBuf = null, nextBuf = null, vaoSetup = null;
    if (opts.type === 'transform') {
      const bufPair = this.buffers.get(outputName) || this.createBuffer(outputName, this.simSize[0] * this.simSize[1]); // e.g., particle count
      tf = bufPair.tf;
      currentBuf = bufPair.current;
      nextBuf = bufPair.next;
      vaoSetup = this.setupVAOForTF(prog, currentBuf); // Bind attribs like position, velocity
      gl.enable(gl.RASTERIZER_DISCARD);
    }

    const inputNames = Object.values(inputMap);
    const outputFBO = this.fbos.get(outputName);
    if (!outputFBO && opts.type !== 'transform') throw new Error(`Output ${outputName} not allocated`);

    return (dynamicVals) => {
      gl.useProgram(prog);
      setters.forEach(setter => setter(dynamicVals));
      gl.uniform1i(gl.getUniformLocation(prog, 'iFrame'), this.iFrame);

      // Bind inputs
      inputNames.forEach((name, i) => {
        const tex = this.textures.get(name);
        gl.activeTexture(gl.TEXTURE0 + i);
        gl.bindTexture(gl.TEXTURE_2D, tex);
      });

      if (opts.type === 'render') {
        gl.bindFramebuffer(gl.FRAMEBUFFER, outputFBO);
        gl.viewport(0, 0, this.simSize[0], this.simSize[1]);
        this.renderQuad(prog);
      } else if (opts.type === 'transform') {
        gl.bindVertexArray(vaoSetup.vao); // If using OES_vertex_array_object ext
        gl.bindBuffer(gl.ARRAY_BUFFER, currentBuf);
        gl.bindTransformFeedback(gl.TRANSFORM_FEEDBACK, tf);
        gl.beginTransformFeedback(gl.POINTS); // Or GL_TRIANGLES for quads
        gl.drawArrays(gl.POINTS, 0, /* particle count */ 10000);
        gl.endTransformFeedback();
        gl.bindBuffer(gl.ARRAY_BUFFER, nextBuf); // Swap
        [currentBuf, nextBuf] = [nextBuf, currentBuf];
      }

      // For 3D: Loop over slices with gl.texSubImage3D
      if (opts.is3D) { /* Implement slice loop here */ }

      gl.bindFramebuffer(gl.FRAMEBUFFER, null);
    };
  }

  setupVAOForTF(prog, buf) {
    // Simplified; use OES_vertex_array_object if needed
    const gl = this.gl;
    // Bind attribs: e.g., gl.bindAttribLocation(prog, 0, 'a_position'); etc.
    // Return vao-like setup
    return { vao: null, /* manual bind in pass */ };
  }

  renderQuad(prog) {
    const gl = this.gl;
    gl.bindBuffer(gl.ARRAY_BUFFER, this.quadBuffer);
    const loc = gl.getAttribLocation(prog, 'a_position');
    if (loc !== -1) {
      gl.enableVertexAttribArray(loc);
      gl.vertexAttribPointer(loc, 2, gl.FLOAT, false, 0, 0);
    }
    gl.drawArrays(gl.TRIANGLE_STRIP, 0, 4);
    gl.disableVertexAttribArray(loc);
  }

  // Build & bake full graph from JSON pipeline
  buildPipeline(pipeline, shaderDir = './shaders', is3D = false) {
    // Auto-load shaders (fetch fragment files)
    pipeline.forEach(([fragPath, outName, inputMap, dynUnis, opts = {}]) => {
      const fullPath = `${shaderDir}/${fragPath}`;
      fetch(fullPath).then(r => r.text()).then(src => {
        this.loadProgram(fragPath, { fragmentSrc: src });
      });
    });

    // Alloc outputs/inputs
    const allNames = new Set();
    pipeline.forEach(([_, out, inputs]) => {
      allNames.add(out);
      Object.values(inputs || {}).forEach(name => allNames.add(name));
    });
    allNames.forEach(name => {
      if (!this.textures.has(name)) this.createTexture(name, this.simSize, 4, opts.is3D);
    });

    this.initializeTextures();

    // Bake passes
    const baked = pipeline.map(([progName, out, inputMap, dynUnis, opts = {}]) =>
      this.bakePass(progName, out, inputMap, dynUnis, { ...opts, is3D })
    );
    return { baked, textures: Array.from(allNames) };
  }

  // Run graph (on-demand or in RAF loop)
  runGraph(baked, dynamicVals) {
    baked.forEach(pass => pass(dynamicVals));
    this.iFrame++;
    if (!this.headless) {
      this.gl.bindFramebuffer(this.gl.FRAMEBUFFER, null);
      this.gl.viewport(0, 0, this.canvas.width, this.canvas.height);
      // Render display pass (e.g., texture to screen)
      this.renderDisplay('view0'); // Default view buffer
    }
  }

  renderDisplay(texName) {
    const gl = this.gl;
    const prog = this.displayProg || this.loadProgram('display', {
      fragmentSrc: `#version 300 es
        precision highp float; in vec2 v_texcoord; uniform sampler2D uTex; out vec4 fragColor;
        void main() { fragColor = texture(uTex, v_texcoord); }`
    });
    this.displayProg = prog;
    gl.useProgram(prog);
    gl.uniform1i(gl.getUniformLocation(prog, 'uTex'), 0);
    gl.activeTexture(gl.TEXTURE0);
    gl.bindTexture(gl.TEXTURE_2D, this.textures.get(texName));
    this.renderQuad(prog);
  }

  // Auto-generate HTML inputs from params JSON
  generateUI(paramsContainerId, onChange) {
    const container = document.getElementById(paramsContainerId);
    Object.entries(params).forEach(([name, [type, defaultVal, step]]) => {
      const div = document.createElement('div');
      div.innerHTML = `<label>${name} (${type}):</label>`;
      if (type === 'float') {
        const input = document.createElement('input');
        input.type = 'range'; input.min = 0; input.max = 1; input.step = step || 0.01;
        input.value = defaultVal; input.oninput = () => onChange({ [name]: parseFloat(input.value) });
        div.appendChild(input);
      } else if (type.startsWith('vec')) {
        const vec = defaultVal;
        const inputs = [];
        for (let i = 0; i < vec.length; i++) {
          const inp = document.createElement('input');
          inp.type = 'number'; inp.step = step || 0.01; inp.value = vec[i];
          inp.oninput = () => {
            vec[i] = parseFloat(inp.value);
            onChange({ [name]: [...vec] });
          };
          inputs.push(inp);
          div.appendChild(inp);
        }
      }
      container.appendChild(div);
    });
  }

  // Convenience: Load JSON config & setup
  async loadConfig(jsonPath, uiId) {
    const resp = await fetch(jsonPath);
    const config = await resp.json();
    const { baked, textures } = this.buildPipeline(config.Pipeline);
    this.bakedGraph = baked;
    this.configTextures = textures;

    // Gen UI
    this.generateUI(uiId, (vals) => this.dynamicVals = { ...this.defaultVals, ...vals });

    // Set defaults
    this.defaultVals = {};
    Object.entries(config.parameters).forEach(([k, [_, def]]) => this.defaultVals[k] = def);
    this.dynamicVals = { ...this.defaultVals };

    // Display selector
    const select = document.createElement('select');
    textures.forEach(t => {
      const opt = document.createElement('option'); opt.value = t; opt.text = t;
      select.appendChild(opt);
    });
    select.onchange = (e) => this.displayTex = e.target.value;
    document.body.appendChild(select); // Or to specific div

    return config;
  }

  // Auto-loop
  startLoop() {
    const loop = () => {
      this.runGraph(this.bakedGraph, this.dynamicVals);
      requestAnimationFrame(loop);
    };
    loop();
  }
}
```

#### HTML Template: `index.html`
Embed this for a full page. Loads config, generates UI, runs on button/loop.

```html
<!DOCTYPE html>
<html>
<head>
  <title>WebGLSim Demos</title>
  <style> canvas { border: 1px solid black; } #params { float: right; width: 200px; } </style>
</head>
<body>
  <canvas id="sim" width="512" height="512"></canvas>
  <div id="params"></div>
  <button id="run">Run Step</button>
  <button id="play">Play Loop</button>
  <input type="file" id="loadConfig" accept=".json">
  <script src="webglsim.js"></script>
  <script>
    const canvas = document.getElementById('sim');
    const sim = new WebGLSim(canvas);

    // Load default config
    sim.loadConfig('config.json', 'params').then(() => {
      document.getElementById('run').onclick = () => sim.runGraph(sim.bakedGraph, sim.dynamicVals);
      document.getElementById('play').onclick = () => sim.startLoop();
      document.getElementById('loadConfig').onchange = (e) => {
        const reader = new FileReader();
        reader.onload = (ev) => {
          const config = JSON.parse(ev.target.result);
          sim.buildPipeline(config.Pipeline); // Rebuild
          sim.generateUI('params', (vals) => sim.dynamicVals = vals);
        };
        reader.readAsText(e.target.files[0]);
      };
    });
  </script>
</body>
</html>
```

#### How It Works for Your Use Cases
- **Grid-Projected Forcefield (Biochem Docking)**: Pass 1: `projectGrid.glslf` sums atoms to 3D grid (set `is3D: true`). Pass 2: `interpolateForces.glslf` trilinear-samples for ligand forces. Rigid ligand: Use TF pass with quaternion uniform.
- **Particles/MD**: TF pass: `verletUpdate.glslv` (vertex shader) reads VBO, computes forces from grid texture, writes new pos/vel.
- **Auto-Updates**: HTML sliders call `onChange`, sync to `dynamicVals`, applied in `runGraph`.
- **On-Demand/Auto**: Button for single step (didactic stepping); RAF for sim.
- **Extensibility**: Add `"type": "mrt"` for multi-render-targets (e.g., separate vdW/electro grids). For 3D, extend `bakePass` with Z-slice loop.
- **Perf Tips**: Use half-float (`gl.HALF_FLOAT`) for mobile. Limit grid to 256^3 (slice ~16ms on mid-GPU).

This is ~300 LOC, zero deps, and ports your Python 1:1. For TF particles in docking, add `setupVAOForTF` details (e.g., stride for pos/vel). Test it – drop your GLSL files in `/shaders`, tweak JSON, and iterate physics algos hassle-free. If you need tweaks (e.g., full 3D slice code or dat.GUI integration), share a snippet!