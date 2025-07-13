# pySimGLSL: Physics Simulation with GLSL and Render Graphs

This project aims to build a flexible and performant framework for physics simulations, particularly focusing on grid-based methods and partial differential equations (PDEs), by leveraging the power of GLSL shaders and GPU acceleration. The core idea is to treat simulation states as 2D images (textures) and use GLSL shaders to perform computations on these textures, creating a feedback loop for time-stepping the simulation.

## Core Concepts:

1.  **Grid-Based Simulations:** Many physics simulations, such as fluid dynamics, heat transfer, or wave propagation, can be discretized on a grid. Each grid cell's state (e.g., velocity, pressure, temperature) can be mapped to a pixel in a 2D texture.

2.  **GLSL Shaders for Computation:** The GPU's highly parallel architecture is ideal for performing computations on large grids. GLSL (OpenGL Shading Language) fragment shaders are used to implement the simulation logic. Each pixel (fragment) processed by the shader corresponds to a grid cell, and the shader calculates the new state of that cell based on its current state and the states of its neighbors.

3.  **Texture Feedback Loop:** To advance the simulation in time, a feedback loop is employed. The output of a shader (the new state) is rendered to a framebuffer object (FBO) which stores the result as a texture. This output texture then becomes an input texture for the next simulation step, effectively passing the state from one frame to the next. This typically involves ping-ponging between two textures/FBOs to avoid reading and writing to the same texture simultaneously.

4.  **Render Graphs (JSON Configuration):** The sequence of simulation steps, including which shaders to run, which textures to use as inputs and outputs, and what uniforms to set, is defined as a render graph. This graph is configured using JSON files, allowing for flexible and easily modifiable simulation pipelines without recompiling code. This enables the creation of complex multi-pass simulations.

5.  **Mapping to 2D Images:** The natural mapping of grid data to 2D textures (images) makes GLSL a powerful tool for these types of simulations. Visualizing the simulation state is also straightforward as it's inherently an image.

## Project Structure (Anticipated): 

-   `GLSL_GUI.py`: Handles the graphical user interface, displaying the simulation output, and potentially providing controls for parameters.
-   `GLSL_Simulation.py`: Manages the OpenGL context, framebuffer objects, textures, and the execution of the GLSL shaders based on the render graph.
-   `shaders/`: Directory containing GLSL shader files (`.glslf`). These implement the actual physics equations.
-   `pipelines/`: Directory containing JSON files that define the simulation pipelines/render graphs, specifying the order of shader execution, input/output textures, and uniform variables.

## Texture Sampling and Parameters

### Key Differences from OpenCL
- In GLSL, sampler parameters are set on the texture object from host code (Python)
- Cannot declare sampler parameters in shader like OpenCL's `__constant sampler_t`
- All textures share same sampling parameters (no multiple samplers per texture)

### Available Parameters (set via moderngl.Texture):
- `filter`: Tuple of (min_filter, mag_filter) - NEAREST or LINEAR
- `repeat_x/y`: Boolean - True=REPEAT, False=CLAMP_TO_EDGE
- `anisotropy`: Float - Typically 1.0 to 16.0

### Current Defaults in Our System:
```python
tex.filter = (moderngl.NEAREST, moderngl.NEAREST)
tex.repeat_x = False  # CLAMP_TO_EDGE
tex.repeat_y = False  # CLAMP_TO_EDGE
tex.anisotropy = 1.0  # No anisotropic filtering
```

### Coordinate System:
- Uses normalized coordinates [0,1] by default
- Handled by vertex shader: `v_texcoord = (in_position + 1.0) * 0.5`

## Current Issue (Investigated):

There is an observed issue where the output image of the fluid simulation is strangely split into four parts, suggesting a discontinuity or replication. This is likely due to incorrect texture coordinates or sampler parameters (e.g., `GL_TEXTURE_WRAP_S`, `GL_TEXTURE_WRAP_T`, `GL_TEXTURE_MIN_FILTER`, `GL_TEXTURE_MAG_FILTER`) within the GLSL feedback loop, leading to misaligned sampling or rendering.

**Investigation Focus:**

-   How are textures created and configured (e.g., `glTexParameteri`) in `GLSL_Simulation.py`?
-   How are texture coordinates generated and used in `solveFluid.glslf`?
-   How does `fluid.json` define the texture inputs and outputs for each pass?
-   Are there any transformations or scaling applied to the viewport or textures that could cause this splitting?

## Resolution of Initial Issue:

The initial issue of the output image being split into four parts was identified as an incorrect `uv` coordinate calculation in `solveFluid.glslf`. The `uv` coordinates were being calculated as `gl_FragCoord.xy/iResolution.xy - 0.5`, which shifted the sampling range to `[-0.5, 0.5]`. This, combined with the default `CLAMP_TO_EDGE` texture wrapping and `NEAREST` filtering, led to the observed replication. The fix involved changing the calculation to `gl_FragCoord.xy/iResolution.xy`, ensuring `uv` coordinates remain within the `[0,1]` range for proper texture sampling.
