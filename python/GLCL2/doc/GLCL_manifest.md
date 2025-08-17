
# GLCL Framework Documentation

## 1. Overview and Goal

The GLCL framework is a powerful, yet easy-to-use tool for developing and visualizing scientific simulations in Python. Inspired by [Shadertoy](https://www.shadertoy.com/), its primary goal is to abstract away the boilerplate code associated with creating a GUI, managing OpenGL/OpenCL contexts, and handling data buffers.

This allows users—researchers, students, and hobbyists—to focus purely on their core algorithms, which are expressed through three simple components:

1.  **OpenCL kernels** (`.cl` files) for high-performance parallel computation on the GPU.
2.  **GLSL shaders** (`.glslv`, `.glslf` files) for visualizing the results of the computation.
3.  A central **Python configuration script** (e.g., `nbody.py`) that declaratively wires everything together.

The `GLCLBrowser` dynamically loads these user-provided components, providing an interactive environment for rapid prototyping, exploration, and visualization of complex physical systems.

## 2. Core Design Principles

The framework's architecture is guided by several key principles to ensure it is powerful, efficient, and easy to extend.

*   **Generality:** The browser is designed to be completely agnostic about the simulation it is running. It makes no assumptions about the names of buffers (`positions`, `velocities`), uniforms (`dt`), kernels (`nbody_sim`), or shaders. All simulation-specific logic and names are defined entirely within the user's Python script. This allows the same browser to run an N-body simulation, a fluid dynamics simulation, or a reaction-diffusion system without any changes to the framework code.

*   **Baking:** To achieve high performance, the framework uses a "baking" architecture. Instead of interpreting the configuration dictionary and looking up objects by name in every frame, the browser processes the configuration *once* at startup. It builds optimized, pre-packaged execution pipelines for both computation and rendering. This eliminates runtime string comparisons, dictionary lookups, and other overhead from the main simulation loop.

*   **Low Runtime Cost:** The main simulation loop (`update_simulation`) is designed to be as lean as possible. It consists of three simple steps:
    1.  Execute a list of pre-baked OpenCL kernel calls.
    2.  Transfer only the necessary data from OpenCL to OpenGL buffers.
    3.  Trigger a repaint of the OpenGL widget.
    All complex logic is handled during the initial baking phase, not at runtime.

*   **Simplicity and Modularity:** The system is divided into clear, logical components (`OCLSystem`, `OGLSystem`, `GLCLWidget`, `GLCLBrowser`). Each class has a distinct responsibility, making the codebase easier to navigate, debug, and maintain.

## 3. Usage Example

You can run a simulation by pointing the `GLCLBrowser` to a user-defined Python script from the command line.

```bash
# Basic usage, loads the nbody.py script by default
python -m pyBall.GLCL.GLCLBrowser

# Specify a different script
python -m pyBall.GLCL.GLCLBrowser --script path/to/my_simulation.py

# Run with OpenCL debugging (prints buffer data each frame)
python -m pyBall.GLCL.GLCLBrowser --debug-cl --debug-frames 10

# Run with OpenGL debugging (skips kernel execution to test rendering)
python -m pyBall.GLCL.GLCLBrowser --debug-gl
```

## 4. The User Simulation Script (`nbody.py`)

The user defines their entire simulation in a Python script containing a `config` dictionary and an `init()` function.

```python
import numpy as np

# The central configuration dictionary
config = {
    # --- General Info ---
    "simulation_name":  "NBody Simulation",
    "description":      "...",

    # --- Simulation Parameters (Tunable via GUI) ---
    "parameters": {
        # name:           (default_value, type,    step_size)
        "particle_count": (2048,          "int" ,   1      ),
        "dt":             (0.001,         "float",  0.0001 )
    },

    # --- Data Buffers (for OpenCL and/or OpenGL) ---
    "buffers": {
        # name:         (size_expression, stride, numpy_dtype)
        "positions":    ("particle_count", 4,      "f4"),
        "velocities":   ("particle_count", 4,      "f4")
    },

    # --- OpenCL Computation ---
    "opencl_source": ["../cl/nbody.cl"],
    "kernels": {
        # kernel_name: (local_size, global_size_expr, [buffer_args], [scalar_args])
        "nbody_sim" :  ((32,),      "particle_count", ["positions", "velocities"], ["dt", "particle_count"])
    },
    "kernel_pipeline": ["nbody_sim"],

    # --- OpenGL Visualization ---
    "opengl_shaders": {
        # shader_name:   (vertex_shader_path,  fragment_shader_path, [uniforms])
        "nbody_render" : ("../shaders/points.glslv", "../shaders/monocolor.glslf",  ["color"])
    },
    "render_pipeline":   [
        # (shader_name, element_count_expr, vertex_buffer_name, index_buffer_name)
        ("nbody_render", "particle_count",   "positions",        None),
    ]
}

# Function to provide initial buffer data
def init():
    n = config["parameters"]["particle_count"][0]
    pos = (np.random.rand(n, 4) * 2 - 1).astype(np.float32)
    vel = np.zeros_like(pos)
    return {
        "positions": pos,
        "velocities": vel
    }
```

### Format Explanation

*   `parameters`: Defines variables that will be exposed in the GUI as tunable controls. The browser automatically creates a `QSpinBox` for `"int"` types and a `QDoubleSpinBox` for `"float"` types.
*   `buffers`: Declares all data buffers used by the simulation. The `size_expression` can be an integer or a string matching a name from `parameters`.
*   `opencl_source`: A list of `.cl` files to load. All kernels within these files will be compiled and made available.
*   `kernels`: Defines the properties for each OpenCL kernel you intend to use. This includes work sizes (`local_size`, `global_size`) and the ordered lists of buffer and scalar arguments that will be passed to it.
*   `kernel_pipeline`: An ordered list of kernel names to execute sequentially in each simulation step.
*   `opengl_shaders`: Defines named shader programs by linking vertex and fragment shader source files.
*   `render_pipeline`: An ordered list of draw passes to execute for rendering. Each tuple specifies which shader to use, how many elements to draw, and which buffer contains the vertex data.
*   `init()`: A required function that returns a dictionary. The keys must match buffer names from the `buffers` section, and the values must be NumPy arrays containing the initial data for those buffers.

## 5. Technical Implementation (For Developers)

This section provides an overview of the internal architecture for developers looking to debug or extend the framework.

### 5.1 Key Components

*   **`GLCLBrowser.py`**: The central coordinator. It manages the main application window, GUI controls, and the simulation lifecycle (loading, pausing, resetting). It is responsible for orchestrating the "baking" process and running the main simulation loop.
*   **`GLCLWidget.py`**: The rendering view. This `QOpenGLWidget` is responsible for all visualization. It receives a pre-baked render configuration and executes the specified draw passes in its `paintGL` method. It also handles camera controls and other user interactions within the 3D view.
*   **`OCLSystem.py`**: A manager class that encapsulates all PyOpenCL interactions. It handles context creation, program compilation, buffer allocation, and kernel execution. It maintains a global cache of all compiled kernels, which can be retrieved by name.
*   **`OGLSystem.py`**: A manager class for OpenGL resources, primarily shader programs. It is responsible for compiling and linking GLSL shaders and providing them to the `GLCLWidget`.
*   **`BakedKernelCall`**: A simple data class that holds a fully prepared OpenCL kernel for execution. It contains the kernel object itself, the work sizes, and a tuple of pre-resolved arguments (e.g., `cl.Buffer` objects and correctly typed scalars), eliminating any lookup overhead at runtime.

### 5.2 The Application Lifecycle

The framework follows a distinct load-bake-run lifecycle.

#### Phase 1: Initialization and Script Loading (`load_and_apply_script`)

1.  When a script is loaded, the `GLCLBrowser` first clears any existing configuration from the `OCLSystem` and `OGLSystem`.
2.  It uses `importlib` to dynamically load the user's Python script as a module.
3.  It extracts the `config` dictionary and `init` function from the module.
4.  It triggers the main configuration application process.

#### Phase 2: The "Baking" Process (`apply_simulation_config`)

This is the core of the setup phase where the declarative `config` is transformed into imperative, optimized execution pipelines.

1.  **Create GUI**: The `parameters` section is parsed to create the corresponding Qt widgets for user control.
2.  **Load CL Sources**: `OCLSystem` is instructed to load and compile all `.cl` files listed in `opencl_source`. All kernels are compiled and cached globally.
3.  **Initialize Buffers**:
    *   The `init()` function is called to get the initial NumPy arrays.
    *   The `buffers` section is iterated. For each buffer, an OpenCL buffer is created in `OCLSystem` with the correct size, and the initial data from `init()` is immediately copied to the GPU.
    *   The initial data is also stored in the `GLCLBrowser` to be passed to the `GLCLWidget`.
4.  **Bake Kernels**:
    *   The `kernel_pipeline` is iterated. For each kernel name:
    *   The corresponding kernel object is fetched from `OCLSystem`'s cache.
    *   All buffer arguments are resolved to `pyopencl.Buffer` objects.
    *   All scalar arguments are resolved to their correctly typed values (e.g., `np.int32`, `np.float32`).
    *   A `BakedKernelCall` object is instantiated with these pre-resolved objects and stored in `self.baked_kernel_calls`.
5.  **Prepare OpenGL**: The `opengl_shaders` and `render_pipeline` configurations are passed to the `OGLSystem` and `GLCLWidget`, respectively. **No OpenGL objects are created yet**, as the context may not be valid.
6.  **Pre-compute Sync List**: The browser finds the intersection of buffer names used in the `render_pipeline` and those used in the `kernels`. This creates a `self.buffers_to_sync` list, so the runtime loop knows exactly which buffers need to be transferred from OpenCL to OpenGL without checking every frame.

#### Phase 3: OpenGL Context Initialization (`GLCLWidget.initializeGL`)

This method is called by PyQt5 only when a valid OpenGL context has been created. This is the **critical and only safe place** to initialize OpenGL resources.

1.  **Compile Shaders**: `OGLSystem.compile_all_shaders()` is called. It iterates through the stored configurations and compiles/links all GLSL shader programs.
2.  **Bake Render Objects**: `GLCLWidget.bake_render_objects()` is called. It iterates through the buffers used in the `render_pipeline`, creates a `GLobject` (containing a VAO and VBO) for each one, and uploads the initial data.

#### Phase 4: The Runtime Loop (`update_simulation`)

The main loop, triggered by a `QTimer`, is now extremely simple and efficient.

1.  **Execute Kernels**: It iterates through the `self.baked_kernel_calls` list and calls `execute()` on each object.
2.  **Sync Buffers**: It iterates through the pre-computed `self.buffers_to_sync` list. For each buffer name, it downloads the data from OpenCL and uploads it to the corresponding `GLobject`'s VBO in the `GLCLWidget`.
3.  **Request Repaint**: It calls `self.glcl_widget.update()`, which schedules a call to `paintGL`.

### 5.3 Debugging Guide

*   Use the `--debug-cl` flag to print the first element of every OpenCL buffer after each kernel execution. This is useful for checking if your computations are producing the expected results (e.g., non-`NaN` values).
*   Use the `--debug-gl` flag to skip all OpenCL kernel execution. This allows you to verify that your initial data (from `init()`) and your GLSL shaders are set up correctly to render the initial state.
*   The console output is verbose during the "Baking" phase. Read it carefully to see if all your files are being found and if all resources (buffers, kernels, shaders) are being created as you expect.
*   If the application crashes on startup, it is often due to an error in the user script. The most common causes are a misnamed buffer/parameter, incorrect argument lists for kernels, or a syntax error in a shader or kernel file. The console log will typically contain a traceback pointing to the issue.

### 5.4 Fullscreen shader uniforms (FS pipeline)

- __Defaults__
  - `iResolution` is set once per program from the active FS target size in `GLCLWidget.apply_fs_static_uniforms()`.
  - `iFrame` is updated per frame inside `GLCLWidget._execute_fs_pipeline()`.

- __Dynamic uniforms from config__
  - For each FS program listed in `config["opengl_shaders"]`, declared uniform names are taken from the third tuple element (uniform list).
  - For every declared name (excluding defaults and sampler aliases), the current value is looked up in `config["parameters"]` by name and applied based on type (`int`, `float`, `vec2`, `vec3`, `vec4`).

- __Samplers__
  - Inputs of each FS pass bind to texture units in order and are exposed under `iChannelN` (and aliases `texN`, `sN`). These are not expected in `parameters`.

- __Performance__
  - Uniform locations are cached per program; uniforms are applied on GL init and on parameter changes to avoid per-frame overhead.