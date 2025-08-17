# Simple FS-pipeline demo that mimics pySymGLSL fluid ping-pong
# Runs a sequence of fullscreen passes writing to textures fieldA,B,C and then displays

config = {
    "simulation_name":  "GLCL2 Fluid (FS ping-pong)",
    "description":      "Reimplementation of pySymGLSL fluid in GLCL2 using FS ping-pong passes.",

    # Parameters used to size textures and control driver force
    #             value, type , step
    "parameters": {
        #                  value, type, step
        "W":        (512,   "int",   None),
        "H":        (512,   "int",   None),
        # driver = (pos.x,pos.y, force.x, force.y) in normalized coords
        # NOTE: The reference pySymGLSL uses unit-scale forces; previous tiny values led to a static field.
        "driver":   ([0.5, 0.5, 0.0, 10.0], "vec4", 0.1),  # DEBUG: stronger force for visible motion
    },

    # No OpenCL buffers/kernels for this demo
    "buffers": {},
    "opencl_source": [],
    "kernels": {},
    "kernel_pipeline": [],

    # OpenGL shaders: FS quad VS + fluid solver/view FS programs (core profile)
    "opengl_shaders": {
        "solveFluid": ("../shaders/fs_quad.glslv", "../shaders/fluid_solve.glslf", ["iChannel0", "iResolution", "iFrame", "driver"]),
        "viewFluid":  ("../shaders/fs_quad.glslv", "../shaders/fluid_view.glslf",      ["iChannel0", "iResolution", "iFrame"]),
    },

    # Textures and FBOs
    "textures_2d": {
        "fieldA": ("W", "H"),
        "fieldB": ("W", "H"),
        "fieldC": ("W", "H"),
    },
    "fbos": {
        "fboA": "fieldA",
        "fboB": "fieldB",
        "fboC": "fieldC",
    },

    # FS pipeline: 3x solve passes ping-ponging A->B->C->A, then visualize C
    # Tuple format: (shader_name, out_fbo_name or 'default', [input_texture_names])
    "fs_pipeline": [
        ("solveFluid", "fboB", ["fieldA"]),
        ("solveFluid", "fboC", ["fieldB"]),
        ("solveFluid", "fboA", ["fieldC"]),
        ("viewFluid",  "default", ["fieldC"]),
    ],

    # No traditional render pipeline (we present via FS)
    "render_pipeline": [],
}

# No init() needed because we have no CL/GL buffers for this demo
