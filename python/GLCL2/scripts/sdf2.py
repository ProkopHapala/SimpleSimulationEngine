import numpy as np

# SDF2 demo: initialize a small RGBA32F texture (8x8) in Python and view it
config = {
    "simulation_name":  "GLCL2 SDF2 (init texture)",
    "description":      "Initialize a tiny float4 texture in Python (no gen pass) and render it.",

    # Parameters
    # TW,TH: texture size
    "parameters": {
        "TW":      (8,    "int",   1),
        "TH":      (8,    "int",   1),
        # driver kept for future zoom/pan; not used by current view shader path
        "driver":  ([0.5, 0.5, 16.0, 0.0], "vec4", 0.1),
    },

    # No OpenCL for this demo
    "buffers": {},
    "opencl_source": [],
    "kernels": {},
    "kernel_pipeline": [],

    # Shaders: fullscreen quad VS + SDF view FS program (core profile)
    "opengl_shaders": {
        "viewSDF": ("../shaders/fs_quad.glslv", "../shaders/sdf/view_core.glslf", ["iChannel0", "iResolution", "iFrame", "driver"]),
    },

    # One small RGBA32F texture to display
    "textures_2d": {
        "sdfTex": ("TW", "TH"),
    },

    # No FBOs needed (we render directly to default framebuffer)
    "fbos": {},

    # FS pipeline: directly view the pre-initialized texture
    # Format: (shader_name, out_fbo_name or 'default', [input_texture_names])
    "fs_pipeline": [
        ("viewSDF", "default", ["sdfTex"]),
    ],

    # No traditional render pipeline
    "render_pipeline": [],
}


def init():
    # Build 8x8 RGBA float32 image with a few colored pixels; rest black
    W = config["parameters"]["TW"][0]
    H = config["parameters"]["TH"][0]
    img = np.zeros((H, W, 4), dtype=np.float32)
    v0=-0.1
    img[:,:,:] = (v0, v0, v0, 1.0)

    # A few colored pixels
    # pts = [
    #     (0, 0, (1.0, 1.0, 1.0, 1.0)),  # white
    #     (1, 1, (1.0, 0.0, 0.0, 1.0)),  # red
    #     (2, 2, (0.0, 1.0, 0.0, 1.0)),  # green
    #     (3, 3, (0.0, 0.0, 1.0, 1.0)),  # blue
    #     (4, 1, (1.0, 1.0, 0.0, 1.0)),  # yellow
    #     (6, 5, (1.0, 0.0, 1.0, 1.0)),  # magenta
    # ]
    # for x, y, c in pts:
    #     if 0 <= x < W and 0 <= y < H:
    #         img[y, x, :] = c

    img[ 1:5,3,   :] = (1.0, 1.0, 0.0, 1.0);
    img[ 1:3,3:5, :] = (1.0, 1.0, 0.0, 1.0);
    img[ 1:4,2:4, :] = (1.0, 0.0, 0.0, 1.0);

    # Return dict with the special key for FS texture uploads
    return {
        "__textures_2d__": {
            "sdfTex": img
        }
    }
