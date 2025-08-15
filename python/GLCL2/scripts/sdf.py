# Minimal SDF demo: generate a small RGBA32F SDF texture (8x8) and zoom it to screen with smooth edges

config = {
    "simulation_name":  "GLCL2 SDF (tiny RGBA32F)",
    "description":      "Generate a tiny float4 SDF texture and render it zoomed with smooth AA edges.",

    # Parameters
    # TW,TH: texture size; driver = (cx,cy, zoom, radius)
    "parameters": {
        "TW":      (8,    "int",   1),
        "TH":      (8,    "int",   1),
        "driver":  ([0.5, 0.5, 16.0, 0.35], "vec4", 0.1),
    },

    # No OpenCL for this demo
    "buffers": {},
    "opencl_source": [],
    "kernels": {},
    "kernel_pipeline": [],

    # Shaders: fullscreen quad VS + SDF gen/view FS programs (core profile)
    "opengl_shaders": {
        "genSDF":  ("../shaders/fs_quad.glslv", "../shaders/sdf/gen_core.glslf",  ["iResolution", "iFrame", "driver"]),
        "viewSDF": ("../shaders/fs_quad.glslv", "../shaders/sdf/view_core.glslf", ["iChannel0", "iResolution", "iFrame", "driver"]),
    },

    # One small RGBA32F texture to store the SDF
    "textures_2d": {
        "sdfTex": ("TW", "TH"),
    },
    # FBO that renders into sdfTex
    "fbos": {
        "fboSDF": "sdfTex",
    },

    # FS pipeline: generate SDF into sdfTex, then view it on the default framebuffer
    # Format: (shader_name, out_fbo_name or 'default', [input_texture_names])
    "fs_pipeline": [
        ("genSDF",  "fboSDF", []),
        ("viewSDF", "default", ["sdfTex"]),
    ],

    # No traditional render pipeline (present via FS pipeline)
    "render_pipeline": [],
}
