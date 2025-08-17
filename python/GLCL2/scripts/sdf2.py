import numpy as np
import os, sys

# SDF2 demo: initialize a small RGBA32F texture (8x8) in Python and view it
config = {
    "simulation_name":  "GLCL2 SDF2 (init texture)",
    "description":      "Initialize a tiny float4 texture in Python (no gen pass) and render it.",

    # Parameters
    # TW,TH: texture size
    "parameters": {
        #                  value, type, step
        "TW":      (10,    "int",   None),
        "TH":      (10,    "int",   None),
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
        "viewSDF": ("../shaders/fs_quad.glslv", "../shaders/sdf_view.glslf", ["iChannel0", "iResolution", "iFrame", "driver"]),
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
    #from ..subPixelContour import subPixelContour as spc  # type: ignore
    import os, sys, importlib.util
    

    this_dir = os.path.dirname(__file__)
    root     = os.path.abspath(os.path.join(this_dir, '..', '..'))
    spc_path = os.path.join(root, 'subPixelContour', 'subPixelContour.py')
    print( "!!!!!!!!! spc_path ", spc_path)
    spec = importlib.util.spec_from_file_location('subPixelContour', spc_path)
    spc_module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(spc_module)
    spc = spc_module
    if spc is None:
        raise Exception("ERROR: Failed to import subPixelContour, spc is None")

    # Build low-res signed field via subPixelContour bilinear fit and upload as RGBA32F
    W = int(config["parameters"]["TW"][0])
    H = int(config["parameters"]["TH"][0])

    # Robust import of subPixelContour regardless of loader context
    # DEBUG: do not overwrite the dynamically imported module object
    # spc = None

    # Basis nodes on integer grid 0..W-1 x 0..H-1
    basis = spc.generate_basis_points(nx=W, ny=H)

    # Dense sample grid across cell centers domain [-0.5, W-0.5] x [-0.5, H-0.5]
    res = max(64, 8 * max(W, H))
    xmin, xmax = -0.5, (W - 0.5)
    ymin, ymax = -0.5, (H - 0.5)
    _, _, pts = spc.generate_grid_points(xmin, xmax, ymin, ymax, res=res)

    # Example shapes
    circles = [
        (W*0.35, H*0.35, 1.8),
        (W*0.65, H*0.65, 1.2),
    ]
    polygon = np.array([
        [0.0,   0.0],
        [W-2.0, 1.0],
        [W-2.0, H-2.0],
        [1.0,   H-2.0],
    ], dtype=float)  # CCW

    b = spc.evaluate_reference(pts, circles=circles, polygon=polygon)  # +/-1 values
    opt_log = {}
    # Build bilinear A and fit coefficients with clamped error descent
    A = spc.build_A(pts, basis, kind='bilinear')
    c = spc.fit_coeffs(A, b, learning_rate=0.1, n_iterations=1000, log=opt_log, ftol=1e-5, )
    C = c.reshape(H, W).astype(np.float32)

    #from matplotlib import pyplot as plt; spc.plot_opt_log(opt_log); plt.show()

    # Pack into RGBA32F; shader reads red channel. Set alpha=1.
    img = np.zeros((H, W, 4), dtype=np.float32)
    img[:, :, 0] = C
    img[:, :, 3] = 1.0

    # v0=-0.1
    # img[:,:,:] = (v0, v0, v0, 1.0)
    # img[ 1:5,3,   :] = (1.0, 1.0, 0.0, 1.0);
    # img[ 1:3,3:5, :] = (1.0, 1.0, 0.0, 1.0);
    # img[ 1:4,2:4, :] = (1.0, 0.0, 0.0, 1.0);

    # Return dict with the special key for FS texture uploads
    return { "__textures_2d__": { "sdfTex": img } }
