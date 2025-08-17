import numpy as np

config = {
    "simulation_name":  "Soldiers 2D (GLCL2)",
    "description":      "Boids-like 2D soldiers with lateral line formation and enemy-facing dynamics.",

    # GUI parameters (name: (value, type, step))
    "parameters": {
        "particle_count": (20,   "int",   1      ),
        "formation_count":(2,    "int",   1      ),
        "formation_point_count": (4, "int", 1),
        "dt":             (0.01, "float", 0.0005 ),
        # neighborhoods / spacings
        "r_cut":          (0.25, "float", 0.005  ),
        "d_front":        (0.20, "float", 0.005  ),
        "d_side":         (0.06, "float", 0.002  ),
        # weights
        "w_sep":          (1.00, "float", 0.02   ),
        "w_align":        (0.80, "float", 0.02   ),
        "w_coh":          (0.30, "float", 0.02   ),
        "w_enemy":        (1.20, "float", 0.02   ),
        # kinematics
        "max_speed":      (0.50, "float", 0.01   ),
        "max_turn":       (2.00, "float", 0.05   ),  # rad/s
        "friction":       (0.50, "float", 0.01   ),
        # anisotropy shaping
        "side_bias":      (1.50, "float", 0.05   ),
        "front_bias":     (0.50, "float", 0.05   ),
        # desire to formations
        "desire_gain":    (1.00, "float", 0.02   ),
        "w_goal":         (1.00, "float", 0.02   ),
        "beta":           (0.20, "float", 0.02   ),
        "max_accel":      (3.00, "float", 0.05   ),
    },

    # Buffers: (size, stride, type)
    "buffers": {
        "state_pos_dir":   ("particle_count", 4, "f4"),
        "state_vel_mr":    ("particle_count", 4, "f4"),
        "state_team_type": ("particle_count", 4, "i4"),
        # per-formation: (o.x,o.y,d.x,d.y,w_par,w_perp,k_par,k_perp)
        "formation_params": ("formation_count", 8, "f4"),
        # for visualization: two endpoints per formation as points (px,py,ux,uy) using same VS
        "formation_points": ("formation_point_count", 4, "f4"),
        # interleaved line vertices for formations: pos.xy + color.rgba
        "formation_lines": ("formation_point_count", 6, "f4"),
        # local quad geometry (4 vertices, 2 components) for instanced rendering
        "quad_verts": (4, 2, "f4"),
    },

    # OpenCL
    "opencl_source": ["../cl/soldiers.cl"],
    "kernels": {
        #            local_size  global_size        buffers                                                           parameters (in kernel order)
        "soldiers_step": (None, ("particle_count"), ["state_pos_dir","state_vel_mr","state_team_type","formation_params"], [
            "dt","particle_count","r_cut","d_front","d_side","w_sep","w_align","w_coh","w_enemy","max_speed","max_turn","friction","side_bias","front_bias",
            "desire_gain","w_goal","beta","max_accel"
        ])
    },
    "kernel_pipeline": ["soldiers_step"],

    # OpenGL shaders
    "opengl_shaders": {
        "soldier_render": ("../shaders/soldier_points.glslv", "../shaders/monocolor.glslf", ["color"]),  # color is set by viewer
        "line_color":     ("../shaders/line_color.glslv",     "../shaders/color_passthrough.glslf", []),
        # Instanced rectangle for formations: uses formation_params split into two vec4 attributes [pose, misc]
        #"pose2d_rect":    ("../shaders/pose2d_rect.glslv",    "../shaders/monocolor.glslf", ["color"]),
        # NOTE: soldier-oriented quads (instanced_quad) kept for later; disabled in render_pipeline for clarity
        "instanced_quad": ("../shaders/instanced_quad.glslv", "../shaders/monocolor.glslf", ["color"]),
    },

    # Rendering pipeline: (shader, element_count, vertex_buffer, index_buffer)
    # NOTE: Soldier-oriented instanced quads are temporarily disabled to avoid confusion while we introduce
    #       formation instanced rectangles driven by formation_params. Soldiers still render as points.
    "render_pipeline": [
        ("line_color",    "formation_point_count", "formation_lines", None, {"mode":"LINES", "attribs":[2,4]}),
        # Formation rectangles: one instance per formation, geometry from quad_verts
        ("pose2d_rect",  "formation_count",       "quad_verts",       None, {"mode":"TRIANGLE_STRIP", "attribs":[2], "instance_buffer":"formation_params", "instance_attribs":[4,4], "instance_divisor":1}),
        # Soldier-oriented quads (instanced) using state_pos_dir as (px,py,ux,uy); draw before points so points overlay as fallback
        ("instanced_quad", "particle_count",       "quad_verts",       None, {"mode":"TRIANGLE_STRIP", "attribs":[2], "instance_buffer":"state_pos_dir", "instance_attribs":[4], "instance_divisor":1}),
        ("soldier_render", "formation_point_count", "formation_points", None),
        ("soldier_render", "particle_count", "state_pos_dir", None),
    ],
}


def init():
    n  = config["parameters"]["particle_count"][0]
    nf = config["parameters"]["formation_count"][0]
    n0 = n // 2
    n1 = n - n0

    # Buffers
    pos_dir   = np.zeros((n, 4), dtype=np.float32)
    vel_mr    = np.zeros((n, 4), dtype=np.float32)
    team_type = np.zeros((n, 4), dtype=np.int32)
    fpars     = np.zeros((nf, 8), dtype=np.float32)
    fpoints   = np.zeros((nf*2, 4), dtype=np.float32)
    flines    = np.zeros((nf*2, 6), dtype=np.float32)
    quad      = np.array([[-1.0,-1.0],[1.0,-1.0],[-1.0,1.0],[1.0,1.0]], dtype=np.float32)

    # Team 0: y = -0.4, facing +x
    if n0 > 0:
        x0 = np.linspace(-0.6, -0.1, n0, dtype=np.float32)
        y0 = np.full(n0, -0.4, dtype=np.float32)
        jitter0 = (np.random.rand(n0,2).astype(np.float32) - 0.5) * 0.02
        pos_dir[:n0, 0] = x0 + jitter0[:,0]
        pos_dir[:n0, 1] = y0 + jitter0[:,1]
        pos_dir[:n0, 2] = 1.0  # ux
        pos_dir[:n0, 3] = 0.0  # uy
        team_type[:n0, 0] = 0    # team 0
        team_type[:n0, 1] = 0    # type id 0
        team_type[:n0, 2] = 0    # formation 0

    # Team 1: y = +0.4, facing -x
    if n1 > 0:
        x1 = np.linspace(0.1, 0.6, n1, dtype=np.float32)
        y1 = np.full(n1,  0.4, dtype=np.float32)
        jitter1 = (np.random.rand(n1,2).astype(np.float32) - 0.5) * 0.02
        pos_dir[n0:, 0] = x1 + jitter1[:,0]
        pos_dir[n0:, 1] = y1 + jitter1[:,1]
        pos_dir[n0:, 2] = -1.0  # ux
        pos_dir[n0:, 3] =  0.0  # uy
        team_type[n0:, 0] = 1    # team 1
        team_type[n0:, 1] = 0    # type id 0
        team_type[n0:, 2] = 1    # formation 1

    # Velocities, mass, radius
    vel_mr[:, 0:2] = 0.0
    vel_mr[:, 2]   = 1.0   # mass
    vel_mr[:, 3]   = 0.015 # radius

    # Formation params and points
    # fpar = (o.x,o.y,d.x,d.y,w_par,w_perp,k_par,k_perp)
    # Formation 0: line along +x near y=-0.4
    if nf > 0:
        fpars[0,0:2] = np.array([-0.8, -0.4], dtype=np.float32)
        fpars[0,2:4] = np.array([ 1.0,  0.0], dtype=np.float32)
        fpars[0,4:8] = np.array([0.6, 0.08, 0.5, 2.0], dtype=np.float32)
        # endpoints for visualization
        L = 1.0
        o = fpars[0,0:2]; d = fpars[0,2:4] / np.linalg.norm(fpars[0,2:4])
        fpoints[0,0:2] = o - L*d
        fpoints[1,0:2] = o + L*d
        # line vertices with color (red)
        flines[0,0:2] = fpoints[0,0:2]; flines[0,2:6] = np.array([1.0,0.0,0.0,0.9], dtype=np.float32)
        flines[1,0:2] = fpoints[1,0:2]; flines[1,2:6] = np.array([1.0,0.0,0.0,0.9], dtype=np.float32)
    if nf > 1:
        fpars[1,0:2] = np.array([ 0.8,  0.4], dtype=np.float32)
        fpars[1,2:4] = np.array([-1.0,  0.0], dtype=np.float32)
        fpars[1,4:8] = np.array([0.6, 0.08, 0.5, 2.0], dtype=np.float32)
        L = 1.0
        o = fpars[1,0:2]; d = fpars[1,2:4] / np.linalg.norm(fpars[1,2:4])
        fpoints[2,0:2] = o - L*d
        fpoints[3,0:2] = o + L*d
        # line vertices with color (cyan)
        flines[2,0:2] = fpoints[2,0:2]; flines[2,2:6] = np.array([0.0,1.0,1.0,0.9], dtype=np.float32)
        flines[3,0:2] = fpoints[3,0:2]; flines[3,2:6] = np.array([0.0,1.0,1.0,0.9], dtype=np.float32)

    return {
        "state_pos_dir":   pos_dir,
        "state_vel_mr":    vel_mr,
        "state_team_type": team_type,
        "formation_params": fpars,
        "formation_points": fpoints,
        "formation_lines":  flines,
        "quad_verts":       quad,
    }
