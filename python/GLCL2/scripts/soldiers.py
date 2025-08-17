import numpy as np

config = {
    "simulation_name":  "Soldiers 2D (GLCL2)",
    "description":      "Boids-like 2D soldiers with lateral line formation and enemy-facing dynamics.",

    # GUI parameters (name: (value, type, step))
    "parameters": {
        "particle_count": (20,   "int",   1      ),
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
    },

    # Buffers: (size, stride, type)
    "buffers": {
        "state_pos_dir":   ("particle_count", 4, "f4"),
        "state_vel_mr":    ("particle_count", 4, "f4"),
        "state_team_type": ("particle_count", 4, "f4"),
    },

    # OpenCL
    "opencl_source": ["../cl/soldiers.cl"],
    "kernels": {
        #            local_size  global_size        buffers                                   parameters (in kernel order)
        "soldiers_step": (None, ("particle_count"), ["state_pos_dir","state_vel_mr","state_team_type"], [
            "dt","particle_count","r_cut","d_front","d_side","w_sep","w_align","w_coh","w_enemy","max_speed","max_turn","friction","side_bias","front_bias"
        ])
    },
    "kernel_pipeline": ["soldiers_step"],

    # OpenGL shaders
    "opengl_shaders": {
        "soldier_render": ("../shaders/soldier_points.glslv", "../shaders/monocolor.glslf", ["color"])  # color is set by viewer
    },

    # Rendering pipeline: (shader, element_count, vertex_buffer, index_buffer)
    "render_pipeline": [
        ("soldier_render", "particle_count", "state_pos_dir", None),
    ],
}


def init():
    n = config["parameters"]["particle_count"][0]
    n0 = n // 2
    n1 = n - n0

    # Buffers
    pos_dir   = np.zeros((n, 4), dtype=np.float32)
    vel_mr    = np.zeros((n, 4), dtype=np.float32)
    team_type = np.zeros((n, 4), dtype=np.float32)

    # Team 0: y = -0.4, facing +x
    if n0 > 0:
        x0 = np.linspace(-0.6, -0.1, n0, dtype=np.float32)
        y0 = np.full(n0, -0.4, dtype=np.float32)
        jitter0 = (np.random.rand(n0,2).astype(np.float32) - 0.5) * 0.02
        pos_dir[:n0, 0] = x0 + jitter0[:,0]
        pos_dir[:n0, 1] = y0 + jitter0[:,1]
        pos_dir[:n0, 2] = 1.0  # ux
        pos_dir[:n0, 3] = 0.0  # uy
        team_type[:n0, 0] = 0.0  # team 0
        team_type[:n0, 1] = 0.0  # type id 0

    # Team 1: y = +0.4, facing -x
    if n1 > 0:
        x1 = np.linspace(0.1, 0.6, n1, dtype=np.float32)
        y1 = np.full(n1,  0.4, dtype=np.float32)
        jitter1 = (np.random.rand(n1,2).astype(np.float32) - 0.5) * 0.02
        pos_dir[n0:, 0] = x1 + jitter1[:,0]
        pos_dir[n0:, 1] = y1 + jitter1[:,1]
        pos_dir[n0:, 2] = -1.0  # ux
        pos_dir[n0:, 3] =  0.0  # uy
        team_type[n0:, 0] = 1.0  # team 1
        team_type[n0:, 1] = 0.0  # type id 0

    # Velocities, mass, radius
    vel_mr[:, 0:2] = 0.0
    vel_mr[:, 2]   = 1.0   # mass
    vel_mr[:, 3]   = 0.015 # radius

    return {
        "state_pos_dir":   pos_dir,
        "state_vel_mr":    vel_mr,
        "state_team_type": team_type,
    }
