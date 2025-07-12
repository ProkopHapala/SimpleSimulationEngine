"""pySymGLSL – Lightweight convenience layer around moderngl for GPGPU style
physics simulations using the traditional rendering pipeline (ping-pong
framebuffers, multiple render passes).

The public API is intentionally minimal.  Everything is exposed through the
`Simulation` class which lives in `simulation.py`.

Typical usage::

    from pySymGLSL import Simulation

    sim = Simulation(sim_size=(512, 512))
    sim.create_texture("velocity")          # creates velocity_a & velocity_b
    sim.load_program("advect", fragment_path="shaders/advect.glsl")
    # ... bake render graph and execute ...
"""

from .GLSL_Simulation import GLSL_Simulation

__all__ = ["GLSL_Simulation"]
