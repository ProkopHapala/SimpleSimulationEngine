# === AUTO-DOC BEGIN ===
# @file EulerianImpacFluid.py
# @brief PyOpenCL wrapper for 2D Eulerian multi-material compressible flow solver with level-set interface tracking.
#
# Thin Python layer over OpenCL kernels in EulerianImpacFluid.cl. Manages buffer
# allocation, ping-pong state swapping, and host-device data transfer. The solver
# implements a 5-equation diffuse-interface model with stiffened-gas EOS for
# high-velocity impact simulations (e.g. uranium projectile into liquid hydrogen).
#
# **Key design**: U (rho_a1, rho_a2, ru, rv) packed as float4, E (energy, phi) as float2.
# Two buffer pairs for ping-pong — no read-after-write hazards. NVIDIA GPU preferred
# via _pick_nvidia_gpu() to avoid interactive context prompts on multi-GPU systems.
# === AUTO-DOC END ===

import os
import numpy as np
import pyopencl as cl
import math


def _pick_nvidia_gpu():
    """Select NVIDIA GPU OpenCL context, falling back to any GPU or default device."""
    plats = cl.get_platforms()
    for p in plats:
        if "nvidia" in p.name.lower() or "cuda" in p.name.lower():
            gpus = [d for d in p.get_devices() if d.type & cl.device_type.GPU]
            if gpus:
                return cl.Context(devices=[gpus[0]])
    # Fallback: first GPU, else first device
    for p in plats:
        gpus = [d for d in p.get_devices() if d.type & cl.device_type.GPU]
        if gpus:
            return cl.Context(devices=[gpus[0]])
    return cl.create_some_context(interactive=False)


class EulerianImpactFluid:
    """OpenCL-accelerated 2D Eulerian fluid solver with level-set multi-material interface tracking."""

    def __init__(self, width, height, dx):
        self.width = np.int32(width)
        self.height = np.int32(height)
        self.dx = np.float32(dx)
        
        # OpenCL Setup: prefer NVIDIA GPU, avoid interactive prompts
        self.ctx = _pick_nvidia_gpu()
        self.queue = cl.CommandQueue(self.ctx)
        
        # Load kernels
        cl_path = os.path.join(os.path.dirname(__file__), "EulerianImpacFluid.cl")
        with open(cl_path, "r") as f:
            self.program = cl.Program(self.ctx, f.read()).build()
        # Cache kernels to avoid repeated retrieval warnings
        self.k_update_fluid = cl.Kernel(self.program, "update_fluid")
        self.k_redistance_phi = cl.Kernel(self.program, "redistance_phi")
            
        self.allocate_buffers()

    def allocate_buffers(self):
        """Allocate ping-pong U/E buffers plus diagnostic pressure, sound speed, and volume accumulator."""
        mf = cl.mem_flags
        # in_U_vec: float4(rho_a1, rho_a2, ru, rv)
        # in_E_vec: float2(E, phi)
        self.buf_U1 = cl.Buffer(self.ctx, mf.READ_WRITE, size=self.width * self.height * 16) # float4
        self.buf_U2 = cl.Buffer(self.ctx, mf.READ_WRITE, size=self.width * self.height * 16)
        self.buf_E1 = cl.Buffer(self.ctx, mf.READ_WRITE, size=self.width * self.height * 8)  # float2
        self.buf_E2 = cl.Buffer(self.ctx, mf.READ_WRITE, size=self.width * self.height * 8)
        self.buf_E_orig = cl.Buffer(self.ctx, mf.READ_WRITE, size=self.width * self.height * 8)
        # RK2: buffers to save U^n for corrector stage
        self.buf_U_old = cl.Buffer(self.ctx, mf.READ_ONLY, size=self.width * self.height * 16)
        self.buf_E_old = cl.Buffer(self.ctx, mf.READ_ONLY, size=self.width * self.height * 8)
        # Pressure diagnostic buffer (float per cell)
        self.buf_P = cl.Buffer(self.ctx, mf.READ_WRITE, size=self.width * self.height * 4)  # [Diagnostic] Pressure
        # Sound speed diagnostic buffer (float per cell)
        self.buf_C = cl.Buffer(self.ctx, mf.READ_WRITE, size=self.width * self.height * 4)  # [Diagnostic] Sound speed
        
        # For volume calculation / mass conservation
        self.buf_vol = cl.Buffer(self.ctx, mf.READ_WRITE, size=4) # int (atomic accumulator)

    def initialize(self, rho_a1, rho_a2, ru, rv, E, phi):
        """Upload initial conserved fields to GPU, mirroring to both ping-pong buffers to avoid uninitialized reads on first swap."""
        # host arrays should be numpy arrays of appropriate shape and type
        mf = cl.mem_flags
        
        # Pack data into float4 and float2
        U_host = np.zeros((self.height, self.width, 4), dtype=np.float32)
        U_host[..., 0] = rho_a1
        U_host[..., 1] = rho_a2
        U_host[..., 2] = ru
        U_host[..., 3] = rv
        
        E_host = np.zeros((self.height, self.width, 2), dtype=np.float32)
        E_host[..., 0] = E
        E_host[..., 1] = phi
        
        cl.enqueue_copy(self.queue, self.buf_U1, U_host.flatten())
        cl.enqueue_copy(self.queue, self.buf_E1, E_host.flatten())
        # Mirror to ping-pong buffers to avoid reading uninitialized data on first swap
        cl.enqueue_copy(self.queue, self.buf_U2, U_host.flatten())
        cl.enqueue_copy(self.queue, self.buf_E2, E_host.flatten())
        cl.enqueue_copy(self.queue, self.buf_E_orig, E_host.flatten())
        self.queue.finish()

    def step(self, dt, b_diagnostics=True):
        """Advance one timestep via SSP-RK2 (predictor + corrector), each using MUSCL 2nd-order flux.

        SSP-RK2 (Shu-Osher):
          Stage 0 (predictor):  U* = U^n + dt * L(U^n)
          Stage 1 (corrector): U^{n+1} = 0.5*U^n + 0.5*(U* + dt*L(U*))

        Buffer flow:
          Before: buf_U1/buf_E1 = U^n
          Save U^n → buf_U_old/buf_E_old
          Stage 0: read buf_U1 → write buf_U2 (U*). Swap: buf_U1=U*, buf_U2=free
          Stage 1: read buf_U1(U*) + buf_U_old(U^n) → write buf_U2 (U^{n+1}). Swap.
          After: buf_U1/buf_E1 = U^{n+1}
        """
        dt = np.float32(dt)
        b_diag = np.int32(1 if b_diagnostics else 0)
        stage0 = np.int32(0)
        stage1 = np.int32(1)

        tiles_x = math.ceil(int(self.width) / 8)
        tiles_y = math.ceil(int(self.height) / 8)
        global_size = (tiles_x * 64, tiles_y)
        local_size = (64, 1)

        # Save U^n for corrector
        cl.enqueue_copy(self.queue, self.buf_U_old, self.buf_U1)
        cl.enqueue_copy(self.queue, self.buf_E_old, self.buf_E1)

        # Stage 0 (predictor): U* = U^n + dt*L(U^n)
        # Read buf_U1 (U^n), write buf_U2 (U*)
        self.k_update_fluid.set_args(
            self.buf_U1, self.buf_E1,
            self.buf_U2, self.buf_E2, self.buf_P, self.buf_C,
            self.buf_U_old, self.buf_E_old,
            self.width, self.height, self.dx, dt, b_diag, stage0
        )
        cl.enqueue_nd_range_kernel(self.queue, self.k_update_fluid, global_size, local_size)
        # Swap: buf_U1 = U*, buf_U2 = old (will be overwritten)
        self.buf_U1, self.buf_U2 = self.buf_U2, self.buf_U1
        self.buf_E1, self.buf_E2 = self.buf_E2, self.buf_E1

        # Stage 1 (corrector): U^{n+1} = 0.5*U^n + 0.5*(U* + dt*L(U*))
        # Read buf_U1 (U*) + buf_U_old (U^n), write buf_U2 (U^{n+1})
        self.k_update_fluid.set_args(
            self.buf_U1, self.buf_E1,
            self.buf_U2, self.buf_E2, self.buf_P, self.buf_C,
            self.buf_U_old, self.buf_E_old,
            self.width, self.height, self.dx, dt, b_diag, stage1
        )
        cl.enqueue_nd_range_kernel(self.queue, self.k_update_fluid, global_size, local_size)
        # Swap: buf_U1 = U^{n+1}
        self.buf_U1, self.buf_U2 = self.buf_U2, self.buf_U1
        self.buf_E1, self.buf_E2 = self.buf_E2, self.buf_E1
        
    def redistance(self, iterations=5):
        """Reinitialize phi to signed-distance property via iterated redistance_phi kernel — prevents interface smearing."""
        dtau = np.float32(0.5 * self.dx)
        
        # Cache original state for sign function
        cl.enqueue_copy(self.queue, self.buf_E_orig, self.buf_E1)
        
        tiles_x = math.ceil(int(self.width) / 8)
        tiles_y = math.ceil(int(self.height) / 8)
        global_size = (tiles_x * 64, tiles_y)
        local_size = (64, 1)
        
        for _ in range(iterations):
            self.k_redistance_phi.set_args(
                self.buf_E1, self.buf_E_orig, self.buf_E2,
                self.width, self.height, self.dx, dtau
            )
            cl.enqueue_nd_range_kernel(self.queue, self.k_redistance_phi, global_size, local_size)
            self.buf_E1, self.buf_E2 = self.buf_E2, self.buf_E1
            
    def get_data(self):
        """Read back all fields (rho_a1, rho_a2, ru, rv, E, phi, pressure, sound speed) as a dict of 2D numpy arrays."""
        U_host = np.empty((self.height, self.width, 4), dtype=np.float32)
        E_host = np.empty((self.height, self.width, 2), dtype=np.float32)
        P_host = np.empty((self.height, self.width), dtype=np.float32)
        C_host = np.empty((self.height, self.width), dtype=np.float32)
        
        cl.enqueue_copy(self.queue, U_host, self.buf_U1)
        cl.enqueue_copy(self.queue, E_host, self.buf_E1)
        cl.enqueue_copy(self.queue, P_host, self.buf_P)
        cl.enqueue_copy(self.queue, C_host, self.buf_C)
        self.queue.finish()
        
        return {
            'rho_a1': U_host[..., 0],
            'rho_a2': U_host[..., 1],
            'ru': U_host[..., 2],
            'rv': U_host[..., 3],
            'E': E_host[..., 0],
            'phi': E_host[..., 1],
            'p': P_host,
            'c': C_host
        }

    def get_max_wave_speed(self):
        """Read max wave speed (|u|+c) from diagnostic buffers for CFL-based dt computation."""
        U_host = np.empty((self.height, self.width, 4), dtype=np.float32)
        C_host = np.empty((self.height, self.width), dtype=np.float32)
        cl.enqueue_copy(self.queue, U_host, self.buf_U1)
        cl.enqueue_copy(self.queue, C_host, self.buf_C)
        self.queue.finish()
        rho = U_host[..., 0] + U_host[..., 1]
        mask = rho > 1e-4
        u = np.where(mask, U_host[..., 2] / (rho + 1e-6), 0.0)
        v = np.where(mask, U_host[..., 3] / (rho + 1e-6), 0.0)
        vel_mag = np.sqrt(u**2 + v**2)
        return float((vel_mag + C_host).max())
