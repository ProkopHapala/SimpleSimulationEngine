import os
import numpy as np
import pyopencl as cl
import math


def _pick_nvidia_gpu():
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
        mf = cl.mem_flags
        # in_U_vec: float4(rho_a1, rho_a2, ru, rv)
        # in_E_vec: float2(E, phi)
        self.buf_U1 = cl.Buffer(self.ctx, mf.READ_WRITE, size=self.width * self.height * 16) # float4
        self.buf_U2 = cl.Buffer(self.ctx, mf.READ_WRITE, size=self.width * self.height * 16)
        self.buf_E1 = cl.Buffer(self.ctx, mf.READ_WRITE, size=self.width * self.height * 8)  # float2
        self.buf_E2 = cl.Buffer(self.ctx, mf.READ_WRITE, size=self.width * self.height * 8)
        self.buf_E_orig = cl.Buffer(self.ctx, mf.READ_WRITE, size=self.width * self.height * 8)
        # Pressure diagnostic buffer (float per cell)
        self.buf_P = cl.Buffer(self.ctx, mf.READ_WRITE, size=self.width * self.height * 4)  # [Diagnostic] Pressure
        # Sound speed diagnostic buffer (float per cell)
        self.buf_C = cl.Buffer(self.ctx, mf.READ_WRITE, size=self.width * self.height * 4)  # [Diagnostic] Sound speed
        
        # For volume calculation / mass conservation
        self.buf_vol = cl.Buffer(self.ctx, mf.READ_WRITE, size=4) # int (atomic accumulator)

    def initialize(self, rho_a1, rho_a2, ru, rv, E, phi):
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

    def step(self, dt):
        dt = np.float32(dt)
        
        # 2D NDRange: x-dim covers tiles (local 32 threads), y-dim covers tile rows
        tiles_x = math.ceil(int(self.width) / 8)
        tiles_y = math.ceil(int(self.height) / 8)
        global_size = (tiles_x * 32, tiles_y)
        local_size = (32, 1)
        
        # Update Fluid
        self.k_update_fluid.set_args(
            self.buf_U1, self.buf_E1,
            self.buf_U2, self.buf_E2, self.buf_P, self.buf_C,
            self.width, self.height, self.dx, dt
        )
        cl.enqueue_nd_range_kernel(self.queue, self.k_update_fluid, global_size, local_size)
        
        # Ping-pong swap
        self.buf_U1, self.buf_U2 = self.buf_U2, self.buf_U1
        self.buf_E1, self.buf_E2 = self.buf_E2, self.buf_E1
        
    def redistance(self, iterations=5):
        dtau = np.float32(0.5 * self.dx)
        
        # Cache original state for sign function
        cl.enqueue_copy(self.queue, self.buf_E_orig, self.buf_E1)
        
        tiles_x = math.ceil(int(self.width) / 8)
        tiles_y = math.ceil(int(self.height) / 8)
        global_size = (tiles_x * 32, tiles_y)
        local_size = (32, 1)
        
        for _ in range(iterations):
            self.k_redistance_phi.set_args(
                self.buf_E1, self.buf_E_orig, self.buf_E2,
                self.width, self.height, self.dx, dtau
            )
            cl.enqueue_nd_range_kernel(self.queue, self.k_redistance_phi, global_size, local_size)
            self.buf_E1, self.buf_E2 = self.buf_E2, self.buf_E1
            
    def get_data(self):
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
