import numpy as np
import pyopencl as cl
import os

class GPUTerrainSolver:
    def __init__(self, width, height):
        self.width = np.int32(width)
        self.height = np.int32(height)
        self.ctx = cl.create_some_context()
        self.queue = cl.CommandQueue(self.ctx)
        
        # Load Kernel
        current_dir = os.path.dirname(__file__) if __file__ != "" else "."
        with open(os.path.join(current_dir, "terrain.cl"), "r") as f:
            self.prg = cl.Program(self.ctx, f.read()).build()

        self.mf = cl.mem_flags
        self.num_tiles_x = width // 16
        self.num_tiles_y = height // 16

    def solve(self, h_map, K=0.01):
        # Buffers
        h_gpu = cl.Buffer(self.ctx, self.mf.READ_ONLY | self.mf.COPY_HOST_PTR, hostbuf=h_map.astype(np.float32))
        parent_gpu = cl.Buffer(self.ctx, self.mf.READ_WRITE, 4 * self.width * self.height)
        cost_gpu = cl.Buffer(self.ctx, self.mf.READ_WRITE, 4 * self.width * self.height)
        sinks_gpu = cl.Buffer(self.ctx, self.mf.READ_WRITE, 4 * self.num_tiles_x * self.num_tiles_y)

        # Run local solver
        self.prg.solve_tiles(self.queue, (self.width, self.height), (16, 16),
                            h_gpu, parent_gpu, cost_gpu, sinks_gpu,
                            self.width, self.height, np.float32(K))

        # Retrieve
        parent_res = np.empty((self.height, self.width), dtype=np.uint32)
        cost_res = np.empty((self.height, self.width), dtype=np.float32)
        sinks_res = np.empty(self.num_tiles_y * self.num_tiles_x, dtype=np.uint32)

        cl.enqueue_copy(self.queue, parent_res, parent_gpu)
        cl.enqueue_copy(self.queue, cost_res, cost_gpu)
        cl.enqueue_copy(self.queue, sinks_res, sinks_gpu)

        return parent_res, cost_res, sinks_res

def unpack(val):
    return val & 0xFFFF, (val >> 16) & 0xFFFF

def pack(x, y):
    return np.uint32(x | (y << 16))