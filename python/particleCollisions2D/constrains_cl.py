import numpy as np
import pyopencl as cl
import os

def select_platform_device():
    """Select OpenCL platform and device based on environment variables or defaults"""
    # Check environment variables
    platform_pref = os.environ.get('OCL_PLATFORM', 'NVIDIA')  # prefer NVIDIA by default
    device_pref   = os.environ.get('OCL_DEVICE', '0')        # first device by default
    
    # Get all platforms
    platforms = cl.get_platforms()
    
    # Find preferred platform
    platform = None
    for p in platforms:
        if platform_pref.lower() in p.name.lower():
            platform = p
            break
    if platform is None:
        platform = platforms[0]  # fallback to first platform
    
    # Get devices for the platform
    devices = platform.get_devices()
    device_idx = int(device_pref)
    if device_idx >= len(devices):
        device_idx = 0
    device = devices[device_idx]
    
    print(f"Using OpenCL platform: {platform.name}")
    print(f"Using device: {device.name}")
    
    return cl.Context([device])

class CLConstrains:
    def __init__(self, max_points=1000, max_neighs=32, local_size=32):
        # OpenCL setup
        self.ctx = select_platform_device()
        self.queue = cl.CommandQueue(self.ctx)
        self.local_size = local_size
        
        # Load kernel
        kernel_path = os.path.join(os.path.dirname(__file__), 'constrains.cl')
        with open(kernel_path, 'r') as f:
            kernel_src = f.read()
        self.prg = cl.Program(self.ctx, kernel_src).build()
        
        # Initialize parameters
        self.max_points = max_points
        self.max_neighs = max_neighs
        
        # Initialize buffers
        self.allocate_buffers()
    
    def allocate_buffers(self):
        mf = cl.mem_flags
        
        # Main particle data (float4: x,y,z,mass)
        self.ps_buf = cl.Buffer(self.ctx, mf.READ_WRITE, size=self.max_points * 4 * 4)
        self.ps_out_buf = cl.Buffer(self.ctx, mf.READ_WRITE, size=self.max_points * 4 * 4)
        
        # Neighbor data
        total_neighs = self.max_points * self.max_neighs
        self.neighs_buf = cl.Buffer(self.ctx, mf.READ_WRITE, size=total_neighs * 4)
        self.kngs_buf = cl.Buffer(self.ctx, mf.READ_WRITE, size=total_neighs * 4)
        self.l0s_buf = cl.Buffer(self.ctx, mf.READ_WRITE, size=total_neighs * 4)
    
    def updateJacobi_neighs(self, pos2D, masses, neighs, kngs, l0s, inv_dt2, Rd):
        """
        Run one iteration of Jacobi solver with direct data conversion and upload
        
        Args:
            pos2D: numpy array of shape (n,2) with positions
            masses: numpy array of shape (n,) with masses
            neighs: list of lists with neighbor indices
            kngs: list of lists with neighbor stiffnesses
            l0s: list of lists with rest lengths
            inv_dt2: inverse of dt squared
            Rd: relaxation distance
        Returns:
            Updated 2D positions
        """
        n = len(pos2D)
        assert n <= self.max_points, f"Too many particles: {n} > {self.max_points}"
        
        # Convert positions and masses to float4
        pos_mass = np.zeros((n, 4), dtype=np.float32)
        pos_mass[:, :2] = pos2D.astype(np.float32)
        pos_mass[:, 3] = masses.astype(np.float32)
        
        # Convert neighbor data to arrays
        neighs_arr = np.full((n, self.max_neighs), -1, dtype=np.int32)
        kngs_arr = np.zeros((n, self.max_neighs), dtype=np.float32)
        l0s_arr = np.zeros((n, self.max_neighs), dtype=np.float32)
        
        for i in range(n):
            ni = len(neighs[i])
            if ni > self.max_neighs:
                raise ValueError(f"Too many neighbors for particle {i}: {ni} > {self.max_neighs}")
            neighs_arr[i,:ni] = neighs[i]
            kngs_arr[i,:ni] = kngs[i]
            l0s_arr[i,:ni] = l0s[i]
        
        # Upload data
        cl.enqueue_copy(self.queue, self.ps_buf, pos_mass)
        cl.enqueue_copy(self.queue, self.neighs_buf, neighs_arr)
        cl.enqueue_copy(self.queue, self.kngs_buf, kngs_arr)
        cl.enqueue_copy(self.queue, self.l0s_buf, l0s_arr)
        
        # Run kernel
        self.prg.updateJacobi_neighs(
            self.queue, (n,), None,
            np.int32(n),
            np.int32(self.max_neighs),
            self.ps_buf,
            self.ps_out_buf,
            self.neighs_buf,
            self.kngs_buf,
            self.l0s_buf,
            np.float32(inv_dt2),
            np.float32(Rd)
        )
        
        # Read back results
        pos_mass = np.zeros((n, 4), dtype=np.float32)
        cl.enqueue_copy(self.queue, pos_mass, self.ps_out_buf)
        return pos_mass[:, :2].copy()  # return only x,y coordinates
    
    def solve_constraints(self, pos2D, masses, neighs, kngs, l0s, inv_dt2, Rd, niter=10):
        """Run multiple iterations of Jacobi solver"""
        pos = pos2D
        for _ in range(niter):
            pos = self.updateJacobi_neighs(pos, masses, neighs, kngs, l0s, inv_dt2, Rd)
        return pos

if __name__ == "__main__":
    # Initialize
    solver = CLConstrains(max_points=1000, max_neighs=32)

    # Update particle data
    pos2D = np.random.rand(100, 2).astype(np.float32)
    masses = np.ones(100).astype(np.float32)
    neighs = [[1, 2], [0, 2], [0, 1]] + [[-1] * 32 for _ in range(97)]
    kngs = [[1.0, 1.0], [1.0, 1.0], [1.0, 1.0]] + [[0.0] * 32 for _ in range(97)]
    l0s = [[1.0, 1.0], [1.0, 1.0], [1.0, 1.0]] + [[0.0] * 32 for _ in range(97)]

    # Run simulation step
    inv_dt2 = 1.0
    Rd = 0.05
    niter = 10
    final_positions = solver.solve_constraints(pos2D, masses, neighs, kngs, l0s, inv_dt2, Rd, niter)