import numpy as np
import pyopencl as cl
import os

class CLConstrains:
    def __init__(self, max_points=1000, max_neighs=32, local_size=32):
        # OpenCL setup
        self.ctx = cl.create_some_context()
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
        
        # Main particle data
        self.ps_buf = cl.Buffer(self.ctx, mf.READ_WRITE, size=self.max_points * 4 * 4)  # float4 positions
        self.ps_out_buf = cl.Buffer(self.ctx, mf.READ_WRITE, size=self.max_points * 4 * 4)  # float4 output positions
        self.Rs_buf = cl.Buffer(self.ctx, mf.READ_ONLY, size=self.max_points * 4)  # float collision radii
        
        # Neighbor data
        total_neighs = self.max_points * self.max_neighs
        self.neighs_buf = cl.Buffer(self.ctx, mf.READ_WRITE, size=total_neighs * 4)  # int neighbor indices
        self.kngs_buf = cl.Buffer(self.ctx, mf.READ_WRITE, size=total_neighs * 4)    # float stiffnesses
        self.l0s_buf = cl.Buffer(self.ctx, mf.READ_WRITE, size=total_neighs * 4)     # float rest lengths
        
        # Group data for collision detection
        self.max_groups = self.max_points // 32 + 1
        self.granges_buf = cl.Buffer(self.ctx, mf.READ_ONLY, size=self.max_groups * 4 * 4)  # int4 group ranges
        self.inds_buf = cl.Buffer(self.ctx, mf.READ_ONLY, size=self.max_points * 4)  # int indices
        
    def update_positions(self, positions, masses):
        """Update particle positions and masses in GPU memory"""
        n = len(positions)
        assert n <= self.max_points, f"Too many particles: {n} > {self.max_points}"
        
        # Convert to float4 array (x,y,z,mass)
        pos_mass = np.zeros((n, 4), dtype=np.float32)
        pos_mass[:, :2] = positions  # 2D positions
        pos_mass[:, 3] = masses
        
        cl.enqueue_copy(self.queue, self.ps_buf, pos_mass)
    
    def update_collision_radii(self, radii):
        """Update particle collision radii in GPU memory"""
        cl.enqueue_copy(self.queue, self.Rs_buf, radii.astype(np.float32))
    
    def update_group_data(self, granges, indices):
        """Update spatial grouping data for collision detection"""
        cl.enqueue_copy(self.queue, self.granges_buf, granges.astype(np.int32))
        cl.enqueue_copy(self.queue, self.inds_buf, indices.astype(np.int32))
    
    def detect_collisions(self, Kcol, Rp):
        """Run collision detection kernel"""
        ngroups = self.max_groups  # or actual number if known
        
        self.prg.updateCollisonNeighbors(
            self.queue, (ngroups * self.local_size,), (self.local_size,),
            np.int32(ngroups),
            np.int32(self.max_points),
            np.int32(self.max_neighs),
            self.granges_buf,
            self.inds_buf,
            self.ps_buf,
            self.Rs_buf,
            self.neighs_buf,
            self.kngs_buf,
            self.l0s_buf,
            np.float32(Kcol),
            np.float32(Rp)
        )
        self.queue.finish()
    
    def solve_constraints(self, inv_dt2, Rd):
        """Run one iteration of Jacobi solver"""
        self.prg.updateJacobi_neighs(
            self.queue, (self.max_points,), None,
            np.int32(self.max_points),
            np.int32(self.max_neighs),
            self.ps_buf,
            self.ps_out_buf,
            self.neighs_buf,
            self.kngs_buf,
            self.l0s_buf,
            np.float32(inv_dt2),
            np.float32(Rd)
        )
        self.queue.finish()
        
        # Swap buffers
        self.ps_buf, self.ps_out_buf = self.ps_out_buf, self.ps_buf
    
    def get_positions(self):
        """Read back particle positions from GPU"""
        pos_mass = np.zeros((self.max_points, 4), dtype=np.float32)
        cl.enqueue_copy(self.queue, pos_mass, self.ps_buf)
        return pos_mass[:, :2]  # return only 2D positions
    
    def get_neighbors(self):
        """Read back neighbor data from GPU for debugging"""
        neighs = np.zeros((self.max_points, self.max_neighs), dtype=np.int32)
        cl.enqueue_copy(self.queue, neighs, self.neighs_buf)
        return neighs



if __name__ == "__main__":
    # Initialize
    solver = CLConstrains(max_points=1000, max_neighs=32)

    # Update particle data
    solver.update_positions(positions, masses)
    solver.update_collision_radii(radii)
    solver.update_group_data(granges, indices)

    # Run simulation step
    solver.detect_collisions(Kcol=1000.0, Rp=0.1)
    for _ in range(n_iterations):
        solver.solve_constraints(inv_dt2=1.0/dt**2, Rd=0.05)

    # Get results
    final_positions = solver.get_positions()