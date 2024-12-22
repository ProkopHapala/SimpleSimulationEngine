import numpy as np
import pyopencl as cl
import os

from sparse import (
    build_neighbor_list, neigh_stiffness, make_Aii, neighs_to_dense_arrays,
    linsolve_Jacobi_sparse
)
from projective_dynamics import make_pd_rhs, make_pd_Aii0
from truss import Truss

class TrussOpenCLSolver:
    def __init__(self):
        # Initialize OpenCL
        platforms = cl.get_platforms()
        self.ctx = cl.Context(
            dev_type=cl.device_type.ALL,
            properties=[(cl.context_properties.PLATFORM, platforms[0])]
        )
        self.queue = cl.CommandQueue(self.ctx)
        
        # Load and build the OpenCL kernel
        kernel_path = os.path.join(os.path.dirname(__file__), 'truss.cl')
        with open(kernel_path, 'r') as f:
            kernel_src = f.read()
        self.prg = cl.Program(self.ctx, kernel_src).build()
        
        # Initialize buffers to None
        self.x_buf = None
        self.x_out_buf = None
        self.b_buf = None
        self.neighs_buf = None
        self.kngs_buf = None
        self.Aii_buf = None
        self.r_buf = None
        self.n_points = None
        self.n_max = None

    def _init_pd_system(self, truss, dt, gravity, fixed_points):
        """Initialize the projective dynamics system"""
        # Get truss data
        bonds, points, masses, ks, fixed, l0s, neighbs = truss.get_pd_quantities()
        neighbs = build_neighbor_list(bonds, len(points))
        
        # Compute initial state
        velocity = points * 0 + np.array(gravity)[None, :] * dt
        pos_pred = points + velocity * dt
        if fixed_points is not None:
            pos_pred[fixed_points] = points[fixed_points]
        
        # Prepare sparse data
        neighs, kngs, n_max = neigh_stiffness(neighbs, bonds, ks)
        Aii0 = make_pd_Aii0(masses, dt)
        Aii  = make_Aii(neighs, kngs, Aii0)
        b    = make_pd_rhs(neighbs, bonds, masses, dt, ks, points, l0s, pos_pred)
        
        self.n_max = n_max
        return pos_pred, masses, neighs, kngs, Aii, b

    def _create_ocl_buffers(self, pos_pred, masses, neighs, kngs, Aii, b):
        """Create and initialize OpenCL buffers"""
        # Convert and pad data for OpenCL
        x = pos_pred.astype(np.float32)
        b = b.astype(np.float32)
        Aii = Aii.astype(np.float32)
        
        # Use neighs_to_dense_arrays to get padded arrays
        padded_neighs, padded_kngs, _ = neighs_to_dense_arrays(neighs, kngs, self.n_max)
        
        # Create float4 arrays for positions and forces
        self.n_points = len(pos_pred)
        x_padded = np.zeros((self.n_points, 4), dtype=np.float32)
        b_padded = np.zeros((self.n_points, 4), dtype=np.float32)
        x_padded[:, :3] = x
        x_padded[:,  3] = masses  # Store mass in w component
        b_padded[:, :3] = b
        
        # Create OpenCL buffers
        mf = cl.mem_flags
        self.x_buf      = cl.Buffer(self.ctx, mf.READ_WRITE | mf.COPY_HOST_PTR, hostbuf=x_padded)
        self.x_out_buf  = cl.Buffer(self.ctx, mf.READ_WRITE, x_padded.nbytes)
        self.b_buf      = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=b_padded)
        self.neighs_buf = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=padded_neighs)
        self.kngs_buf   = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=padded_kngs)
        self.Aii_buf    = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=Aii)
        self.r_buf      = cl.Buffer(self.ctx, mf.READ_WRITE, x_padded.nbytes)

    def solve_pd(self, truss, dt, gravity, fixed_points=None, niter=10, tol=1e-6, nsub=1):
        # Initialize system and create buffers
        pos_pred, masses, neighs, kngs, Aii, b = self._init_pd_system(truss, dt, gravity, fixed_points)
        self._create_ocl_buffers(pos_pred, masses, neighs, kngs, Aii, b)
        
        # Solve using OpenCL
        x_out = np.zeros((self.n_points, 4), dtype=np.float32)
        r = np.zeros_like(x_out)

        itr = 0
        while itr < niter:
            # python functions are costly therefore we call multiple times the kernell before we read the results
            for itr_sub in range(nsub):
                self.prg.jacobi_iteration_sparse(
                    self.queue, (self.n_points,), None,
                    self.x_buf, self.b_buf, self.neighs_buf, self.kngs_buf, self.Aii_buf,
                    self.x_out_buf, self.r_buf,
                    np.int32(self.n_max), np.int32(self.n_points)
                )
            itr += nsub
            
            # Read results back
            cl.enqueue_copy(self.queue, r, self.r_buf)
            cl.enqueue_copy(self.queue, x_out, self.x_out_buf)
            
            # Check convergence
            err = np.linalg.norm(r[:, :3])
            if err < tol:
                break
            
            # Update x for next iteration
            cl.enqueue_copy(self.queue, self.x_buf, x_out)
        
        return x_out[:, :3]  # Return only xyz coordinates

    def solve_pd_event(self, truss, dt, gravity, fixed_points=None, niter=10, tol=1e-6, nsub=1):
        # Initialize system and create buffers
        pos_pred, masses, neighs, kngs, Aii, b = self._init_pd_system(truss, dt, gravity, fixed_points)
        self._create_ocl_buffers(pos_pred, masses, neighs, kngs, Aii, b)
        
        # Solve using OpenCL with event-based execution
        x_out = np.zeros((self.n_points, 4), dtype=np.float32)
        r = np.zeros_like(x_out)

        itr = 0
        while itr < niter:
            # Queue multiple kernel executions with events
            events = []
            prev_event = None
            for _ in range(nsub):
                event = self.prg.jacobi_iteration_sparse(
                    self.queue, (self.n_points,), None,
                    self.x_buf, self.b_buf, self.neighs_buf, self.kngs_buf, self.Aii_buf,
                    self.x_out_buf, self.r_buf,
                    np.int32(self.n_max), np.int32(self.n_points),
                    wait_for=[prev_event] if prev_event else None
                )
                events.append(event)
                # Swap buffers for next iteration
                self.x_buf, self.x_out_buf = self.x_out_buf, self.x_buf
                prev_event = event
            
            # Wait for all kernels to complete
            cl.wait_for_events(events)
            itr += nsub
            
            # Read results back (from the last output buffer used)
            cl.enqueue_copy(self.queue, r, self.r_buf)
            cl.enqueue_copy(self.queue, x_out, self.x_buf if nsub % 2 == 0 else self.x_out_buf)
            
            # Check convergence
            err = np.linalg.norm(r[:, :3])
            if err < tol:
                break

        return x_out[:, :3]  # Return only xyz coordinates

def test_against_reference():
    print("\nTesting against reference implementation:")
    dt = 1.0
    gravity = np.array([0., -9.81, 0.0])
    fixed_points = [0]

    # Create test truss (grid)
    truss = Truss()
    truss.build_grid_2d(nx=5, ny=5, m=1.0, m_end=1000.0, l=1.0, k=10000.0, k_diag=1000.0)
    
    # Get initial data
    bonds, points, masses, ks, fixed, l0s, neighbs = truss.get_pd_quantities()
    neighbs = build_neighbor_list(bonds, len(points))
    velocity = points * 0 + gravity[None,:] * dt
    pos_pred = points + velocity * dt
    if fixed_points is not None:
        pos_pred[fixed_points] = points[fixed_points]
    
    # Reference solution
    print("\nReference solution:")
    neighs, kngs, n_max = neigh_stiffness(neighbs, bonds, ks)
    Aii0 = make_pd_Aii0(masses, dt)
    b = make_pd_rhs(neighbs, bonds, masses, dt, ks, points, l0s, pos_pred)
    
    err_cpu = []
    def callback(itr, x, r):
        err = np.linalg.norm(r)
        err_cpu.append(err)
        print(f"Reference itr: {itr}, err: {err}")
    
    x_ref = linsolve_Jacobi_sparse(b, neighs, kngs, x0=pos_pred, Aii0=Aii0, niter=100, tol=1e-6, callback=callback)
    
    # OpenCL solution
    print("\nOpenCL solution:")
    solver = TrussOpenCLSolver()
    
    # Modify solve_pd to track errors
    err_gpu = []
    x_ocl = np.zeros_like(pos_pred)
    r = np.zeros_like(pos_pred)
    
    # Initialize system and create buffers
    pos_pred, masses, neighs, kngs, Aii, b = solver._init_pd_system(truss, dt, gravity, fixed_points)
    solver._create_ocl_buffers(pos_pred, masses, neighs, kngs, Aii, b)
    
    # Solve using OpenCL
    x_out = np.zeros((solver.n_points, 4), dtype=np.float32)
    r = np.zeros_like(x_out)

    for itr in range(100):  # Same number of iterations as CPU
        solver.prg.jacobi_iteration_sparse(
            solver.queue, (solver.n_points,), None,
            solver.x_buf, solver.b_buf, solver.neighs_buf, solver.kngs_buf, solver.Aii_buf,
            solver.x_out_buf, solver.r_buf,
            np.int32(solver.n_max), np.int32(solver.n_points)
        )
        
        # Read results back
        cl.enqueue_copy(solver.queue, r, solver.r_buf)
        cl.enqueue_copy(solver.queue, x_out, solver.x_out_buf)
        
        # Calculate error
        err = np.linalg.norm(r[:, :3])
        err_gpu.append(err)
        print(f"OpenCL itr: {itr}, err: {err}")
        
        if err < 1e-6:
            break
            
        # Update x for next iteration
        cl.enqueue_copy(solver.queue, solver.x_buf, x_out)
    
    x_ocl = x_out[:, :3]
    
    # Compare results
    diff = np.abs(x_ref - x_ocl)
    max_diff = np.max(diff)
    print(f"\nMax difference between reference and OpenCL: {max_diff}")
    
    # Plot error histories
    import matplotlib.pyplot as plt
    plt.figure(figsize=(10, 6))
    plt.plot(err_cpu, label="Jacobi CPU")
    plt.plot(err_gpu, label="Jacobi GPU")
    plt.yscale('log')
    plt.xlabel('Iteration')
    plt.ylabel('Error (log scale)')
    plt.title('Convergence Comparison: CPU vs GPU Implementation')
    plt.grid(True)
    plt.legend()
    plt.show()
    
    return x_ref, x_ocl, err_cpu, err_gpu

if __name__ == "__main__":
    # Run reference comparison test only
    test_against_reference()
