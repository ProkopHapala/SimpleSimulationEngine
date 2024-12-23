import numpy as np
import pyopencl as cl
import os
import matplotlib.pyplot as plt

from sparse import (
    build_neighbor_list, neigh_stiffness, make_Aii, neighs_to_dense_arrays,
    linsolve_Jacobi_sparse, jacobi_iteration_sparse, linsolve_iterative, color_graph
)
from projective_dynamics import make_pd_rhs, make_pd_Aii0, makeSparseSystem
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
        self.color_group_buf = None
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
        mf = cl.mem_flags
        
        # Pad arrays for OpenCL
        x_padded = np.zeros((len(pos_pred), 4), dtype=np.float32)
        x_padded[:, :3] = pos_pred
        b_padded = np.zeros_like(x_padded)
        b_padded[:, :3] = b
        Aii = Aii.astype(np.float32)
        
        # Create dense arrays from sparse data
        padded_neighs, padded_kngs, n_max = neighs_to_dense_arrays(neighs, kngs, self.n_max if self.n_max is not None else max(len(ng) for ng in neighs))
        
        # Create buffers
        self.x_buf      = cl.Buffer(self.ctx, mf.READ_WRITE | mf.COPY_HOST_PTR, hostbuf=x_padded)
        self.x_out_buf  = cl.Buffer(self.ctx, mf.READ_WRITE, x_padded.nbytes)
        self.b_buf      = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=b_padded)
        self.neighs_buf = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=padded_neighs)
        self.kngs_buf   = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=padded_kngs)
        self.Aii_buf    = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=Aii)
        self.r_buf      = cl.Buffer(self.ctx, mf.READ_WRITE, x_padded.nbytes)

    def _prepare_color_groups(self, neighs):
        """Prepare color groups for Gauss-Seidel"""
        # Color the graph
        colors, color_groups = color_graph(neighs)
        
        # Convert color groups to flat array format
        color_counts = [len(group) for group in color_groups]
        flat_groups = np.concatenate(color_groups).astype(np.int32)
        
        # Create color group buffer
        mf = cl.mem_flags
        self.color_group_buf = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=flat_groups)
        
        return color_counts

    def solve_pd_arrays(self, pos_pred, neighs, kngs, Aii, b, n_max, niter=10, tol=1e-6, nsub=1, errs=None, bPrint=False):
        """Solve projective dynamics using OpenCL with pre-computed arrays"""
        # Create OpenCL buffers from input arrays
        self._create_ocl_buffers(pos_pred, None, neighs, kngs, Aii, b)
        self.n_max = n_max
        self.n_points = len(pos_pred)
        
        # Prepare output arrays
        x_out = np.zeros((self.n_points, 4), dtype=np.float32)
        r = np.zeros_like(x_out)
        
        # Main iteration loop
        itr = 0
        while itr < niter:
            # Run nsub iterations without checking convergence
            events = []
            last_write_buf = None  # Keep track of which buffer has the latest results
            
            for sub_iter in range(nsub):
                # Run kernel
                if sub_iter % 2 == 0:  # Even iterations
                    event = self.prg.jacobi_iteration_sparse(
                        self.queue, (self.n_points,), None,
                        self.x_buf, self.b_buf, self.neighs_buf, self.kngs_buf, self.Aii_buf,
                        self.x_out_buf, self.r_buf,
                        np.int32(self.n_max), np.int32(self.n_points)
                    )
                    last_write_buf = self.x_out_buf
                else:  # Odd iterations
                    event = self.prg.jacobi_iteration_sparse(
                        self.queue, (self.n_points,), None,
                        self.x_out_buf, self.b_buf, self.neighs_buf, self.kngs_buf, self.Aii_buf,
                        self.x_buf, self.r_buf,
                        np.int32(self.n_max), np.int32(self.n_points)
                    )
                    last_write_buf = self.x_buf
                events.append(event)
            
            cl.wait_for_events(events)
            itr += nsub
            
            # Read results back (from the last output buffer used)
            cl.enqueue_copy(self.queue, r, self.r_buf)
            cl.enqueue_copy(self.queue, x_out, last_write_buf)
            
            # Update x_buf for next iteration with the latest results
            if last_write_buf != self.x_buf:
                cl.enqueue_copy(self.queue, self.x_buf, x_out)
            
            # Check convergence
            #err = np.linalg.norm(r[:, :3])  # Use all components for now
            err = np.sum( r*r, axis=0 )
            if errs is not None:
                errs.append( np.sqrt(err) )  # Store error value
            if bPrint:
                print(f"TrussOpenCLSolver::solve_pd() itr: {itr}, err: {err}")
            if np.sqrt(err.sum()) < tol:
                break

        return x_out[:, :3]  # Return only xyz coordinates

    def solve_pd_gauss_seidel(self, pos_pred, neighs, kngs, Aii, b, n_max, niter=10, tol=1e-6, errs=None, bPrint=False):
        """Solve projective dynamics using OpenCL with colored Gauss-Seidel method"""
        # Create OpenCL buffers
        self._create_ocl_buffers(pos_pred, None, neighs, kngs, Aii, b)
        self.n_max = n_max
        self.n_points = len(pos_pred)
        
        # Prepare color groups
        color_counts = self._prepare_color_groups(neighs)
        color_starts = np.cumsum([0] + color_counts[:-1])
        
        # Prepare output arrays
        x_out = np.zeros((self.n_points, 4), dtype=np.float32)
        r = np.zeros_like(x_out)
        
        # Main iteration loop
        for itr in range(niter):
            # Process each color group in sequence
            for color_idx, (start, count) in enumerate(zip(color_starts, color_counts)):
                event = self.prg.gauss_seidel_iteration_colored(
                    self.queue, (count,), None,
                    self.x_buf, self.b_buf, self.neighs_buf, self.kngs_buf, self.Aii_buf,
                    self.color_group_buf, self.r_buf,
                    np.int32(self.n_max), np.int32(self.n_points),
                    np.int32(start), np.int32(count)
                )
                event.wait()
            
            # Check convergence after processing all colors
            cl.enqueue_copy(self.queue, r, self.r_buf)
            cl.enqueue_copy(self.queue, x_out, self.x_buf)
            
            err = np.sum(r*r, axis=0)
            if errs is not None:
                errs.append(np.sqrt(err))
            if bPrint:
                print(f"TrussOpenCLSolver::solve_pd_gauss_seidel() itr: {itr}, err: {err}")
            if np.sqrt(err.sum()) < tol:
                break
        
        return x_out[:, :3]  # Return only xyz coordinates

    def solve_pd_gauss_seidel_truss(self, truss, dt, gravity, fixed_points=None, niter=10, tol=1e-6, errs=None, bPrint=False):
        """Solve projective dynamics for a truss using colored Gauss-Seidel"""
        pos_pred, masses, neighs, kngs, Aii, b = self._init_pd_system(truss, dt, gravity, fixed_points)
        return self.solve_pd_gauss_seidel(pos_pred, neighs, kngs, Aii, b, self.n_max, niter, tol, errs, bPrint)

    def solve_pd(self, truss, dt, gravity, fixed_points=None, niter=10, tol=1e-6, nsub=1, errs=None, bPrint=False):
        """Original solve_pd that initializes arrays internally"""
        pos_pred, masses, neighs, kngs, Aii, b = self._init_pd_system(truss, dt, gravity, fixed_points)
        return self.solve_pd_arrays(pos_pred, neighs, kngs, Aii, b, self.n_max, niter, tol, nsub, errs, bPrint)

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
    truss.build_grid_2d(nx=2, ny=2, m=1.0, m_end=1000.0, l=1.0, k=10000.0, k_diag=1000.0)
    
    # Get initial data
    bonds, points, masses, ks, fixed, l0s, neighbs = truss.get_pd_quantities()
    neighbs = build_neighbor_list(bonds, len(points))
    velocity = points * 0 + gravity[None,:] * dt
    pos_pred = points + velocity * dt
    if fixed_points is not None:
        pos_pred[fixed_points] = points[fixed_points]
    
    # Reference solution
    print("\nReference solution:")
    neighs, kngs, Aii, b, pos_pred, velocity, n_max = makeSparseSystem( dt, bonds, points, masses, ks, fixed, l0s, neighbs )

    # Run one iteration of CPU solver
    print("\nCPU Iteration:")
    x_ref, r_ref = jacobi_iteration_sparse(pos_pred[:,0], b[:,0], neighs, kngs, Aii)  
    
    # OpenCL solution
    print("\nGPU Iteration:")
    solver = TrussOpenCLSolver()

    solver._create_ocl_buffers(pos_pred, masses, neighs, kngs, Aii, b)
    
    solver.n_points = len(pos_pred)
    solver.n_max = n_max

    # Run one iteration
    x_out = np.zeros((solver.n_points, 4), dtype=np.float32)
    r = np.zeros_like(x_out)

    solver.prg.jacobi_iteration_sparse(
        solver.queue, (solver.n_points,), None,
        solver.x_buf, solver.b_buf, solver.neighs_buf, solver.kngs_buf, solver.Aii_buf,
        solver.x_out_buf, solver.r_buf,
        np.int32(solver.n_max), np.int32(solver.n_points)
    )
    
    # Read results back
    cl.enqueue_copy(solver.queue, r, solver.r_buf)
    cl.enqueue_copy(solver.queue, x_out, solver.x_out_buf)
    x_ocl = x_out[:,0]
    
    # Compare results
    print("\nComparison:")
    diff = np.abs(x_ref - x_ocl)
    max_diff = np.max(diff)
    print(f"Max difference between reference and OpenCL: {max_diff}")

    if max_diff > 1e-6:
        print("\nReference solution:")
        print(x_ref)
        print("\nOpenCL solution:")
        print(x_ocl)
        print("\nError:")
        print(diff)
        print("\nMax difference between reference and OpenCL: ", max_diff)
        print("EROOR in test_against_reference() => exit()")
        exit()
    
    return x_ref, x_ocl

if __name__ == "__main__":
    from projective_dynamics import makeSparseSystem
    from sparse import linsolve_iterative
    import matplotlib.pyplot as plt

    # First run reference comparison test
    test_against_reference()


    #exit()
    
    print("\nRunning full solver comparison...")
    # Create test truss (grid)
    truss = Truss()
    truss.build_grid_2d(nx=5, ny=5, m=1.0, m_end=1000.0, l=1.0, k=10000.0, k_diag=1000.0)
    
    # Get system arrays
    dt = 1.0
    gravity = np.array([0., -9.81, 0.0])
    fixed_points = [0]
    bonds, points, masses, ks, fixed, l0s, neighbs = truss.get_pd_quantities()
    neighbs = build_neighbor_list(bonds, len(points))
    neighs, kngs, Aii, b, pos_pred, velocity, n_max = makeSparseSystem(dt, bonds, points, masses, ks, fixed, l0s, neighbs)
    
    # Run both solvers
    errs_cpu_x = []
    errs_cpu_y = []
    errs_gpu = []
    niter = 100

    bPrint = False
    
    # CPU solver (just x component)
    x_cpu = linsolve_iterative( lambda A, b, x: jacobi_iteration_sparse(x, b, neighs, kngs, Aii),  b[:,0] , None, pos_pred[:,0], niter=niter, tol=1e-6, bPrint=bPrint, errs=errs_cpu_x  )
    y_cpu = linsolve_iterative( lambda A, b, x: jacobi_iteration_sparse(x, b, neighs, kngs, Aii),  b[:,1] , None, pos_pred[:,1], niter=niter, tol=1e-6, bPrint=bPrint, errs=errs_cpu_y  )
    
    # GPU solver (all components)
    solver = TrussOpenCLSolver()
    x_gpu = solver.solve_pd_arrays( pos_pred, neighs, kngs, Aii, b, n_max,  niter=niter, tol=1e-6, errs=errs_gpu, bPrint=bPrint  )
    errs_gpu = np.array(errs_gpu)
    
    # Plot convergence comparison
    plt.figure(figsize=(10, 6))
    plt.semilogy(errs_cpu_x, 'r:', label='CPU x')
    plt.semilogy(errs_cpu_y, 'g:', label='CPU y')
    plt.semilogy(errs_gpu[:,0], 'r-', lw=0.5, label='GPU x')
    plt.semilogy(errs_gpu[:,1], 'g-', lw=0.5, label='GPU y')
    #plt.semilogy(errs_gpu[:,2], 'b-', label='GPU z')
    plt.grid(True)
    plt.xlabel('Iteration')
    plt.ylabel('Error (log scale)')
    plt.title('Convergence Comparison: CPU vs GPU Solver')
    plt.legend()
    plt.show()