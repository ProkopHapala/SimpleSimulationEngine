import numpy as np
import pyopencl as cl
import os
import matplotlib.pyplot as plt

from sparse import (
    build_neighbor_list, neigh_stiffness, make_Aii, neighs_to_dense_arrays,
    linsolve_Jacobi_sparse, jacobi_iteration_sparse, linsolve_iterative, color_graph,
    gauss_seidel_iteration_colored
)
from projective_dynamics import make_pd_rhs, make_pd_Aii0, makeSparseSystem
from truss import Truss


VBD_WORKGROUP_SIZE = 32
VBD_NEIGHBOR_LIMIT = 128


def build_rest_length_dense(neigh_lists, bonds, l0s, n_max):
    """Dense rest-length table matching the padded neighbor layout."""
    assert len(bonds) == len(l0s), "bonds and l0s must have identical length"
    pair_rest = {}
    for (ia, ib), rest in zip(bonds, l0s):
        pair_rest[(int(ia), int(ib))] = float(rest)
        pair_rest[(int(ib), int(ia))] = float(rest)
    rest_dense = np.zeros((len(neigh_lists), n_max), dtype=np.float32)
    for vid, neigh in enumerate(neigh_lists):
        for jj, nj in enumerate(neigh):
            rest_dense[vid, jj] = pair_rest[(vid, nj)]
    return rest_dense


def build_vbd_chunks(padded_neighs, color_groups, chunk_size=VBD_WORKGROUP_SIZE, neighbor_limit=VBD_NEIGHBOR_LIMIT):
    """Split each color into workgroup batches respecting vertex and neighbor limits."""
    _, n_max = padded_neighs.shape
    chunks_per_color = []
    for group in color_groups:
        color_vertices = list(group)
        cursor = 0
        color_chunks = []
        while cursor < len(color_vertices):
            vertices = []
            neighbors = []
            neighbor_map = {}
            while cursor < len(color_vertices) and len(vertices) < chunk_size:
                vid = int(color_vertices[cursor])
                row = padded_neighs[vid]
                neigh_list = [int(n) for n in row if n >= 0]
                new_neighbors = [n for n in neigh_list if n not in neighbor_map]
                if neighbors and len(neighbors) + len(new_neighbors) > neighbor_limit:
                    break
                if not neighbors and len(new_neighbors) > neighbor_limit:
                    raise ValueError(f"Vertex {vid} uses {len(new_neighbors)} neighbors > limit {neighbor_limit}")
                vertices.append(vid)
                for n in neigh_list:
                    if n not in neighbor_map:
                        neighbor_map[n] = len(neighbors)
                        neighbors.append(n)
                cursor += 1
            if not vertices:
                raise RuntimeError("Unable to pack VBD workgroup; adjust neighbor limit or chunk size")
            slot_table = np.full((chunk_size, n_max), -1, dtype=np.int32)
            for local_idx, vid in enumerate(vertices):
                row = padded_neighs[vid]
                for jj in range(n_max):
                    nid = row[jj]
                    if nid < 0:
                        break
                    slot_table[local_idx, jj] = neighbor_map[nid]
            color_chunks.append({
                "vertices": np.array(vertices, dtype=np.int32),
                "neighbors": np.array(neighbors, dtype=np.int32),
                "slot_table": slot_table,
            })
        chunks_per_color.append(color_chunks)
    return chunks_per_color


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
        self.y_buf = None
        self.l0_buf = None
        self.neighs_dense = None
        self.kngs_dense = None
        self.rest_lengths_dense = None
        self.vbd_chunks = None
        self.color_groups_host = None
        self.last_bonds = None
        self.last_l0s = None
        self.last_neigh_lists = None

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
        self.last_bonds = bonds
        self.last_l0s = l0s
        self.last_neigh_lists = neighs
        return pos_pred, masses, neighs, kngs, Aii, b

    def _create_ocl_buffers(self, pos_pred, masses, neighs, kngs, Aii, b, rest_lengths=None, y=None):
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
        self.neighs_dense = padded_neighs
        self.kngs_dense = padded_kngs
        
        # Create buffers
        self.x_buf      = cl.Buffer(self.ctx, mf.READ_WRITE | mf.COPY_HOST_PTR, hostbuf=x_padded)
        self.x_out_buf  = cl.Buffer(self.ctx, mf.READ_WRITE, x_padded.nbytes)
        self.b_buf      = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=b_padded)
        self.neighs_buf = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=padded_neighs)
        self.kngs_buf   = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=padded_kngs)
        self.Aii_buf    = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=Aii)
        self.r_buf      = cl.Buffer(self.ctx, mf.READ_WRITE, x_padded.nbytes)

        if rest_lengths is not None:
            assert rest_lengths.shape == padded_neighs.shape, "rest_lengths must match neighbor layout"
            self.rest_lengths_dense = rest_lengths.astype(np.float32, copy=False)
            self.l0_buf = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=self.rest_lengths_dense)
        else:
            self.rest_lengths_dense = None
            self.l0_buf = None

        if y is not None:
            y_padded = np.zeros_like(x_padded)
            y_padded[:, :3] = y
            self.y_buf = cl.Buffer(self.ctx, mf.READ_ONLY | mf.COPY_HOST_PTR, hostbuf=y_padded)
        else:
            self.y_buf = None

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
        
        self.color_groups_host = color_groups
        return color_counts

    def build_vbd_work_chunks(self, chunk_size=VBD_WORKGROUP_SIZE, neighbor_limit=VBD_NEIGHBOR_LIMIT):
        """Pack per-color vertex chunks for the VBD kernel."""
        if self.neighs_dense is None:
            raise ValueError("Dense neighbor table not initialized; call _create_ocl_buffers first")
        if self.color_groups_host is None:
            raise ValueError("Color groups not prepared; call _prepare_color_groups first")
        self.vbd_chunks = build_vbd_chunks(self.neighs_dense, self.color_groups_host, chunk_size, neighbor_limit)
        return self.vbd_chunks

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

    def solve_pd_gauss_seidel(self, pos_pred, neighs, kngs, Aii, b, n_max, niter=10, tol=1e-6, errs=None, bPrint=False, *, masses=None, rest_lengths=None, y=None):
        """Solve projective dynamics using OpenCL with colored Gauss-Seidel method"""
        # Create OpenCL buffers
        self._create_ocl_buffers(pos_pred, masses, neighs, kngs, Aii, b, rest_lengths=rest_lengths, y=y)
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
        return self.solve_pd_gauss_seidel(pos_pred, neighs, kngs, Aii, b, self.n_max, niter, tol, errs, bPrint, masses=masses)

    def solve_pd(self, truss, dt, gravity, fixed_points=None, niter=10, tol=1e-6, nsub=1, errs=None, bPrint=False):
        """Original solve_pd that initializes arrays internally"""
        pos_pred, masses, neighs, kngs, Aii, b = self._init_pd_system(truss, dt, gravity, fixed_points)
        return self.solve_pd_arrays(pos_pred, neighs, kngs, Aii, b, self.n_max, niter, tol, nsub, errs, bPrint, masses=masses)

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

    def solve_vbd(self, truss, dt, gravity, fixed_points=None, niter=1, det_eps=1e-6, errs=None, bPrint=False):
        """Run Vertex Block Descent iterations on the GPU."""
        if dt <= 0:
            raise ValueError("dt must be positive")
        pos_pred, masses, neighs, kngs, Aii, b = self._init_pd_system(truss, dt, gravity, fixed_points)
        inv_h2 = np.float32(1.0 / (dt * dt))

        rest_dense = build_rest_length_dense(self.last_neigh_lists, self.last_bonds, self.last_l0s, self.n_max)
        y_pred = pos_pred.astype(np.float32, copy=True)

        self._create_ocl_buffers(pos_pred, masses, neighs, kngs, Aii, b, rest_lengths=rest_dense, y=y_pred)
        self.n_points = len(pos_pred)

        color_counts = self._prepare_color_groups(neighs)
        self.build_vbd_work_chunks()

        if self.l0_buf is None or self.y_buf is None:
            raise RuntimeError("VBD buffers not initialized (l0/y missing)")

        mf = cl.mem_flags

        chunk_vertices_host = np.full(VBD_WORKGROUP_SIZE, -1, dtype=np.int32)
        chunk_neighbors_host = np.full(VBD_NEIGHBOR_LIMIT, -1, dtype=np.int32)
        chunk_slot_host = np.full(VBD_WORKGROUP_SIZE * self.n_max, -1, dtype=np.int32)

        chunk_vertices_buf = cl.Buffer(self.ctx, mf.READ_ONLY, chunk_vertices_host.nbytes)
        chunk_neighbors_buf = cl.Buffer(self.ctx, mf.READ_ONLY, chunk_neighbors_host.nbytes)
        chunk_slot_buf = cl.Buffer(self.ctx, mf.READ_ONLY, chunk_slot_host.nbytes)

        for itr in range(niter):
            if bPrint:
                print(f"TrussOpenCLSolver::solve_vbd() itr {itr}")
            for color_chunks in self.vbd_chunks:
                for chunk in color_chunks:
                    vertices = chunk["vertices"]
                    neighbors = chunk["neighbors"]
                    slot_table = chunk["slot_table"]

                    vcount = int(len(vertices))
                    ncount = int(len(neighbors))

                    chunk_vertices_host.fill(-1)
                    if vcount > 0:
                        chunk_vertices_host[:vcount] = vertices

                    chunk_neighbors_host.fill(-1)
                    if ncount > 0:
                        chunk_neighbors_host[:ncount] = neighbors

                    np.copyto(chunk_slot_host, -1)
                    slot_flat = slot_table.reshape(-1).astype(np.int32, copy=False)
                    chunk_slot_host[:slot_flat.size] = slot_flat

                    cl.enqueue_copy(self.queue, chunk_vertices_buf, chunk_vertices_host)
                    cl.enqueue_copy(self.queue, chunk_neighbors_buf, chunk_neighbors_host)
                    cl.enqueue_copy(self.queue, chunk_slot_buf, chunk_slot_host)

                    event = self.prg.vbd_vertex_chunk(
                        self.queue, (VBD_WORKGROUP_SIZE,), (VBD_WORKGROUP_SIZE,),
                        self.x_buf, self.y_buf,
                        chunk_vertices_buf, chunk_neighbors_buf, chunk_slot_buf,
                        self.neighs_buf, self.kngs_buf, self.l0_buf,
                        np.int32(vcount), np.int32(ncount),
                        np.int32(self.n_max), inv_h2, np.float32(det_eps)
                    )
                    event.wait()

        x_out = np.zeros((self.n_points, 4), dtype=np.float32)
        cl.enqueue_copy(self.queue, x_out, self.x_buf)
        v_out = None
        if self.y_buf is not None:
            y_out = np.zeros_like(x_out)
            cl.enqueue_copy(self.queue, y_out, self.y_buf)
            v_out = (x_out[:, :3] - y_out[:, :3]) / np.float32(dt)
        if errs is not None:
            errs.append(0.0)
        if v_out is not None:
            return x_out[:, :3], v_out
        return x_out[:, :3], None

def setup_test_problem(nx=2, ny=2):
    """Common setup for test problems"""
    dt = 1.0
    gravity = np.array([0., -9.81, 0.0])
    fixed_points = [0]

    # Create test truss (grid)
    truss = Truss()
    truss.build_grid_2d(nx=nx, ny=ny, m=1.0, m_end=1000.0, l=1.0, k=10000.0, k_diag=1000.0)
    
    # Get truss data
    bonds, points, masses, ks, fixed, l0s, neighbs = truss.get_pd_quantities()
    neighbs = build_neighbor_list(bonds, len(points))
    
    # Compute initial state
    velocity = points * 0 + gravity[None, :] * dt
    pos_pred = points + velocity * dt
    if fixed_points is not None:
        pos_pred[fixed_points] = points[fixed_points]
    
    # Prepare sparse data
    neighs, kngs, n_max = neigh_stiffness(neighbs, bonds, ks)
    Aii0 = make_pd_Aii0(masses, dt)
    Aii  = make_Aii(neighs, kngs, Aii0)
    b    = make_pd_rhs(neighbs, bonds, masses, dt, ks, points, l0s, pos_pred)
    
    return truss, pos_pred, neighs, kngs, Aii, b, n_max

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

def test_Jacobi_vs_python():
    print("\nRunning Jacobi solver comparison...")
    
    # Setup test problem
    truss, pos_pred, neighs, kngs, Aii, b, n_max = setup_test_problem(nx=2, ny=2)
    
    # Run CPU solver
    errs_cpu = []
    x_cpu = linsolve_iterative(
        lambda A, b, x: jacobi_iteration_sparse(x, b, neighs, kngs, Aii),
        b, None, pos_pred, niter=50, tol=1e-10, errs=errs_cpu, bPrint=True
    )
    
    # Run GPU solver
    errs_gpu = []
    solver = TrussOpenCLSolver()
    x_gpu = solver.solve_pd_arrays(pos_pred, neighs, kngs, Aii, b, n_max, niter=50, tol=1e-10, nsub=1, errs=errs_gpu, bPrint=True)
    
    # Plot convergence
    plt.figure(figsize=(10, 6))
    plt.semilogy(errs_cpu, 'b-', label='CPU Jacobi')
    plt.semilogy(errs_gpu, 'r--', label='GPU Jacobi')
    plt.grid(True)
    plt.xlabel('Iteration')
    plt.ylabel('Error (log scale)')
    plt.title('Convergence Comparison: CPU vs GPU Jacobi Solver')
    plt.legend()
    plt.show()

def test_GaussSeidel_vs_python():
    print("\nRunning Gauss-Seidel solver comparison...")
    
    # Setup test problem
    truss, pos_pred, neighs, kngs, Aii, b, n_max = setup_test_problem(nx=2, ny=2)
    
    # Color the graph for CPU solver
    colors, color_groups = color_graph(neighs)
    
    # Run CPU solver
    errs_cpu = []
    x_cpu = linsolve_iterative(
        lambda A, b, x: gauss_seidel_iteration_colored(x, b, neighs, kngs, Aii, color_groups),
        b, None, pos_pred, niter=50, tol=1e-10, errs=errs_cpu, bPrint=True
    )
    
    # Run GPU solver
    errs_gpu = []
    solver = TrussOpenCLSolver()
    x_gpu = solver.solve_pd_gauss_seidel(pos_pred, neighs, kngs, Aii, b, n_max, niter=50, tol=1e-10, errs=errs_gpu, bPrint=True)
    
    # Plot convergence
    plt.figure(figsize=(10, 6))
    plt.semilogy(errs_cpu, 'b-', label='CPU Gauss-Seidel')
    plt.semilogy(errs_gpu, 'r--', label='GPU Gauss-Seidel')
    plt.grid(True)
    plt.xlabel('Iteration')
    plt.ylabel('Error (log scale)')
    plt.title('Convergence Comparison: CPU vs GPU Gauss-Seidel Solver')
    plt.legend()
    plt.show()

if __name__ == "__main__":
    from projective_dynamics import makeSparseSystem
    from sparse import (
        linsolve_iterative, jacobi_iteration_sparse,
        gauss_seidel_iteration_colored, color_graph
    )
    import matplotlib.pyplot as plt

    # Run tests
    #test_against_reference()
    #test_Jacobi_vs_python()
    test_GaussSeidel_vs_python()
