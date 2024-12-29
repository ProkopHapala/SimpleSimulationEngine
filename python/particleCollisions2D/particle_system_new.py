import numpy as np
import matplotlib.pyplot as plt

def rhs_ij( pi, pj, k, l0 ):
    '''
    b_i = \sum_j (K_{ij} p_{ij})
    where p_{ij} is ideal position of particle i which satisfy the constraint between i and j
    p_{ij} = p_j + (p_i-p_j)/|p_i-p_j| * l0
    '''
    d = pi - pj
    l = np.linalg.norm(d)
    return d * (k * l0/l)

def make_PD_RHS( ps, neighs, kngs, l0s, masses=None, dt=1.0 ):
    '''
    b_i = \sum_j (K_{ij} p_{ij})  + M/dt^2 p'_i
    '''
    n = len(ps)
    b = np.zeros((n,2))
    inv_dt2 = 1.0/(dt*dt);
    for i in range(n):
        pi = ps[i]
        bi = np.zeros(2)
        ngsi  = neighs[i]
        ksi   = kngs  [i]
        l0i   = l0s   [i]
        ni    = len(ngsi) 
        for jj in range(ni):
            k  = ksi [jj]
            l0 = l0i [jj] 
            j  = ngsi[jj]
            pj = ps  [j ]
            bi += rhs_ij( pi, pj, k, l0 )
        if masses is not None:
            bi += masses[i] * inv_dt2    # for the moment we neglect the inertial term
        b[i] = bi
    return b

def update_PD_matrix( apos, neighs, kngs, dt=1.0 ):
    """
    Update the diagonal matrix for Projective Dynamics
    """
    n = len(apos)
    Aii = np.zeros(n)
    for ip in range(n):
        # Inertial term: m_i/dt^2
        Aii[ip] = 1.0 / (dt**2)  # Using unit mass
        # Sum of stiffness for all neighbors
        for j, k in zip(neighs[ip], kngs[ip]):
            Aii[ip] += k
    return Aii

def jacobi_iteration_sparse(x, b, neighs, kngs, Aii ):
    """One iteration of Jacobi method using sparse operations"""
    n = len(x)
    x_out = np.zeros_like(x)
    r     = np.zeros_like(x)
    for i in range(n):
        #print("CPU i: %i Aii[i]: %f b[i]: %f " %(i, Aii[i], b[i]) );
        sum_j  = 0  # RHS term
        ngsi = neighs[i]
        ksi  = kngs[i] 
        ni = len(ngsi)
        for jj in range(ni):
            j      = ngsi[jj]
            k      = ksi[jj]
            sum_j -= k * x[j]   # Off-diagonal contribution
        x_out[i] =  (b[i] - sum_j) / Aii[i]   # solution x_new = (b - sum_(j!=i){ Aij * x[j] } ) / Aii
        y_i = Aii[i]*x[i] + sum_j              # This is Ax_i
        Api_ = b[i] - sum_j                    # Api_ = sum_j K_ij ( d_ij + p_j ) 
        r[i] = b[i] - y_i                      # Residual r = b - Ax
        print("i: %i r_i: %f pi_: %f x_i: %f " %(i, r[i], Api_/Aii[i], x[i] ) )
    return x_out, r

def jacobi_iteration_sparse_debug(ps, b, neighs, kngs, Aii, bPlot=True, l0s=None):
    """Debug version of Jacobi iteration that plots p'_ij points"""
    n = len(ps)
    ps_out = np.zeros_like(ps)
    r = np.zeros_like(ps)
    
    if bPlot: 
        plt.scatter(ps[:,0], ps[:,1], c='blue', label='Current positions')
    
    # For each particle
    for i in range(n):
        pi    = ps[i]
        bi    = b [i]    #   bi = \sum_j (K_{ij} d_{ij})
        ngsi  = neighs[i]
        ksi   = kngs[i]
        ni    = len(ngsi)
        sum_j = np.zeros(2)  

        sum_pij = np.zeros(2)  
        for jj in range(ni):
            j  = ngsi[jj]
            k  = ksi [jj]
            pj = ps  [j]
            sum_j += k * pj
            if bPlot:
                l0 = l0s[i][jj]
                d_ij = rhs_ij( pi, pj, 1.0, l0 )
                p_ij = pj + d_ij

                sum_pij += p_ij * k
                plt.plot(p_ij[0], p_ij[1], '.r', alpha=0.5)                     # Plot p'_ij point
                plt.plot([pj[0], p_ij[0]], [pj[1], p_ij[1]], 'k--', lw=0.5, alpha=0.3 ) # Draw line from current position to p'_ij
                
        # Note 
        # d_{ij}  = (p_i - p_j)/|p_i-p_j| * l0
        # bi      = \sum_j K_{ij} d_{ij}
        # sum_j   = \sum_j K_{ij} p_j
        # p'_{ij} = p_j + d_{ij}
        # A_{ii}  = \sum_j K_{ij}
        # p^{new}_i = (bi + sum_j) / A_{ii} 
        #           = \sum_j K_{ij} (d_{ij} + p_j) / A_{ii} 
        #           = K_{ij} p'_{ij} / A_{ii} 
        #           = {sum_j K_{ij} p'_{ij} } / { sum_j K_{ij} }
        pi_new = (bi + sum_j) / Aii[i]

        pi_ij = sum_pij / Aii[i]

        if bPlot:
            plt.plot(pi_new[0], pi_new[1], '+k', alpha=0.5)                     # Plot p'_ij point
            
            plt.plot([pi[0], pi_new[0]], [pi[1], pi_new[1]], 'k--', lw=0.5, alpha=1.0 ) # Draw line from current position to p'_ij

            plt.plot(pi_ij[0], pi_ij[1], '+r', alpha=0.5)  
            plt.plot([pi[0], pi_ij[0]], [pi[1], pi_ij[1]], 'r--', lw=0.5, alpha=1.0 ) 
        
        ps_out[i] = pi_new

        # Calculate residual
        y_i = Aii[i]*ps[i] + sum_j  # This is Ax_i
        r[i] = b[i] - y_i  # Residual r = b - Ax
        #print(f"i: {i} r_i: {np.linalg.norm(r[i]):.6f} pi_: {np.linalg.norm(b[i]/Aii[i]):.6f}")
    
    # Plot bonds between particles
    # for i in range(n):
    #     for j in neighs[i]:
    #         plt.plot([ps[i,0], ps[j,0]], [ps[i,1], ps[j,1]], 'k-', lw=0.5, alpha=1.0 )
    
    plt.axis('equal')
    plt.grid(True)
    plt.legend()
    plt.title('Debug visualization of Jacobi iteration')
    #plt.show()
    
    return ps_out, r

class Particle:
    def __init__(self, pos, vel=None, mass=1.0):
        self.pos = np.array(pos, dtype=np.float64)    # position vector [x,y]
        self.vel = np.zeros(2) if vel is None else np.array(vel, dtype=np.float64)   # velocity vector [vx,vy]
        self.m = mass                                 # mass
        self.f = np.zeros(2)                         # force accumulator
        self.neighbors = []                          # list of neighbor particles for collision detection
        
class Bond:
    """A bond between two particles"""
    def __init__(self, i, j, k=1.0, r0=1.0):
        self.i = i   # index of first particle
        self.j = j   # index of second particle
        self.k = k   # spring constant
        self.r0 = r0 # rest length
    
    def get_error(self, pos):
        """Get the error in the bond constraint"""
        d = pos[self.j] - pos[self.i]
        r = np.linalg.norm(d)
        return abs(r - self.r0)

class Molecule:
    def __init__(self, types, pos, bonds):
        self.natom = len(pos)
        self.types = types.copy()  # types of each atom
        self.particles = [Particle(pos[i], mass=types[i]) for i in range(len(pos))]
        self.pos   = pos   # list of integers of particle objects
        self.bonds = bonds           # list of Bond objects
        # -- to be done
        self.cog   = None            # np.array[3]
        self.inds  = np.zeros(len(pos), dtype=np.int32)   # list of indices of this p
        
    def update(self, bUpdateCenter):
        """Update center of geometry and radius of the molecule"""
        # Convert pos to numpy array if it's a list
        pos = np.array(self.pos)
        
        # Compute center of geometry
        if bUpdateCenter or self.cog is None:
            self.cog = np.mean(pos, axis=0)
        self.Rg = np.max(np.linalg.norm(pos - self.cog, axis=1))

class World:
    def __init__(self, dt=0.1, D=0.1):
        nmaxneighs = 8
        natoms     = 0
        # dynamicsl varialbes (for all particles in the system)
        self.apos      = None  # [natoms,2] {x,y} 
        self.pos_pred  = None  # [natoms,2] predicted positions
        self.vels      = None  # [natoms,2]
        self.forces    = None  # [natoms,2]
        # topplogy and parameters
        self.molecules    = []   # list of grups of atoms 
        
        self.neighs_bonds = None # [natoms,set()] # bond neighbors of each atom `i`; It is a set for fast 'in' query
        self.neighs       = None # [natoms,[]]    # all neighbors from both bonds and collisions 
        self.kngs         = None # [natoms,[]]    # stiffness   for each of the neighbors
        self.l0s          = None # [natoms,[]]    # rest length for each of the neighbors
        self.Aii          = None # [natoms]       # diagonal of the PD matrix
        
        # simulation parameters
        self.dt = dt               # time step
        self.D = D                 # collision offset
        self.Rc     = 1.0          # cutoff radius for non-bonded interactions
        self.Rg     = 2.0          # grouping radius for neighbor search
        self.K_col  = 10.0        # collision stiffness
        self.Kbond  = 10.0        # bond stiffness
        self.n_iter = 10           # number of iterations for constraint solver
        self.damping = 0.1         # velocity damping

    def allocate(self, natoms):
        self.apos      = np.zeros((natoms,2))
        self.pos_pred  = np.zeros((natoms,2))
        self.vels      = np.zeros((natoms,2))
        self.forces    = np.zeros((natoms,2))
        self.neighs_bonds = [set() for _ in range(natoms)]
        self.neighs       = [[]    for _ in range(natoms)]
        self.kngs         = [[]    for _ in range(natoms)]
        self.l0s          = [[]    for _ in range(natoms)]
        self.Aii          = np.zeros(natoms)
        self.masses       = np.ones(natoms)
    
    def add_molecule(self, molecule):
        """Add a molecule to the world"""
        self.molecules.append(molecule)
        self.from_molecules(self.molecules)

    def from_molecules(self, molecules):
        """
        Take all molecules from list of molecules and create world arrays 
        like apos, vels, forces, and bond-neighbors
        
        Args:
            molecules (list): List of Molecule objects to add to the world
        """
        # Set molecules in the world
        self.molecules = molecules
        
        # Count total number of atoms
        natoms = sum(len(mol.pos) for mol in molecules)
        
        # Allocate arrays for positions, velocities, and forces
        self.allocate(natoms)
        
        # Populate world arrays and bond neighbors
        atom_index = 0
        for mol in molecules:
            # Update molecule's indices
            mol.inds = np.arange(atom_index, atom_index + len(mol.pos))
            
            # Add molecule's atoms to world arrays
            for i, pos in enumerate(mol.pos):
                # Add position (x, y, z, mass)
                # For 2D simulation, set z=0 and use mass as 4th component
                self.apos[atom_index] = np.array([pos[0], pos[1]])
                
                # Add initial velocities (assuming zero initial velocity)
                self.vels[atom_index] = np.zeros(2)
                
                # Add initial forces (assuming zero initial force)
                self.forces[atom_index] = np.zeros(2)
                
                # Add bond neighbors for this atom
                for bond in mol.bonds:
                    if bond.i == i:
                        j_world = mol.inds[bond.j]  # Convert local to world index
                        self.neighs_bonds[atom_index].add((j_world, bond.r0))  # Store both neighbor index and rest length
                        #print(f"Added bond: {atom_index}-{j_world} with r0={bond.r0}")  # Debug
                    elif bond.j == i:
                        i_world = mol.inds[bond.i]  # Convert local to world index
                        self.neighs_bonds[atom_index].add((i_world, bond.r0))  # Store both neighbor index and rest length
                        #print(f"Added bond: {atom_index}-{i_world} with r0={bond.r0}")  # Debug
                
                atom_index += 1
        
        # Update groups (center of geometry and radius)
        self.update_groups()
    
    def update_groups(self, bUpdateCenter=False):
        '''
        Goes through atoms of molecule and find the center (optionaly) and group radius Rg
        '''
        for mol in self.molecules: mol.update(bUpdateCenter)
    
    def update_group_neighbors(self):
        """
        this is to accelerate the neighbor search. Each molecule keep its list of atoms which are closer than r<Rcg where Rcg=Rc+Rd and Rg is radius of the group (distance of furthers member atom from center of the group)  
          - this list is update less ofthen for performace
        """
        for im,mol in enumerate(self.molecules):
            cog  = mol.cog
            gngs = mol.mol_neighs
            for jm,molj in enumerate(self.molecules):
                if im == jm: continue
                for jp in molj:
                    if np.linalg.norm(self.apos[jp] - cog) < molj.Rg:
                        gngs[jp].add(im)

    def update_collision_neighbors(self, bCollisions=True):
        """
        This function goes through all the molecules and update the neighs_colls by adding all atoms which are closer than Rc from both the same molecule and tho other molecules
        - output is merged list of all neighbors ( both bonds and collisions) and their stiffness
        - this is accelerated by using the list  
        """
        for im, mol in enumerate(self.molecules):
            for ip in mol.inds:
                ngs  = []  # neighbors (both bonds and collisions)
                kngs = [] # stiffness for each neighbor
                l0s  = []  # rest length for each neighbor
                Kii = 0.0
                
                # First add bonded neighbors
                print(f"Bond neighbors: {self.neighs_bonds[ip]}")  # Debug
                for j, r0 in self.neighs_bonds[ip]:
                    ngs .append(j)
                    k = self.Kbond
                    kngs.append(k)
                    Kii += k
                    l0s .append(r0)
                    print(f"  Added bond neighbor {ip}-{j} with r0={r0}")  # Debug
                
                if bCollisions:
                    for jm, molj in enumerate(self.molecules):
                        for jp in molj.inds:
                            if ip == jp:
                                continue
                            if im == jm and jp in [j for j, _ in self.neighs_bonds[ip]]:
                                continue
                            if np.linalg.norm(self.apos[ip]- self.apos[jp]) < self.Rc:
                                ngs.append(jp)
                                k = self.K_col
                                kngs.append(k)
                                Kii += k
                                l0s.append(self.Rc)  # For collisions, rest length is Rc
                                print(f"  Added collision neighbor {ip}-{jp} with r0={self.Rc}")  # Debug

                # Store all neighbors and their parameters
                self.neighs[ip] = ngs
                self.kngs[ip]   = kngs
                self.l0s[ip]    = l0s
                self.Aii[ip]    = Kii   # A_ii = \sum_j K_{ij} + m_i /dt^2, for the moment we neglect the mass (inertial term m_i/dt^2)

        return self.neighs, self.kngs, self.l0s


    # --------- Iterative Solver  -------------
    def step(self):
        """Perform one time step of the simulation"""
        # Damp velocities
        self.vels[:] *= (1.0 - self.damping)
        self.predict_positions()
        self.solve_constraints()
        self.update_velocities()

    def predict_positions(self):
        """Step 1: Predict positions using current velocities"""
        for mol in self.molecules:
            for ip in mol.inds:
                self.pos_pred[ip] = self.apos[ip] + self.vels[ip] * self.dt

    def update_velocities(self):
        """Step 3: Update velocities from position changes"""
        for mol in self.molecules:
            for ip in mol.inds:
                self.vels[ip] = (self.apos[ip] - self.pos_pred[ip]) / self.dt

    def solve_constraints(self, niter=10):
        """Solve position constraints using Jacobi iteration"""
        # Update neighbors and constraints
        self.update_collision_neighbors()
        ps  = self.pos_pred.copy()
        b   = self.update_rhs(ps)
        Aii = self.update_Aii(ps)
        # Solve using Jacobi iteration
        for it in range(niter):
            ps, r = jacobi_iteration_sparse(ps, b, self.neighs, self.kngs, Aii)
            err_x = np.linalg.norm(r[:, 0])
            err_y = np.linalg.norm(r[:, 1])
            print(f"  errors: x={err_x:.3e} y={err_y:.3e}")
        self.pos_pred = ps.copy()
        return ps

    def update_rhs(self, pos_pred ):
        return make_PD_RHS( pos_pred, self.neighs, self.kngs, self.l0s, masses=self.masses, dt=1.0 )

    def update_Aii(self, pos_pred ):
        return update_PD_matrix( self.pos_pred, self.neighs, self.kngs, dt=1.0 )

    def print_bond_lengths(self, end=" "):
        for i in range(len(self.apos)):
            for j in self.neighs[i]:
                if j > i:  # Only print each distance once
                    d = np.linalg.norm(self.apos[j] - self.apos[i])
                    print(f"{i}-{j}: {d:.6f},", end=end )
        #print()

# ========================================================
# ============ Create initial configurations  ============
# ========================================================

def create_triangle_system(l=1.0):
    """Create three particles in equilateral triangle with bonds"""
    # Create three particles at vertices of equilateral triangle
    pos = [
        [0.0, 0.0],                # First particle at origin
        [l, 0.0],                  # Second particle at distance l along x-axis
        [l/2, l*np.sqrt(3)/2]      # Third particle to form equilateral triangle
    ]
    types = [1.0, 1.0, 1.0]  # Equal masses
    bonds = [
        Bond(0, 1, r0=l),  # Bond between particles 0 and 1
        Bond(1, 2, r0=l),  # Bond between particles 1 and 2
        Bond(2, 0, r0=l)   # Bond between particles 2 and 0
    ]
    return Molecule(types, pos, bonds)

def create_water(pos=[0,0], l0=1.0, ang=np.pi/2, vel=[0,0] ):
    """Create a water molecule (H2O) at given position with given velocity"""
    types = [16.0, 1.0, 1.0]  # masses of O, H1, H2

    ca = np.cos(ang)
    sa = np.sin(ang)
    positions = [
        [pos[0],       pos[1]      ],  # O
        [pos[0]+l0,    pos[1]      ],  # H1
        [pos[0]+l0*ca, pos[1]+l0*sa]   # H2
    ]
    bonds = [
        Bond(0, 1, k=10.0, r0=1.0),  # O-H1 bond
        Bond(0, 2, k=10.0, r0=1.0),  # O-H2 bond
    ]
    return Molecule(types, positions, bonds)

def create_water_grid(nx=3, ny=3, spacing=2.0, l0=1.0, ang=np.pi/2):
    """Create a grid of water molecules"""
    molecules = []
    for i in range(nx):
        for j in range(ny):
            pos = [i*spacing, j*spacing]
            molecules.append(create_water(pos, l0=l0, ang=ang))
    return molecules

def init_triangle_world(l=1.0):
    """Initialize world with three particles in triangle"""
    world = World(dt=0.05, D=0.1)  # Create world with small time step
    mol = create_triangle_system(l)
    world.add_molecule(mol)
    return world

def init_world( nx=2, ny=2, spacing=2.0, l0=1.0, ang=np.pi/2 ):
    """Initialize world with water molecules"""
    world = World(dt=0.01)
    mols = create_water_grid(nx, ny, spacing=spacing, l0=l0, ang=ang)
    world.from_molecules(mols)
    return world

# ========================================================
# ============ Create initial configurations  ============
# ========================================================

def test_jacobi_debug(world, n_iter=5, bPrint=True, bPlot=True, sz=5, fRnadom=0.3, bUpdateB=True ):
    """Test convergence with three particles and debug visualization"""
    print("\nTesting three-particle system convergence with debug:")
    
    world.update_collision_neighbors(bCollisions=False)  # Comment this out

    # random perturbation to points
    world.apos[:,0] += np.random.uniform(-fRnadom, fRnadom, len(world.apos))

    #b = make_PD_RHS(world.apos, world.neighs, world.kngs, world.l0s, masses=world.masses, dt=world.dt)
    b = make_PD_RHS(world.apos, world.neighs, world.kngs, world.l0s )

    print("##### Initial world state: \n")
    for i in range(len(world.apos)):
        print( f"i: {i} Aii: {world.Aii[i]} b: {b[i]} ")  

    if bPlot:
        plt.figure(figsize=(n_iter*sz,sz))

    # Solve iteratively the equation Ap = b, where A is the PD matrix, p is the position, and b is the right-hand side from make_PD_RHS
    ps = world.apos  # Only use x,y components
    for it in range(n_iter):

        if bPlot:
            ax = plt.subplot(1, n_iter, it+1)

        if bUpdateB:
            # NOTE: this goes beyond linear-algebra solution, in poper linear algebra we should not update b during iterative solution of Ap=b    
            #       but with this it converge faster, and it is more physical
            b = make_PD_RHS(world.apos, world.neighs, world.kngs, world.l0s )

        ps_new, r = jacobi_iteration_sparse_debug(ps, b, world.neighs, world.kngs, world.Aii, l0s=world.l0s, bPlot=bPlot)
        ps[:] = ps_new[:]
        if bPrint:
            print(f"\nIteration {it} Distances:", end=" " )
            world.print_bond_lengths( end=" " )

def test_triangle_cl():
    """Compare Python and OpenCL implementations for triangle system"""
    # Initialize world with three particles
    world = init_triangle_world()
    n = len(world.apos)

    import constrains_cl
    
    solver = constrains_cl.CLConstrains(max_points=n, max_neighs=2)
    
    # Prepare neighbor data (only bonds, no collisions)
    world.update_collision_neighbors(bCollisions=False)
    
    # Print and perturb initial positions
    print("Initial positions:")
    print(world.apos)
    world.apos += np.random.rand(n,2) * 0.3
    print("\nPerturbed positions:")
    print(world.apos)
    
    # Parameters for both solvers
    dt = 0.1
    inv_dt2 = 1.0/(dt*dt)
    Rd = 0.1
    n_iter = 10
    
    # Run both solvers
    print("\nRunning Python solver...")
    world.solve_constraints(niter=n_iter)
    pos_py = world.apos.copy()
    
    print("\nRunning OpenCL solver...")
    pos_cl = solver.solve_constraints(world.apos, np.ones(n), world.neighs, world.kngs, world.l0s, inv_dt2, Rd, niter=n_iter)
    
    # Compare results
    diff = np.linalg.norm(pos_py - pos_cl)
    print(f"\nDifference between solutions: {diff}")
    
    # Visualize results
    plt.figure(figsize=(12,4))
    plt.subplot(131);  pu.plot_molecules(world, ax=plt.gca(), color='blue',  label='Initial'); plt.title('Initial')
    world.apos = pos_py; plt.subplot(132);  pu.plot_molecules(world, ax=plt.gca(), color='red',   label='Python');  plt.title('Python Solution')
    world.apos = pos_cl; plt.subplot(133);  pu.plot_molecules(world, ax=plt.gca(), color='green', label='OpenCL');  plt.title('OpenCL Solution')
    plt.tight_layout()

def test_1(world, n_iter = 10):
    """Test convergence of constraint solver"""
    print("\nTesting constraint solver convergence:")
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    pu.plot_molecules(world, ax=ax1, color='b', label='Initial', markersize=10)
    errors = []
    
    # Plot intermediate states with different colors
    colors = plt.cm.viridis(np.linspace(0, 1, n_iter))
    def plot_callback(i, x, y, err_x, err_y):
        pu.plot_molecules(world, x=x, y=y, ax=ax1, color=colors[i], alpha=0.5, 
                         label=f'iter {i}', markersize=5)
        ax1.set_xlim(-0.5, 1.5)
        ax1.set_ylim(-0.5, 1.5)
    
    world.solve_constraints(n_iter, errs=errors, callback=plot_callback)
    pu.plot_molecules(world, ax=ax1, color='r', label='Final', markersize=10)
    ax1.legend()
    ax1.set_title('Molecule Positions')
    ax1.grid(True)
    
    pu.plot_convergence(errors, ax=ax2)
    ax2.set_title('Constraint Error')
    ax2.grid(True)
    plt.tight_layout()
    plt.show()

def test_2(world):
    """Test dynamics over multiple steps"""
    print("\nTesting dynamics:")
    
    # Setup visualization
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    errors = []
    
    # Plot initial state
    pu.plot_molecules(world, ax=ax1, color='b', label='Initial')
    
    # Run simulation steps
    n_steps = 10
    for i in range(n_steps):
        # Step simulation
        world.step()
        err = world.solve_constraints()  # Get error from last constraint solve
        errors.append(err)
        print(f"Step {i}: error = {err}")
        
        # Plot intermediate state
        if i % 2 == 0:
            pu.plot_molecules(world, ax=ax1, color='k', alpha=0.2)
    
    # Plot final state and convergence
    pu.plot_molecules(world, ax=ax1, color='r', label='Final')
    ax1.legend()
    ax1.set_title('Molecule Positions')
    
    pu.plot_convergence(errors, ax=ax2)
    ax2.set_title('Constraint Error')
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import plot_utils as pu
    
    # Initialize world with triangle system
    world = init_triangle_world()
    
    # Run test with debug visualization
    test_jacobi_debug(world)
    plt.show()