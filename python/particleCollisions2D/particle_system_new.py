import numpy as np

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
        d = pos[self.j,:2] - pos[self.i,:2]
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
            sum_j += k * x[j]   # Off-diagonal contribution
        x_out[i] =  (b[i] + sum_j) / Aii[i]   # solution x_new = (b - sum_(j!=i){ Aij * x[j] } ) / Aii
        r[i]     = b[i] +sum_j - Aii[i]*x[i] # Residual r = b - Ax ;  Ax = Aii * x[i] + sum_(j!=i){ Aij * x[j] }
        #print("CPU i: %i Aii[i]: %f b[i]: %f sum_j: %f x_out[i]: %f r[i]: %f" %(i, Aii[i], b[i], sum_j, x_out[i], r[i]) );
    return x_out, r

class World:
    def __init__(self, dt=0.1, D=0.1):
        nmaxneighs = 8
        natoms     = 0
        # dynamicsl varialbes (for all particles in the system)
        self.apos      = None  # [natoms,4] {x,y,z,m} 
        self.pos_pred  = None  # [natoms,4] predicted positions
        self.vels      = None  # [natoms,4]
        self.forces    = None  # [natoms,4]
        # topplogy and parameters
        self.molecules    = []   # list of grups of atoms 
        
        self.neighs_bonds = None # [ set() ]*natoms 
        self.neighs      = None # [ [] ]*natoms   # include neighbors from bonds and collisions 
        self.kngs       = None # [ [] ]*natoms   # stiffness for each of the neighbors
        self.r0s        = None # [ [] ]*natoms   # rest length for each of the neighbors
        self.Kbond      = 100.0  # bond stiffness

        # simulation parameters
        self.dt = dt               # time step
        self.D = D                 # collision offset
        self.Rc     = 1.0          # cutoff radius for non-bonded interactions
        self.Rg     = 2.0          # grouping radius for neighbor search
        self.K_col  = 100.0        # collision stiffness
        self.n_iter = 10           # number of iterations for constraint solver
        self.damping = 0.1         # velocity damping

    def allocate(self, natoms):
        self.apos      = np.zeros((natoms,4))
        self.pos_pred  = np.zeros((natoms,4))
        self.vels      = np.zeros((natoms,4))
        self.forces    = np.zeros((natoms,4))
        self.neighs_bonds = [set() for _ in range(natoms)]
        self.neighs       = [[] for _ in range(natoms)]
        self.kngs        = [[] for _ in range(natoms)]
        self.r0s         = [[] for _ in range(natoms)]
    
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
                self.apos[atom_index] = np.array([pos[0], pos[1], 0.0, mol.types[i]])
                
                # Add initial velocities (assuming zero initial velocity)
                self.vels[atom_index] = np.zeros(4)
                
                # Add initial forces (assuming zero initial force)
                self.forces[atom_index] = np.zeros(4)
                
                # Add bond neighbors for this atom
                for bond in mol.bonds:
                    if bond.i == i:
                        j_world = mol.inds[bond.j]  # Convert local to world index
                        self.neighs_bonds[atom_index].add((j_world, bond.r0))  # Store both neighbor index and rest length
                        print(f"Added bond: {atom_index}-{j_world} with r0={bond.r0}")  # Debug
                    elif bond.j == i:
                        i_world = mol.inds[bond.i]  # Convert local to world index
                        self.neighs_bonds[atom_index].add((i_world, bond.r0))  # Store both neighbor index and rest length
                        print(f"Added bond: {atom_index}-{i_world} with r0={bond.r0}")  # Debug
                
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
                    if np.linalg.norm(self.apos[jp][:3] - cog[:3]) < molj.Rg:
                        gngs[jp].add(im)

    def update_collision_neighbors(self):
        """
        This function goes through all the molecules and update the neighs_colls by adding all atoms which are closer than Rc from both the same molecule and tho other molecules
        - output is merged list of all neighbors ( both bonds and collisions) and their stiffness
        - this is accelerated by using the list  
        """
        for im, mol in enumerate(self.molecules):
            for ip in mol.inds:
                ngs = []  # neighbors (both bonds and collisions)
                kngs = [] # stiffness for each neighbor
                r0s = []  # rest length for each neighbor
                
                # First add bonded neighbors
                print(f"\nProcessing atom {ip}:")  # Debug
                print(f"Bond neighbors: {self.neighs_bonds[ip]}")  # Debug
                for j, r0 in self.neighs_bonds[ip]:
                    ngs.append(j)
                    kngs.append(self.Kbond)
                    r0s.append(r0)
                    print(f"  Added bond neighbor {j} with r0={r0}")  # Debug
                
                # Then find collision neighbors
                for jm, molj in enumerate(self.molecules):
                    for jp in molj.inds:
                        if ip == jp:
                            continue
                        if im == jm and jp in [j for j, _ in self.neighs_bonds[ip]]:
                            continue
                        if np.linalg.norm(self.apos[ip][:2] - self.apos[jp][:2]) < self.Rc:
                            ngs.append(jp)
                            kngs.append(self.K_col)
                            r0s.append(self.Rc)  # For collisions, rest length is Rc
                            print(f"  Added collision neighbor {jp} with r0={self.Rc}")  # Debug

                # Store all neighbors and their parameters
                self.neighs[ip] = ngs
                self.kngs[ip] = kngs
                self.r0s[ip] = r0s

        return self.neighs, self.kngs, self.r0s


    # --------- Iterative Solver  -------------
    def step(self):
        """Perform one time step of the simulation"""
        # Damp velocities
        self.vels[:,:2] *= (1.0 - self.damping)
        
        self.predict_positions()
        self.solve_constraints()
        self.update_velocities()

    def predict_positions(self):
        """Step 1: Predict positions using current velocities"""
        for mol in self.molecules:
            for ip in mol.inds:
                self.pos_pred[ip,:2] = self.apos[ip,:2] + self.vels[ip,:2] * self.dt

    def update_velocities(self):
        """Step 3: Update velocities from position changes"""
        for mol in self.molecules:
            for ip in mol.inds:
                self.vels[ip,:2] = (self.apos[ip,:2] - self.pos_pred[ip,:2]) / self.dt

    def solve_constraints(self):
        """Step 2: Solve position constraints (bonds and collisions)"""
        neighs, kngs, r0s = self.update_collision_neighbors()  # 

        Aii = self.update_PD_matrix( neighs, kngs )
        b   = self.update_rhs      ( neighs, kngs, r0s )
        
        # Separate x and y components
        bx = b[:,0]
        by = b[:,1]
        
        # Initial guess is current positions
        x = self.apos[:,0].copy()
        y = self.apos[:,1].copy()
        
        # Iterate to solve constraints for component x
        for _ in range(self.n_iter):
            x, err_x = jacobi_iteration_sparse(x, bx, neighs, kngs, Aii ) # solve for x component
            y, err_y = jacobi_iteration_sparse(y, by, neighs, kngs, Aii ) # solver for y component
        
        # Update positions
        self.apos[:,0] = x
        self.apos[:,1] = y

        return np.linalg.norm(err_x) + np.linalg.norm(err_y)

    def update_rhs(self, neighs, kngs, r0s):
        """
        Update the right-hand side of the Projective Dynamics equation:
        b_i = (m_i/dt^2)p'_i + sum_j (K_ij * d_ij)
        where d_ij is the displacement vector (length r0 for each neighbor)
        """
        n = len(self.apos)
        b = np.zeros_like(self.apos[:,:2])  # Only x,y coordinates
        
        for ip in range(n):
            # Inertial term: (m_i/dt^2) * p'_i
            b[ip] = (self.apos[ip,3] / (self.dt**2)) * self.pos_pred[ip,:2]
            
            # Sum of constraint terms
            for j, k, r0 in zip(neighs[ip], kngs[ip], r0s[ip]):
                # Current positions
                p_i = self.apos[ip,:2]
                p_j = self.apos[j,:2]
                
                # Vector from i to j
                r_ij = p_j - p_i
                r = np.linalg.norm(r_ij)
                
                if r > 0:  # Avoid division by zero
                    dir_ij = r_ij / r
                    d_ij = dir_ij * r0  # Displacement vector of length r0
                    b[ip] += k * d_ij  # Add K_ij * d_ij to RHS
        
        return b

    def update_PD_matrix(self, neighs, kngs):
        """
        Update the diagonal matrix for Projective Dynamics
        """
        n = len(self.apos)
        Aii = np.zeros(n)
        
        for ip in range(n):
            # Inertial term: m_i/dt^2
            Aii[ip] = self.apos[ip,3] / (self.dt**2)
            
            # Sum of stiffness for all neighbors
            for j, k in zip(neighs[ip], kngs[ip]):
                Aii[ip] += k
        
        return Aii


def create_water(pos=[0,0], vel=[0,0]):
    """Create a water molecule (H2O) at given position with given velocity"""
    types = [16.0, 1.0, 1.0]  # masses of O, H1, H2
    positions = [
        [pos[0], pos[1]],      # O
        [pos[0]+1, pos[1]],    # H1
        [pos[0], pos[1]+1]     # H2
    ]
    bonds = [
        Bond(0, 1, k=10.0, r0=1.0),  # O-H1 bond
        Bond(0, 2, k=10.0, r0=1.0),  # O-H2 bond
    ]
    return Molecule(types, positions, bonds)

def create_water_grid(n=3, spacing=2.0):
    """Create a grid of water molecules"""
    molecules = []
    for i in range(n):
        for j in range(n):
            pos = [i*spacing, j*spacing]
            molecules.append(create_water(pos))
    return molecules

def init_world():
    """Initialize world with water molecules"""
    world = World(dt=0.01)
    mols = create_water_grid(2, 2)
    world.from_molecules(mols)
    return world

def test_1(world):
    """Test convergence of constraint solver"""
    print("\nTesting constraint solver convergence:")
    
    # Setup visualization
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    errors = []
    
    # Plot initial state
    pu.plot_molecules(world, ax=ax1, color='b', label='Initial')
    
    # Run constraint solver iterations
    n_iter = 20
    for i in range(n_iter):
        # Solve constraints
        err = world.solve_constraints()
        errors.append(err)
        print(f"Iteration {i}: error = {err}")
        
        # Plot intermediate state
        if i % 5 == 0:
            pu.plot_molecules(world, ax=ax1, color='k', alpha=0.2)
    
    # Plot final state and convergence
    pu.plot_molecules(world, ax=ax1, color='r', label='Final')
    ax1.legend()
    ax1.set_title('Molecule Positions')
    
    pu.plot_convergence(errors, ax=ax2)
    ax2.set_title('Constraint Error')
    
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
    
    world = init_world()
    test_1(world)  # Test constraint solver convergence

    #world = init_world()  # Re-initialize for test 2
    #test_2(world)  # Test dynamics