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
        self.pos   = pos   # list of position lists [[x,y], ...]
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
        sum_j  = 0  # Off-diagonal contributions
        for j, k in zip(neighs[i], kngs[i]):
            sum_j -= k * x[j]   # Aij = -k_ij
        x_out[i] =  (b[i] - sum_j) / Aii[i]   # solution x_new = (b - sum_(j!=i){ Aij * x[j] } ) / Aii
        y_i = Aii[i]*x[i] + sum_j              # This is Ax_i
        r[i] = b[i] - y_i                      # Residual r = b - Ax
        print(f"i: {i} r_i: {r[i]:.6f} pi_: {b[i]/Aii[i]:.6f} x_i: {x[i]:.6f}")
    return x_out, r

class World:
    def __init__(self, dt=0.1, D=0.1):
        nmaxneighs = 8
        natoms     = 0
        # dynamic variables (for all particles in the system)
        self.apos      = None  # [natoms,4] {x,y,z,m} 
        self.pos_pred  = None  # [natoms,4] predicted positions
        self.vels      = None  # [natoms,4]
        self.forces    = None  # [natoms,4]
        # topology and parameters
        self.molecules    = []   # list of groups of atoms 
        
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

        self.tol = 1e-6            # convergence tolerance

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
        Goes through atoms of molecule and find the center (optionally) and group radius Rg
        '''
        for mol in self.molecules: mol.update(bUpdateCenter)
    
    def update_group_neighbors(self):
        """
        This is to accelerate the neighbor search. Each molecule keeps its list of atoms which are closer than r<Rcg where Rcg=Rc+Rd and Rg is radius of the group (distance of furthest member atom from center of the group)  
          - this list is updated less often for performance
        """
        for im,mol in enumerate(self.molecules):
            cog  = mol.cog
            gngs = mol.mol_neighs
            for jm,molj in enumerate(self.molecules):
                if im == jm: continue
                for jp in molj.inds:
                    if np.linalg.norm(self.apos[jp][:3] - cog[:3]) < molj.Rg:
                        gngs[jp].add(im)

    def update_collision_neighbors(self):
        """
        This function goes through all the molecules and updates the neighs_colls by adding all atoms which are closer than Rc from both the same molecule and other molecules
        - Output is a merged list of all neighbors ( both bonds and collisions) and their stiffness
        """
        for im, mol in enumerate(self.molecules):
            for ip in mol.inds:
                ngs = []   # neighbors (both bonds and collisions)
                kngs = []  # stiffness for each neighbor
                r0s = []   # rest length for each neighbor
                
                # First add bonded neighbors
                print(f"Bond neighbors for atom {ip}: {self.neighs_bonds[ip]}")  # Debug
                for j, r0 in self.neighs_bonds[ip]:
                    # **Removed the `continue` statement to include bond constraints**
                    ngs.append(j)
                    kngs.append(self.Kbond)
                    r0s.append(r0)
                    print(f"  Added bond neighbor {ip}-{j} with r0={r0}")  # Debug
                
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
                            print(f"  Added collision neighbor {ip}-{jp} with r0={self.Rc}")  # Debug

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

    def solve_constraints(self, niter=10, errs=None, callback=None):
        """Step 2: Solve position constraints (bonds and collisions)"""
        neighs, kngs, r0s = self.update_collision_neighbors()  # Update neighbors

        Aii = self.update_PD_matrix(neighs, kngs)  # Compute Aii
        b   = self.update_rhs(neighs, kngs, r0s)    # Compute RHS with corrected terms

        # Separate x and y components
        bx = b[:,0]
        by = b[:,1]

        # Initial guess is current positions
        x = self.apos[:,0].copy()
        y = self.apos[:,1].copy()

        print("\nSolving for x component:")
        # Iterate to solve constraints for component x
        for itr in range(niter):

            print(f"\nIteration {itr}:")
            print(f"  x positions: {x}")
            print(f"  y positions: {y}")
            
            x, err_x = jacobi_iteration_sparse(x, bx, neighs, kngs, Aii ) # solve for x component
            y, err_y = jacobi_iteration_sparse(y, by, neighs, kngs, Aii ) # solve for y component

            print(f"  errors: x={np.linalg.norm(err_x):.3e} y={np.linalg.norm(err_y):.3e}")

            if errs is not None:
                errs.append(np.linalg.norm(err_x) + np.linalg.norm(err_y))

            if callback is not None:
                callback(itr, x,y, err_x, err_y)
            

        # Update positions
        self.apos[:,0] = x
        self.apos[:,1] = y

        return np.linalg.norm(err_x) + np.linalg.norm(err_y)

    def update_rhs(self, neighs, kngs, r0s):
        """
        Update the right-hand side of the Projective Dynamics equation:
        b_i = (m_i/dt^2)p'_i + sum_j (K_ij * (p_j + d_ij))
        where d_ij = (p_i - p_j) / |p_i - p_j| * r0
        """
        n = len(self.apos)
        b = np.zeros_like(self.apos[:,:2])  # Only x,y coordinates

        for ip in range(n):
            # Inertial term: (m_i/dt^2) * p'_i
            b[ip] = (self.apos[ip,3] / (self.dt**2)) * self.pos_pred[ip,:2]
            print(f"\nRHS for atom {ip}:")
            print(f"  Inertial term: {b[ip]}")
            
            # Sum of constraint terms
            for j, k, r0 in zip(neighs[ip], kngs[ip], r0s[ip]):
                # Current positions
                p_i = self.apos[ip,:2]
                p_j = self.apos[j,:2]
                
                # Vector from i to j
                r_ij = p_j - p_i
                r = np.linalg.norm(r_ij)
                
                if r > 1e-12:  # Avoid division by zero
                    dir_ij = (p_i - p_j) / r  # **Reversed direction**
                    d_ij = dir_ij * r0  # Desired displacement
                    
                    # Include p_j in the RHS
                    p_ij_desired = p_j + d_ij
                    force = k * p_ij_desired
                    b[ip] += force
                    
                    # Debugging output
                    print(f"  Neighbor {j}: r={r:.3f} r0={r0:.3f} k={k:.1f}")
                    print(f"    p_j={p_j} d_ij={d_ij} p_ij_desired={p_ij_desired} force={force}")
        
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

def create_two_particle_system(l=2.0):
    """Create a simple two-particle system bonded with r0=1.0"""
    types = [1.0, 1.0]  # masses
    positions = [
        [0.0, 0.0],  # Particle 0
        [2.0, 0.0]   # Particle 1
    ]
    bonds = [
        Bond(0, 1, k=10.0, r0=1.0)  # Bond between Particle 0 and 1
    ]
    return Molecule(types, positions, bonds)

def init_two_particle_world( l=2.0):
    """Initialize world with a two-particle bonded system"""
    world = World(dt=0.1)
    two_particle = create_two_particle_system()
    world.from_molecules([two_particle])
    return world

def test_two_particle(world, n_iter=10):
    """Test convergence with a two-particle system"""
    print("\nTesting two-particle system convergence:")
    
    for itr in range(n_iter):
        print(f"\nIteration {itr}:")
        world.solve_constraints(niter=1)  # One iteration at a time
        
        # Print positions
        for ip, p in enumerate(world.apos[:,:2]):
            print(f"Particle {ip}: Position = {p}")
        
        # Check distance
        p0 = world.apos[0,:2]
        p1 = world.apos[1,:2]
        dist = np.linalg.norm(p1 - p0)
        print(f"Distance between Particle 0 and 1: {dist:.6f}")
        
        if dist >= 1.0 - 1e-4 and dist <= 1.0 + 1e-4:
            print("Converged to desired bond length.")
            break

if __name__ == "__main__":
    world = init_two_particle_world( l=3.123432)
    test_two_particle(world, n_iter=10)
