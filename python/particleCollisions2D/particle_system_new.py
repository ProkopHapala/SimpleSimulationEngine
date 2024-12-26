import numpy as np

        
class Bond:
    def __init__(self, i, j, k, r0):
        self.i = i          # index of first particle
        self.j = j          # index of second particle
        self.k = k          # spring constant
        self.r0 = r0        # equilibrium length

class Molecule:
    def __init__(self, types, pos, bonds):
        self.natom = len(pos)
        self.types =  # types of each atom
        self.pos   = pos   # list of integers of particle objects
        self.bonds = bonds           # list of Bond objects
        # -- to be done
        self.cog   = None            # np.array[3]
        self.inds  = np.zeros(len(pos), dtype=np.int32)   # list of indices of this p
        
    def update(self, bUpdateCenter ):
        """Update center of geometry and radius of the molecule"""
        # Compute center of geometry
        if bUpdateCenter:
            self.cog = np.mean(self.pos, axis=0)
        self.Rg = np.max(np.linalg.norm(self.pos - self.cog, axis=1))

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
        self.apos   = None  # [natoms,4] {x,y,z,m} 
        self.vels   = None  # [natoms,4]
        self.forces = None  # [natoms,4]
        # topplogy and parameters
        self.molecules    = []   # list of grups of atoms 
        
        self.neighs_bonds = None # [ set() ]*natoms 
        #self.neighs_bonds = [natoms, nmaxneighs]   # neighbor lists for each atom due to bonds ( not changed during the simulation)
        #self.neighs_colls = [natoms, nmaxneighs]   # neighbor lists for collisions ( collision changes during the simulation)
        self.neighs_all   = None #[natoms, nmaxneighs]   # include neighbors from bonds and collisions 
        self.kngs         = None #[natoms, nmaxneighs]   # stiffness for each of the neighbors

        # simulation parameters
        self.dt = dt               # time step
        self.D = D                 # collision offset
        self.Rc     = 1.0          # cutoff radius for non-bonded interactions
        self.Rg     = 2.0          # grouping radius for neighbor search
        self.K_col  = 100.0        # collision stiffness
        self.n_iter = 10           # number of iterations for constraint solver
        
    def allocate(self, natoms):
        self.apos   = np.zeros((natoms,4))
        self.vels   = np.zeros((natoms,4))
        self.forces = np.zeros((natoms,4))
    
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
        
        # Initialize bond neighbors list
        self.neighs_bonds = [set() for _ in range(natoms)]
        
        # Populate world arrays and bond neighbors
        atom_index = 0
        for mol in molecules:
            # Update molecule's indices
            mol.inds = np.arange(atom_index, atom_index + len(mol.pos))
            
            # Add molecule's atoms to world arrays
            for i, pos in enumerate(mol.pos):
                # Add position (x, y, z, mass)
                self.apos[atom_index] = np.concatenate([pos, [mol.types[i]]])
                
                # Add initial velocities (assuming zero initial velocity)
                self.vels[atom_index] = np.zeros(4)
                
                # Add initial forces (assuming zero initial force)
                self.forces[atom_index] = np.zeros(4)
                
                # Add bond neighbors for this atom
                for bond in mol.bonds:
                    if bond.i == i:
                        self.neighs_bonds[atom_index].add(mol.inds[bond.j])
                    elif bond.j == i:
                        self.neighs_bonds[atom_index].add(mol.inds[bond.i])
                
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
                ngs = [] 
                kngs = []
                bngs = self.neighs_bonds[ip]
                # find collision neighbors of atom ip (which are not bonded to it)
                for jm, molj in enumerate(self.molecules):
                    for jp in molj.inds:
                        if ip == jp:
                            continue
                        if im == jm and jp in bngs:
                            continue
                        if np.linalg.norm(self.apos[ip][:2] - self.apos[jp][:2]) < self.Rc:
                            ngs.append(jp)
                            kngs.append(self.K_col)

                # store both bond and collision neighbors
                self.neighs[ip] = list(bngs) + ngs 
                self.kngs[ip] = [self.Kbond]*len(bngs) + kngs

        return self.neighs, self.kngs


    # --------- Iterative Solver  -------------

    def step(self):
        """Perform one time step of the simulation"""
        self._predict_positions()
        self._solve_constraints()
        self._update_velocities()

    def _solve_constraints(self):
        """Step 2: Solve position constraints (bonds and collisions)"""
        neighs, kngs = self.update_collision_neighbors()  # 

        Aii = self.update_PD_matrix( neighs, kngs )
        b   = self.update_rhs      ( neighs, kngs )
        
        bx = b[:,0]
        by = b[:,1]
        # Iterate to solve constraints for component x
        for _ in range(self.n_iter):
            x, err_x = jacobi_iteration_sparse(x, bx, neighs, kngs, Aii ) # solve for x component
            y, err_y = jacobi_iteration_sparse(y, by, neighs, kngs, Aii ) # solver for y component
        self.apos[:,:2] = np.column_stack((x, y))

    def _update_velocities(self):
        """Step 3: Update velocities from position changes"""
        for mol in self.molecules:
            for p in mol.particles:
                p.vel = (p.pos - p.pos_pred) / self.dt
                p.pos_pred = None  # cleanup
                


    def _predict_positions(self):
        """Step 1: Predict positions using current velocities"""
        for mol in self.molecules:
            for p in mol.particles:
                p.pos_pred = p.pos + p.vel * self.dt
                p.pos      = p.pos_pred.copy()  # Initialize position for constraint solving
                
    def update_PD_matrix(self, neighs, kngs):
        """
        Update the diagonal matrix for Projective Dynamics
        
        Args:
            neighs (list): List of neighbor indices for each particle
            kngs (list): List of corresponding spring/collision stiffnesses
        
        Returns:
            numpy array: Diagonal matrix of coefficients
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

    def update_rhs(self, neighs, kngs):
        """
        Update the right-hand side of the Projective Dynamics equation
        
        Args:
            neighs (list): List of neighbor indices for each particle
            kngs (list): List of corresponding spring/collision stiffnesses
        
        Returns:
            numpy array: Right-hand side vector for each particle
        """
        n = len(self.apos)
        b = np.zeros_like(self.apos[:,:2])  # Only x,y coordinates
        
        for ip in range(n):
            # Inertial term: (m_i/dt^2) * p'_i
            b[ip] = (self.apos[ip,3] / (self.dt**2)) * self.apos[ip,:2]
            
            # Sum of constraint displacements
            for j, k in zip(neighs[ip], kngs[ip]):
                # Compute displacement to satisfy constraint
                d_ij = self._compute_constraint_displacement(ip, j)
                b[ip] += k * d_ij
        
        return b


# Example usage:
if __name__ == "__main__":
    #import numpy as np
    import matplotlib.pyplot as plt
    #from particle_system import Particle, Bond, Molecule, World
    from visualizer import Visualizer

    def create_water(pos=[0,0], vel=[0,0]):
        """Create a water molecule (H2O) at given position with given velocity"""
        p1 = Particle([pos[0], pos[1]], vel=vel, mass=16.0)      # O
        p2 = Particle([pos[0]+1, pos[1]], vel=vel, mass=1.0)     # H1
        p3 = Particle([pos[0], pos[1]+1], vel=vel, mass=1.0)     # H2
        bonds = [
            Bond(0, 1, k=100.0, r0=1.0),  # O-H1 bond
            Bond(0, 2, k=100.0, r0=1.0),  # O-H2 bond
        ]
        return Molecule([p1, p2, p3], bonds)

    def create_water_grid(n=3, spacing=2.0):
        """Create a grid of water molecules"""
        molecules = []
        for i in range(n):
            for j in range(n):
                pos = [i*spacing, j*spacing]
                vel = [0.1*(i-n/2), 0.1*(j-n/2)]  # velocity depends on position
                molecules.append(create_water(pos, vel))
        return molecules

    # Create simulation world with multiple water molecules
    dt = 0.01
    world = World(dt=dt, D=0.1)
    world.K_col = 1000.0  # Increase collision stiffness
    molecules = create_water_grid(n=3, spacing=1.3)
    for mol in molecules:
        world.add_molecule(mol)
    
    # Setup visualization
    vis = Visualizer(world)
    
    def callback(itr, x, r):
        """Callback to update visualization during iteration"""
        err = np.linalg.norm(r)
        print(f"Iteration {itr:4d} residual: {err:10.6f}")
        if itr % 5 == 0:  # Update visualization every 5 iterations
            vis.update()
    
    # Solve one step of projective dynamics
    print("\nStarting Jacobi iteration:")
    errs = world.solve_pd_step(niter_outer=10, niter_inner=10, tol=1e-8, callback=callback)
    
    # Plot convergence
    plt.figure(figsize=(10,6))
    plt.semilogy(errs, 'b-', label='Jacobi')
    plt.grid(True)
    plt.xlabel('Iteration')
    plt.ylabel('Error (log scale)')
    plt.title('Convergence of Jacobi Method for Water Molecules')
    plt.legend()
    plt.show()