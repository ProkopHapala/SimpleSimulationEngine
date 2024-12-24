import numpy as np

class Particle:
    def __init__(self, pos, vel=None, mass=1.0):
        self.pos = np.array(pos, dtype=np.float64)    # position vector [x,y]
        self.vel = np.zeros(2) if vel is None else np.array(vel, dtype=np.float64)   # velocity vector [vx,vy]
        self.m = mass                                 # mass
        self.f = np.zeros(2)                         # force accumulator
        self.neighbors = []                          # list of neighbor particles for collision detection
        
class Bond:
    def __init__(self, i, j, k, r0):
        self.i = i          # index of first particle
        self.j = j          # index of second particle
        self.k = k          # spring constant
        self.r0 = r0        # equilibrium length

class Molecule:
    def __init__(self, particles, bonds):
        self.particles = particles    # list of Particle objects
        self.bonds = bonds           # list of Bond objects
        self._update_com()           # center of mass
        
    def _update_com(self):
        """Update center of mass of the molecule"""
        total_mass = sum(p.m for p in self.particles)
        self.com = np.sum([p.pos * p.m for p in self.particles], axis=0) / total_mass

class World:
    def __init__(self, dt=0.1, D=0.1):
        self.molecules = []          # list of molecules
        self.dt = dt                # time step
        self.D = D                  # collision offset
        self.Rc = 1.0              # cutoff radius for non-bonded interactions
        self.Rg = 2.0              # grouping radius for neighbor search
        self.K_col = 100.0         # collision stiffness
        self.n_iter = 10           # number of iterations for constraint solver
        
    def add_molecule(self, molecule):
        self.molecules.append(molecule)

    def _update_neighbors(self):
        """Update neighbor lists for collision detection using molecule-based grouping"""
        # Clear old neighbors
        for mol in self.molecules:
            for p in mol.particles:
                p.neighbors = []

        # First check collisions within each molecule (non-bonded atoms)
        for mol in self.molecules:
            n = len(mol.particles)
            # Create set of bonded pairs for quick lookup
            bonded_pairs = set((b.i, b.j) for b in mol.bonds) | set((b.j, b.i) for b in mol.bonds)
            # Check all pairs in molecule
            for i in range(n):
                for j in range(i+1, n):
                    # Skip if atoms are bonded
                    if (i,j) not in bonded_pairs:
                        p1, p2 = mol.particles[i], mol.particles[j]
                        r = np.linalg.norm(p1.pos - p2.pos)
                        if r < self.Rc:
                            p1.neighbors.append(p2)
                            p2.neighbors.append(p1)

        # Then check between different molecules
        # Check molecule pairs
        n = len(self.molecules)
        for i in range(n):
            mol_i = self.molecules[i]
            for j in range(i+1, n):
                mol_j = self.molecules[j]
                # Check if molecules are close enough
                d = np.linalg.norm(mol_i.com - mol_j.com)
                if d < self.Rg + self.Rc:
                    # Check all particle pairs
                    for p1 in mol_i.particles:
                        for p2 in mol_j.particles:
                            r = np.linalg.norm(p1.pos - p2.pos)
                            if r < self.Rc:
                                p1.neighbors.append(p2)
                                p2.neighbors.append(p1)
        
    def _solve_bond_constraint(self, bond, p1, p2):
        """Solve single bond constraint"""
        r = p2.pos - p1.pos
        d = np.linalg.norm(r)
        if d > 0:
            dr = r * (1 - bond.r0/d)  # displacement vector
            w1 = 1/(p1.m + p2.m)      # weight factor
            dp1 = w1 * p1.m * dr      # displacement for p1
            dp2 = -w1 * p2.m * dr     # displacement for p2
            return dp1, dp2
        return np.zeros(2), np.zeros(2)

    def _solve_collision_constraint(self, p1, p2):
        """Solve single collision constraint"""
        r = p2.pos - p1.pos
        d = np.linalg.norm(r)
        if d < self.Rc - self.D and d > 0:  # Check if collision constraint applies
            dr = r * (1 - (self.Rc - self.D)/d)  # displacement vector
            w1 = 1/(p1.m + p2.m)      # weight factor
            dp1 = w1 * p1.m * dr      # displacement for p1
            dp2 = -w1 * p2.m * dr     # displacement for p2
            return dp1, dp2
        return np.zeros(2), np.zeros(2)

    def prepare_pd_system(self):
        """Prepare system for projective dynamics solver"""
        from sparse import build_neighbor_list, neigh_stiffness, make_Aii
        
        # Collect all particles, bonds, and their properties
        all_bonds = []    # collect all bonds
        all_ks = []       # collect all spring constants
        all_masses = []   # collect all masses
        points = []       # collect all positions
        velocities = []   # collect all velocities
        all_particles = []  # collect all particles for collision handling
        offset = 0
        for mol in self.molecules:
            for p in mol.particles:
                points.append(p.pos[0])  # x-coordinate only for now
                velocities.append(p.vel[0])
                all_masses.append(p.m)
                all_particles.append(p)
            for b in mol.bonds:
                all_bonds.append([b.i + offset, b.j + offset])
                all_ks.append(b.k)
            offset += len(mol.particles)
        
        points = np.array(points)
        velocities = np.array(velocities)
        all_masses = np.array(all_masses)
        
        # Build neighbor lists and stiffness matrices for bonds
        n_points = offset
        neighbs = build_neighbor_list(all_bonds, n_points)
        neighs, kngs, _ = neigh_stiffness(neighbs, all_bonds, all_ks)
        
        # Update collision neighbors and add them to the system
        self._update_neighbors()
        for i, p1 in enumerate(all_particles):
            for p2 in p1.neighbors:
                j = all_particles.index(p2)
                # Add collision neighbor if not already in bond neighbors
                if j not in neighs[i]:
                    neighs[i].append(j)
                    kngs[i].append(self.K_col)
        
        # Make diagonal terms with mass/dt^2 and all stiffness terms
        Aii0 = all_masses / (self.dt * self.dt)
        Aii = make_Aii(neighs, kngs, Aii0)  # Now includes both bond and collision stiffness
        
        # Predict positions using current velocities
        pos_pred = points + velocities * self.dt
        
        return points, pos_pred, velocities, neighs, kngs, Aii, all_masses, all_particles

    def make_pd_rhs(self, neighs, kngs, Aii, masses, points, pos_pred):
        """Build RHS according to projective dynamics, including collisions"""
        n_points = len(points)
        b = np.zeros(n_points)
        idt2 = 1.0 / (self.dt * self.dt)
        
        # Mass and constraint terms for all particles
        for i in range(n_points):
            # Mass term (inertial prediction)
            bi = pos_pred[i] * (masses[i] * idt2)
            
            # Both bond and collision terms (unified in neighs and kngs)
            for j, k in zip(neighs[i], kngs[i]):
                d = pos_pred[i] - pos_pred[j]
                d_norm = abs(d)  # 1D case
                if d_norm > 1e-10:  # Avoid division by zero
                    bi += k * d/d_norm  # Unit vector in 1D is just sign
            
            b[i] = bi
        
        return b

    def solve_pd_step(self, niter=100, tol=1e-8, callback=None):
        """Solve one step of projective dynamics"""
        from sparse import linsolve_iterative, jacobi_iteration_sparse
        
        # Prepare system
        points, pos_pred, velocities, neighs, kngs, Aii, masses, particles = self.prepare_pd_system()
        
        # Build RHS
        b = self.make_pd_rhs(neighs, kngs, Aii, masses, points, pos_pred)
        
        # Define Jacobi iteration function
        def update_func(A, b, x):
            return jacobi_iteration_sparse(x, b, neighs, kngs, Aii)
        
        # Solve system
        errs = []
        x_sol = linsolve_iterative(update_func, b, None, x0=points, niter=niter, tol=tol, errs=errs, callback=callback)
        
        # Update particle positions and velocities
        offset = 0
        for mol in self.molecules:
            for p in mol.particles:
                p.pos[0] = x_sol[offset]  # Update x-coordinate
                p.vel[0] = (x_sol[offset] - points[offset]) / self.dt  # Update velocity
                offset += 1
        
        return errs  # Return convergence history

    def _predict_positions(self):
        """Step 1: Predict positions using current velocities"""
        for mol in self.molecules:
            for p in mol.particles:
                p.pos_pred = p.pos + p.vel * self.dt
                p.pos = p.pos_pred.copy()  # Initialize position for constraint solving
                
    def _solve_constraints(self):
        """Step 2: Solve position constraints (bonds and collisions)"""
        self._update_neighbors()  # Update neighbor list before solving constraints
        
        # Iterate to solve constraints
        for _ in range(self.n_iter):
            # Solve bond constraints
            for mol in self.molecules:
                for bond in mol.bonds:
                    p1, p2 = mol.particles[bond.i], mol.particles[bond.j]
                    dp1, dp2 = self._solve_bond_constraint(bond, p1, p2)
                    p1.pos += dp1
                    p2.pos += dp2
            
            # Solve collision constraints
            for mol in self.molecules:
                for p1 in mol.particles:
                    for p2 in p1.neighbors:
                        dp1, dp2 = self._solve_collision_constraint(p1, p2)
                        p1.pos += dp1
                        p2.pos += dp2

    def _update_velocities(self):
        """Step 3: Update velocities from position changes"""
        for mol in self.molecules:
            for p in mol.particles:
                p.vel = (p.pos - p.pos_pred) / self.dt
                p.pos_pred = None  # cleanup
                
    def step(self):
        """Perform one time step of the simulation"""
        self._predict_positions()
        self._solve_constraints()
        self._update_velocities()

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
    errs = world.solve_pd_step(niter=100, tol=1e-8, callback=callback)
    
    # Plot convergence
    plt.figure(figsize=(10,6))
    plt.semilogy(errs, 'b-', label='Jacobi')
    plt.grid(True)
    plt.xlabel('Iteration')
    plt.ylabel('Error (log scale)')
    plt.title('Convergence of Jacobi Method for Water Molecules')
    plt.legend()
    plt.show()