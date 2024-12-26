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

def to_1d_coords(points_2d):
    """Convert 2D points array to 1D array [x1,y1,x2,y2,...]"""
    return points_2d.flatten()

def to_2d_coords(points_1d):
    """Convert 1D points array [x1,y1,x2,y2,...] to 2D array [[x1,y1],[x2,y2],...]"""
    return points_1d.reshape(-1, 2)

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
        points = []       # collect all positions (x,y)
        velocities = []   # collect all velocities (x,y)
        all_particles = []  # collect all particles for collision handling
        offset = 0
        for mol in self.molecules:
            for p in mol.particles:
                points.extend([p.pos[0], p.pos[1]])  # both x and y coordinates
                velocities.extend([p.vel[0], p.vel[1]])
                all_masses.extend([p.m, p.m])  # mass repeated for x and y
                all_particles.append(p)
            for b in mol.bonds:
                # Each bond creates two constraints (x and y)
                all_bonds.extend([[2*(b.i + offset), 2*(b.j + offset)],     # x coordinates
                                [2*(b.i + offset)+1, 2*(b.j + offset)+1]])  # y coordinates
                all_ks.extend([b.k, b.k])  # spring constant for both x and y
            offset += len(mol.particles)
        
        points = np.array(points)
        velocities = np.array(velocities)
        all_masses = np.array(all_masses)
        
        # Build neighbor lists and stiffness matrices for bonds
        n_points = len(points)
        neighbs = build_neighbor_list(all_bonds, n_points)
        neighs, kngs, _ = neigh_stiffness(neighbs, all_bonds, all_ks)
        
        # Update collision neighbors and add them to the system
        self._update_neighbors()
        for i, p1 in enumerate(all_particles):
            for p2 in p1.neighbors:
                j = all_particles.index(p2)
                # Add collision neighbors for both x and y coordinates
                ix, iy = 2*i, 2*i+1
                jx, jy = 2*j, 2*j+1
                # Add x coordinate neighbor if not already present
                if jx not in neighs[ix]:
                    neighs[ix].append(jx)
                    kngs[ix].append(self.K_col)
                # Add y coordinate neighbor if not already present
                if jy not in neighs[iy]:
                    neighs[iy].append(jy)
                    kngs[iy].append(self.K_col)
        
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
        for i in range(0, n_points, 2):  # Process x,y coordinates together
            # Mass terms (inertial prediction) for x,y
            bi_x = pos_pred[i] * (masses[i] * idt2)
            bi_y = pos_pred[i+1] * (masses[i+1] * idt2)
            
            
            # Process neighbors for x coordinate
            for j, k in zip(neighs[i], kngs[i]):
                dx = pos_pred[i] - pos_pred[j]
                dy = pos_pred[i+1] - pos_pred[j+1]
                d = np.sqrt(dx*dx + dy*dy)  # true distance in 2D
                if d > 1e-10:  # Avoid division by zero
                    bi_x += k * dx/d  # x component of force
            
            # Process neighbors for y coordinate
            for j, k in zip(neighs[i+1], kngs[i+1]):
                dx = pos_pred[i] - pos_pred[j-1]
                dy = pos_pred[i+1] - pos_pred[j]
                d = np.sqrt(dx*dx + dy*dy)  # true distance in 2D
                if d > 1e-10:  # Avoid division by zero
                    bi_y += k * dy/d  # y component of force
            
            b[i] = bi_x
            b[i+1] = bi_y
        
        return b

    def prepare_constraints(self, points_1d):
        """Prepare combined constraints from bonds and collisions"""
        from sparse import build_neighbor_list, neigh_stiffness
        
        points_2d = to_2d_coords(points_1d)
        n_particles = len(points_2d)
        
        # Get static bond constraints
        bond_bonds = []
        bond_ks = []
        offset = 0
        for mol in self.molecules:
            for b in mol.bonds:
                # Each bond creates two constraints (x and y)
                bond_bonds.extend([[2*(b.i + offset), 2*(b.j + offset)],      # x coordinates
                                 [2*(b.i + offset)+1, 2*(b.j + offset)+1]])   # y coordinates
                bond_ks.extend([b.k, b.k])  # spring constant for both x and y
            offset += len(mol.particles)
        
        # Get collision constraints based on current positions
        col_bonds = []
        col_ks = []
        all_particles = [p for mol in self.molecules for p in mol.particles]
        
        # Update collision neighbors using current positions
        self._update_neighbors()
        for i, p1 in enumerate(all_particles):
            for p2 in p1.neighbors:
                j = all_particles.index(p2)
                # Add collision bonds for both x and y coordinates
                col_bonds.extend([[2*i, 2*j], [2*i+1, 2*j+1]])
                col_ks.extend([self.K_col, self.K_col])
        
        # Combine all constraints
        all_bonds = bond_bonds + col_bonds
        all_ks = bond_ks + col_ks
        
        # Build neighbor lists and stiffness matrices
        n_points = len(points_1d)
        neighbs = build_neighbor_list(all_bonds, n_points)
        neighs, kngs, _ = neigh_stiffness(neighbs, all_bonds, all_ks)
        
        return neighs, kngs

    def solve_pd_step(self, niter_outer=10, niter_inner=10, tol=1e-8, callback=None):
        """Solve one step of projective dynamics with periodic collision updates"""
        from sparse import make_Aii, jacobi_iteration_sparse
        
        # Collect initial state
        points_2d = np.array([p.pos for mol in self.molecules for p in mol.particles])
        vels_2d = np.array([p.vel for mol in self.molecules for p in mol.particles])
        masses = np.array([p.m for mol in self.molecules for p in mol.particles])
        
        # Convert to 1D arrays
        points = to_1d_coords(points_2d)
        velocities = to_1d_coords(vels_2d)
        masses_1d = np.repeat(masses, 2)  # repeat for x,y components
        
        # Predict positions
        pos_pred = points + velocities * self.dt
        x = points.copy()
        errs = []
        
        # Outer iteration loop - update neighbors periodically
        for it_outer in range(niter_outer):
            # Prepare constraints with current positions
            neighs, kngs = self.prepare_constraints(x)
            
            # Make diagonal terms
            Aii0 = masses_1d / (self.dt * self.dt)
            Aii = make_Aii(neighs, kngs, Aii0)
            
            # Inner iteration loop - solve with fixed neighbors
            for it_inner in range(niter_inner):
                # Build RHS (b vector)
                b = self.make_pd_rhs(neighs, kngs, Aii, masses_1d, x, pos_pred)
                
                # One step of Jacobi iteration
                x_new, r = jacobi_iteration_sparse(x, b, neighs, kngs, Aii)
                
                # Check convergence
                err = np.linalg.norm(r)
                errs.append(err)
                if callback: callback(it_outer*niter_inner + it_inner, x_new, r)
                if err < tol: break
                x = x_new
        
        # Update particle states
        points_2d_new = to_2d_coords(x)
        for i, p in enumerate([p for mol in self.molecules for p in mol.particles]):
            p.pos = points_2d_new[i]
            p.vel = (points_2d_new[i] - points_2d[i]) / self.dt
        
        return errs

    def _get_bond_neighbors(self):
        """Get static bond topology neighbors"""
        all_bonds = []
        all_ks = []
        offset = 0
        for mol in self.molecules:
            for b in mol.bonds:
                # Each bond creates one constraint between points
                all_bonds.append([b.i + offset, b.j + offset])
                all_ks.append(b.k)
            offset += len(mol.particles)
        return all_bonds, all_ks

    def _get_collision_neighbors(self, pos_pred):
        """Get current collision neighbors based on predicted positions"""
        collision_bonds = []
        collision_ks = []
        
        # Check all pairs of points for collisions
        n = len(pos_pred) // 2  # number of particles (each has x,y)
        for i in range(n):
            for j in range(i+1, n):
                # Get positions
                i2 = 2*i
                j2 = 2*j
                dx = pos_pred[i2] - pos_pred[j2]
                dy = pos_pred[i2+1] - pos_pred[j2+1]
                r = np.sqrt(dx*dx + dy*dy)
                
                # Add collision constraint if particles are close enough
                if r < self.Rc - self.D:
                    # Add one constraint between points
                    collision_bonds.append([i, j])
                    collision_ks.append(self.K_col)
        
        return collision_bonds, collision_ks

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