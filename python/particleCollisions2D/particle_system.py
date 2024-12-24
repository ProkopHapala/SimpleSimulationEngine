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
    # Create a simple water molecule (H2O)
    p1 = Particle([0, 0], mass=16.0)      # O
    p2 = Particle([1, 0], mass=1.0)       # H1
    p3 = Particle([0, 1], mass=1.0)       # H2
    
    bonds = [
        Bond(0, 1, k=100.0, r0=1.0),  # O-H1 bond
        Bond(0, 2, k=100.0, r0=1.0),  # O-H2 bond
    ]
    
    water = Molecule([p1, p2, p3], bonds)
    
    # Create simulation world
    world = World(dt=0.01, D=0.1)
    world.add_molecule(water)
    
    # Run simulation for 100 steps
    for _ in range(100):
        world.step()
