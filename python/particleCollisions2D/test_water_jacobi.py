import numpy as np
import matplotlib.pyplot as plt
from particle_system import Particle, Bond, Molecule, World
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

if __name__ == "__main__":
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
