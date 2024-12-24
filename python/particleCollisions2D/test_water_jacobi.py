import numpy as np
import matplotlib.pyplot as plt
from particle_system import Particle, Bond, Molecule, World
from sparse import linsolve_iterative, jacobi_iteration_sparse, build_neighbor_list, neigh_stiffness, make_Aii
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

def create_water_grid(n=3, spacing=3.0):
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
    world = World(dt=0.01, D=0.1)
    molecules = create_water_grid(n=3, spacing=3.0)
    for mol in molecules:
        world.add_molecule(mol)
    
    # Setup visualization
    vis = Visualizer(world)
    
    # Prepare data for Jacobi solver
    all_bonds = []    # collect all bonds
    all_ks = []       # collect all spring constants
    offset = 0
    for mol in world.molecules:
        for b in mol.bonds:
            all_bonds.append([b.i + offset, b.j + offset])
            all_ks.append(b.k)
        offset += len(mol.particles)
    
    # Convert to numpy arrays
    bonds = np.array(all_bonds)
    ks = np.array(all_ks)
    
    # Build neighbor lists and stiffness matrices
    n_points = offset
    neighbs = build_neighbor_list(bonds, n_points)
    neighs, kngs, _ = neigh_stiffness(neighbs, bonds, ks)
    Aii = make_Aii(neighs, kngs)
    
    # Setup initial state and target
    x0 = np.zeros(n_points)  # initial positions
    b = np.zeros(n_points)   # target positions
    
    # Set initial positions from current system state
    offset = 0
    for mol in world.molecules:
        for p in mol.particles:
            x0[offset] = p.pos[0]  # x-coordinate
            offset += 1
    
    # Set target positions (shift all molecules to the right)
    b = x0 + 2.0  # shift all x-coordinates by 2.0 units
    
    # Run Jacobi iteration and collect convergence data
    err_jacobi = []
    niter = 100
    
    # Setup visualization
    vis = Visualizer(world)
    
    def callback(itr, x, r):
        err = np.linalg.norm(r)
        print(f"Iteration {itr:4d} residual: {err:10.6f}")
        if itr % 5 == 0:  # Update visualization every 5 iterations
            # Update particle positions based on current solution
            offset = 0
            for mol in world.molecules:
                for p in mol.particles:
                    p.pos[0] = x[offset]  # Update x-coordinate
                    offset += 1
            vis.update()
    
    def update_func(A, b, x):
        return jacobi_iteration_sparse(x, b, neighs, kngs, Aii)
    
    print("\nStarting Jacobi iteration:")
    x_sol = linsolve_iterative(update_func, b, None, x0=x0, niter=niter, tol=1e-8, errs=err_jacobi, callback=callback)
    
    # Plot convergence
    plt.figure(figsize=(10,6))
    plt.semilogy(err_jacobi, 'b-', label='Jacobi')
    plt.grid(True)
    plt.xlabel('Iteration')
    plt.ylabel('Error (log scale)')
    plt.title('Convergence of Jacobi Method for Water Molecules')
    plt.legend()
    plt.show()
    
    # Run simulation with visualization
    n_steps = 1000
    for step in range(n_steps):
        world.step()
        if step % 5 == 0:
            vis.update()
        if step % 100 == 0:
            print(f"Step {step}")
            for i, mol in enumerate(world.molecules):
                ke = sum(0.5 * p.m * np.dot(p.vel, p.vel) for p in mol.particles)
                print(f"  Molecule {i} kinetic energy: {ke:.3f}")
