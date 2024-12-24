from particle_system import Particle, Bond, Molecule, World
from visualizer import Visualizer
import numpy as np

def create_water():
    """Create a water molecule (H2O)"""
    p1 = Particle([0, 0], vel=[0.5, 0.0], mass=16.0)      # O
    p2 = Particle([1, 0], vel=[0.5, 0.2], mass=1.0)       # H1
    p3 = Particle([0, 1], vel=[0.5,-0.2], mass=1.0)       # H2
    
    bonds = [
        Bond(0, 1, k=100.0, r0=1.0),  # O-H1 bond
        Bond(0, 2, k=100.0, r0=1.0),  # O-H2 bond
    ]
    
    return Molecule([p1, p2, p3], bonds)

def create_co2():
    """Create a CO2 molecule"""
    p1 = Particle([-1, 0], vel=[-0.5, 0.0], mass=12.0)    # C
    p2 = Particle([-2, 0], vel=[-0.5, 0.1], mass=16.0)    # O1
    p3 = Particle([ 0, 0], vel=[-0.5,-0.1], mass=16.0)    # O2
    
    bonds = [
        Bond(0, 1, k=100.0, r0=1.0),  # C-O1 bond
        Bond(0, 2, k=100.0, r0=1.0),  # C-O2 bond
    ]
    
    return Molecule([p1, p2, p3], bonds)

def create_triangle():
    """Create a triangle molecule with only two bonds (third side unbonded)"""
    p1 = Particle([3, 0], vel=[0.0, 0.0], mass=1.0)
    p2 = Particle([4, 0], vel=[0.0, 0.2], mass=1.0)
    p3 = Particle([3.5, 0.866], vel=[0.0,-0.2], mass=1.0)  # height = sqrt(3)/2 * side
    
    bonds = [
        Bond(0, 1, k=100.0, r0=1.0),  # base
        Bond(1, 2, k=100.0, r0=1.0),  # right side
        # Left side (0-2) is not bonded, should be maintained by collisions
    ]
    
    return Molecule([p1, p2, p3], bonds)

if __name__ == "__main__":
    # Create simulation world
    world = World(dt=0.01, D=0.1)
    
    # Add molecules
    water = create_water()
    co2 = create_co2()
    triangle = create_triangle()
    world.add_molecule(water)
    world.add_molecule(co2)
    world.add_molecule(triangle)
    
    # Setup visualization
    vis = Visualizer(world)
    
    # Run simulation
    n_steps = 1000
    for step in range(n_steps):
        world.step()
        if step % 5 == 0:  # Update visualization every 5 steps
            vis.update()
            
        # Print some debug info occasionally
        if step % 100 == 0:
            print(f"Step {step}")
            for i, mol in enumerate(world.molecules):
                ke = sum(0.5 * p.m * np.dot(p.vel, p.vel) for p in mol.particles)
                print(f"  Molecule {i} kinetic energy: {ke:.3f}")
