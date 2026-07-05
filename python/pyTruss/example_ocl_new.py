#!/usr/bin/env python3
# === AUTO-DOC BEGIN ===
"""
@brief Example script demonstrating the new TrussSolverOCL API with GPU solver switching.

Shows how to construct a Truss, create a TrussSolverOCL instance, and switch between
GPU linear solvers (VBD, Jacobi-diff, Jacobi-fly) via the `get_solver` factory. Runs
a short simulation and prints max displacement. Serves as the minimal usage example
for the refactored GPU backend in `truss_solver_ocl_new.py`.
"""
# === AUTO-DOC END ===

import numpy as np
from truss import Truss
from truss_solver_ocl_new import TrussSolverOCL, get_solver

# Create a simple test truss (hanging bridge)
def make_bridge_truss(n=10, L=1.0, k=100.0):
    points = np.zeros((n, 3))
    for i in range(n): points[i] = [i * L, 0.0, 0.0]
    bonds, masses = [], []
    for i in range(n-1): bonds.append([i, i+1])
    masses = np.ones(n) * 0.1
    truss = Truss()
    truss.points  = points.astype(np.float64)
    truss.bonds   = bonds
    truss.ks      = np.full(len(bonds), k, dtype=np.float64)
    truss.masses  = masses.astype(np.float64)
    truss.fixed   = {0, n-1}  # Fix both ends
    return truss

if __name__ == "__main__":
    # Create truss
    truss = make_bridge_truss(n=20, L=0.5, k=500.0)
    print(f"Truss: {len(truss.points)} points, {len(truss.bonds)} bonds")
    
    # Initialize solver
    dt = 0.01
    gravity = np.array([0.0, 0.0, -9.8])
    solver = TrussSolverOCL(truss, dt=dt, gravity=gravity, nloc=32)
    
    # Test different solvers
    solvers_to_test = [
        ('vbd_serial',   {'niter': 10, 'det_eps': 1e-6}),
        ('jacobi_fly',   {'niter': 20}),
        ('jacobi_diff',  {'niter': 20, 'use_double': True}),
    ]
    
    for solver_name, config in solvers_to_test:
        print(f"\n{'='*60}")
        print(f"Testing solver: {solver_name}")
        print(f"Config: {config}")
        print(f"{'='*60}")
        
        # Reset to initial state
        solver.reset_state(positions=truss.points.copy(), velocities=None)
        
        # Get solver callback
        solver_callback = get_solver(solver_name)
        
        # Run simulation
        nsteps = 100
        track_indices = [len(truss.points) // 2]  # Track middle point
        
        x_final, v_final, trajectory = solver.run(
            nsteps=nsteps,
            solver_callback=solver_callback,
            solver_config=config,
            track_indices=track_indices,
            verbose=True
        )
        
        # Report results
        mid_point = len(truss.points) // 2
        print(f"\nFinal position of middle point: {x_final[mid_point]}")
        print(f"Final vertical displacement: {x_final[mid_point, 2] - truss.points[mid_point, 2]:.4f} m")
        
        if trajectory is not None:
            print(f"Trajectory shape: {trajectory.shape}")
            print(f"Max vertical displacement: {(trajectory[:, 0, 2] - truss.points[mid_point, 2]).min():.4f} m")
    
    print(f"\n{'='*60}")
    print("All tests completed successfully!")
    print(f"{'='*60}")
