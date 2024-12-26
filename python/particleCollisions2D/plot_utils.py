import matplotlib.pyplot as plt
import numpy as np

def plot_molecules(world, x=None,y=None, ax=None, color='b', alpha=1.0, label=None, marker='o', markersize=5):
    """Plot molecules in the world"""
    if ax is None:
        ax = plt.gca()
    
    # Plot particles
    if x is None: x = world.apos[:,0]
    if y is None: y = world.apos[:,1]
    ax.scatter(x, y,  alpha=alpha, color=color, label=label, marker=marker, s=markersize)
    
    # Plot bonds
    for mol in world.molecules:
        for bond in mol.bonds:
            i, j = bond.i, bond.j
            ax.plot([x[i], x[j]], [y[i], y[j]], color=color, alpha=alpha)
    
    ax.set_aspect('equal')
    return ax

def plot_convergence(errors, ax=None):
    """Plot convergence of the solver"""
    if ax is None:
        ax = plt.gca()
    
    ax.plot(errors, '-o')
    ax.set_yscale('log')
    ax.set_xlabel('Iteration')
    ax.set_ylabel('Error')
    ax.grid(True)
    return ax
