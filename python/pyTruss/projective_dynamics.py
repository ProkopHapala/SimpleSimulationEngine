import numpy as np
#from sparse import build_neighbor_list
import sparse

def make_pd_Aii0(  masses, dt ):
    """Create the diagonal terms of the system matrix for projective dynamics"""
    idt2 = 1.0 / (dt * dt)
    return masses * idt2

def make_pd_matrix(neighbs, bonds, masses, dt, ks):
    """Create the system matrix for projective dynamics"""
    np_total = len(masses)
    A    = np.zeros((np_total, np_total))
    Aii0 = make_pd_Aii0( masses, dt )
    for i in range(np_total):
        Aii = Aii0[i]
        for ib in neighbs[i]:
            k = ks[ib]
            i_, j_ = bonds[ib]
            Aii += k
            if j_ > i:
                A[i, j_] = -k
                A[j_, i] = -k
            elif i_ > i:
                A[i, i_] = -k
                A[i_, i] = -k
        A[i, i] = Aii
    return A

def make_pd_rhs(neighbs, bonds, masses, dt, ks, points, l0s, pnew):
    """Build the right-hand side of the system following the Julia implementation"""
    np_total = len(masses)
    b = np.zeros((np_total, 3))
    idt2 = 1.0 / (dt * dt)
    
    for i in range(np_total):
        # Mass term (inertial prediction)
        bi = pnew[i] * (masses[i] * idt2)
        
        # Spring terms
        for ib in neighbs[i]:
            k = ks[ib]
            i_, j_ = bonds[ib]
            j = j_ if i_ == i else i_
            
            # Using predicted positions for better propagation
            d = pnew[i] - pnew[j]
            d_norm = np.linalg.norm(d)
            if d_norm > 1e-10:  # Avoid division by zero
                d *= k * l0s[ib] / d_norm
                bi += d
        
        b[i] = bi
    
    return b

def update_velocity(ps_cor, points, velocity, dt):
    """Update velocity based on position change"""
    return (ps_cor - points) / dt

def solve_pd(points, velocity, bonds, masses, ks, dt=0.1, n_iter=100, gravity=np.array([0, -9.81, 0]), fixed_points=None, call_back=None, damping=0.01 ):
    """Solve the system using projective dynamics following the Julia implementation"""
    n_points = len(points)
    neighbs = sparse.build_neighbor_list(bonds, n_points)
    A = make_pd_matrix(neighbs, bonds, masses, dt, ks)
    
    # Initialize
    pos      = points.copy()
    pos_cor  = points.copy()
    pos_pred = points.copy()
    l0s = np.array([np.linalg.norm(points[j] - points[i]) for i, j in bonds])
    
    # Main simulation loop
    for itr in range(n_iter):
        
        velocity*=(1-damping)
        velocity[:,:] += gravity[None,:] * dt
        pos_pred[:] = pos + velocity * dt
        
        # Apply fixed point constraints
        if fixed_points is not None:
            pos_pred[fixed_points] = points[fixed_points]
        
        # Build right-hand side
        b = make_pd_rhs(neighbs, bonds, masses, dt, ks, pos, l0s, pos_pred)
        
        # Solve system for each coordinate
        for i in range(3):
            pos_cor[:, i] = np.linalg.solve(A, b[:, i])
        
        # Apply fixed point constraints
        if fixed_points is not None:
            pos_cor[fixed_points] = points[fixed_points]
        
        # Update velocity and position
        velocity = update_velocity(pos_cor, pos, velocity, dt)
        pos[:]   = pos_cor[:]
        
        if call_back is not None:
            call_back(pos)
        
        # Print debug info
        v_norm   = np.linalg.norm(velocity)
        pos_diff = np.linalg.norm(pos - pos_pred)
        print(f"iter:{itr} dt={dt:.3e} |v|={v_norm:.3e} dp={pos_diff:.3e}")
        
        # Optional: check bond lengths for debugging
        if False:  # Set to True to enable
            bond_errors = []
            for (i, j), l0 in zip(bonds, l0s):
                current_length = np.linalg.norm(pos[j] - pos[i])
                bond_errors.append(abs(current_length - l0))
            max_bond_error = max(bond_errors)
            print(f"    max bond length error: {max_bond_error:.3e}")
    
    return pos, velocity

# def makeSparseSystem( dt, bonds, points, masses, ks, fixed, l0s, neighbs, gravity = np.array([0., -9.81, 0.0 ]) ):
#     neighs, kngs, n_max = sparse.neigh_stiffness(neighbs, bonds, ks)
#     Aii0     = make_pd_Aii0(  masses, dt )
#     Aii      = sparse.make_Aii(neighs, kngs, Aii0)
#     velocity = points*0 + gravity[None,:] * dt
#     pos_pred = points + velocity * dt
#     b        = make_pd_rhs(neighbs, bonds, masses, dt, ks, points, l0s, pos_pred)
#     return neighs, kngs, Aii, b, pos_pred, velocity, n_max

def makeSparseSystem( dt, bonds, points, masses, ks, fixed, l0s, neighbs, gravity = np.array([0., -9.81, 0.0 ]) ):
    velocity = points*0 + gravity[None,:] * dt
    pos_pred = points + velocity * dt
    if fixed is not None:
        pos_pred[fixed] = points[fixed]
    b = make_pd_rhs(neighbs, bonds, masses, dt, ks, points, l0s, pos_pred)
    neighs, kngs, n_max = sparse.neigh_stiffness(neighbs, bonds, ks)
    Aii0 = make_pd_Aii0(  masses, dt )
    Aii = sparse.make_Aii(neighs, kngs, Aii0)
    return neighs, kngs, Aii, b, pos_pred, velocity, n_max

def apply_pd_matrix(x, neighbs, bonds, masses, dt, ks):
    """Apply system matrix A to vector x without explicitly building A"""
    n = len(masses)
    result = np.zeros_like(x)
    idt2 = 1.0 / (dt * dt)
    
    for i in range(n):
        # Diagonal term (mass + spring coefficients)
        Aii = masses[i] * idt2
        xi = x[i]
        
        # Accumulate spring contributions
        for ib in neighbs[i]:
            k = ks[ib]
            i_, j_ = bonds[ib]
            j = j_ if i_ == i else i_
            
            Aii += k  # Add to diagonal
            result[i] += k * (xi - x[j])  # Off-diagonal contribution
            
        result[i] += (masses[i] * idt2) * xi  # Mass term
    
    return result

# Example usage
if __name__ == "__main__":

    #import numpy as np
    import matplotlib.pyplot as plt
    from truss import Truss
    import plot_utils as pu 

    bWheel = False

    # Create a truss system
    truss = Truss()
    if bWheel:
        p0 = np.array([0., 0., 0.])
        p1 = np.array([1., 0., 0.])
        ax = np.array([0., 0., 1.])
        k_dict = [10000.0, 5000.0, 2000.0, 2000.0]  # [long, perp, zigIn, zigOut]
        truss.wheel(p0, p1, ax, width=0.2, n=8, k_scale=1.0, k_dict=k_dict)
    else:
        truss.build_grid_2d(nx=5, ny=5, m=1.0, m_end=1000.0, l=1.0, k=10000.0, k_diag=1000.0)

    # Get quantities needed for projective dynamics
    bonds, points, masses, ks, fixed, l0s, neighbs = truss.get_pd_quantities()
    velocity = np.zeros_like(points)

    # Visualize results
    ax = pu.plot_truss(points, truss.bonds, edge_color='b', label='Initial Points')

    # Solve using projective dynamics
    new_points, new_velocity = solve_pd( points, velocity, bonds, masses, ks,  dt=0.02, n_iter=20,  fixed_points=fixed,  
                                         call_back=lambda x: pu.plot_truss(x, truss.bonds, ax=ax, edge_color='k', edge_alpha=0.1) )

    pu.plot_truss(new_points, truss.bonds, ax=ax, edge_color='r', label='Final Points')
    plt.show()
