import numpy as np

def jacobi_iteration(A, b, x):
    Aii = np.diag(A)
    r   = b - A @ x
    return ( r + Aii*x ) / Aii, r

def linsolve_Jacobi( b, A, x0=None, niter=10, tol=1e-6, bPrint=False, callback=None ):
    if x0 is None: x0 = np.zeros_like(b)
    x = x0.copy()
    for itr in range(niter):
        x, r = jacobi_iteration(A, b, x)
        err = np.linalg.norm(r)
        if callback is not None:
            callback(itr, x, r)
        if bPrint:
            print(f"linsolve_Jacobi() itr: {itr}, err: {err}")
        if err < tol:
            break
    return x

def linsolve_iterative( update_func, b, A, x0=None, niter=10, tol=1e-6, bPrint=False, callback=None, errs=None, bmix=1.0, niter_mix=100000 ):
    if x0 is None: x0 = np.zeros_like(b)
    x = x0.copy()
    d = x*0
    err = 1 
    for itr in range(niter):
        x_, r = update_func(A, b, x)

        err_prev = err
        err = np.linalg.norm(r)
        if itr>niter_mix:
            # if err > err_prev:
            #     bmix = min(1.0, bmix * 1.1)  # Increase mixing if error decreasing
            # else:
            #     bmix = max(0.7, bmix * 0.97)  # Decrease mixing if error increasing

            #x = bmix*x_ + (1-bmix)*x
            d_ = x_ - x
            x = x + bmix*d_ + (1-bmix)*d
            d = d_
        else:
            x = x_

        
        if errs is not None:
            errs.append(err)
        if callback is not None:
            callback(itr, x, r)
        if bPrint:
            print(f"linsolve_iterative() itr: {itr}, err: {err} bmix: {bmix}")
        if err < tol:
            break
    return x

def build_neighbor_list(bonds, n_points):
    """Build list of neighboring bonds for each point"""
    neighbs = [[] for _ in range(n_points)]
    for i, (i_, j_) in enumerate(bonds):
        neighbs[i_].append(i)
        neighbs[j_].append(i)
    return neighbs

def neigh_stiffness(neighbs, bonds, kbs):
    '''
    convert bond-centerd neighbor-list and stiffness to point-centered neighbor-list and stiffness
    '''
    n = len(neighbs)
    neighs = []
    kngs   = []
    ni_max = 0
    for i in range(n):
        ngi = []
        ngki = []
        ni_max = max(ni_max, len(neighbs[i]))
        for ib in neighbs[i]:
            k = kbs[ib]
            i_, j_ = bonds[ib]
            j = j_ if i_ == i else i_
            ngi .append( j )
            ngki.append( k )
        #print("neigh_stiffness: i: ", i, " ngi: ", len(ngi), " ngki: ", len(ngki) )
        neighs.append( ngi )
        kngs.append( ngki )
    return neighs, kngs, ni_max

def neighs_to_dense_arrays(neighs, kngs, n_max):
    '''
    create dense arrays from point-centered neighbor-list and stiffness which can be efficiently used in C or OpenCL
    '''
    n = len(neighs)
    nis      = np.zeros(n, dtype=np.int32)
    neighs_  = np.zeros((n, n_max), dtype=np.int32  )
    kngs_    = np.zeros((n, n_max), dtype=np.float32 )
    for i in range(n):
        ngi  = neighs[i]
        ni   = len(ngi)
        nis[i] = ni
        neighs_[i,:ni] = ngi
        kngs_  [i,:ni] = kngs[i]
    return neighs_, kngs_, nis

def make_Aii( neighs, kngs, Aii0=None ):
    n = len(neighs)
    Aii = np.zeros(n)
    for i in range(n):
        if Aii0 is None: 
            aii = 0.0
        else:
            aii = Aii0[i] # helps diagonal dominance, use mass[i] / dt^2 for Projective Dynamics
        ngsi = neighs[i]
        ksi  = kngs[i] 
        ni   = len(ngsi)
        for jj in range(ni):
            aii   += ksi[jj]
        Aii[i] = aii
    return Aii

def dot_sparse(x, neighs, kngs, Aii ):
    """Sparse version of y = A @ x """
    n = len(x)
    y = np.zeros_like(x)
    for i in range(n):
        ngsi = neighs[i]
        ksi  = kngs[i] 
        ni   = len(ngsi)
        sum_j  = 0
        for jj in range(ni):
            j      = ngsi[jj]
            k      = ksi[jj]
            sum_j -= k * x[j]   # Off-diagonal contribution
        y[i] = sum_j + Aii[i] * x[i]
    return y

def jacobi_iteration_sparse(x, b, neighs, kngs, Aii ):
    """One iteration of Jacobi method using sparse operations"""
    n = len(x)
    x_out = np.zeros_like(x)
    r     = np.zeros_like(x)
    #print("n = ", n, " len(x) = ", len(x))
    for i in range(n):
        sum_j  = 0  # RHS term
        ngsi = neighs[i]
        ksi  = kngs[i] 
        ni = len(ngsi)
        for jj in range(ni):
            j      = ngsi[jj]
            k      = ksi[jj]
            #sum_j -= k * x[j]   # Off-diagonal contribution
            sum_j += k * x[j]   # Off-diagonal contribution
        x_out[i] =  (b[i] + sum_j) / Aii[i]   # solution x_new = (b - sum_(j!=i){ Aij * x[j] } ) / Aii
        r[i]     = b[i] + sum_j - Aii[i]*x[i] # Residual r = b - Ax ;  Ax = Aii * x[i] + sum_(j!=i){ Aij * x[j] }

        #x_out[i] =  (b[i] - sum_j) / Aii[i]                 # solution x_new = (b - sum_(j!=i){ Aij * x[j] } ) / Aii
        #r[i]     = b[i] - ( sum_j + Aii[i] * x[i]) # Residual r = b - Ax ;  Ax = Aii * x[i] + sum_(j!=i){ Aij * x[j] }
    
    #r = b - dot_sparse(x, neighs, kngs, Aii )
    return x_out, r

# def jacobi_iteration_sparse(x, b, neighs, kngs, Aii ):
#     """One iteration of Jacobi method using sparse operations"""
#     n = len(x)
#     x_out = np.zeros_like(x)
#     r     = np.zeros_like(x)
#     for i in range(n):
#         sum_j  = b[i]  # Start with RHS
#         ngsi = neighs[i]
#         ksi  = kngs[i] 
#         ni = len(ngsi)
#         for jj in range(ni):
#             j = ngsi[jj]
#             k = ksi[jj]
#             #if i != j:  # Skip diagonal terms
#             sum_j += k * x[j]   # Add because A[i,j] is -k, so when moved to RHS it becomes +k
#         x_out[i] = sum_j / Aii[i]  # x_new = (b + sum_(j!=i){ k * x[j] }) / Aii
#     r = b - dot_sparse(x, neighs, kngs, Aii )  
#     r = b - dot_sparse(x, neighs, kngs, Aii )  
#     return x_out, r


def linsolve_Jacobi_sparse( b, neighs, kngs, x0=None, Aii0=None, niter=10, tol=1e-6, bPrint=False, callback=None ):
    if x0 is None: x = np.zeros_like(b)
    x = x0.copy()
    Aii = make_Aii( neighs, kngs, Aii0 )
    for itr in range(niter):
        x, r = jacobi_iteration_sparse(x, b, neighs, kngs, Aii )
        err = np.linalg.norm(r)
        if callback is not None:
            callback(itr, x, r)
        if bPrint:
            print(f"linsolve_Jacobi_sparse() itr: {itr}, err: {err}")
        if err < tol:
            break
    return x

# check diagonal vs ref
def check_diagonal(A, Aii, bPrint=True):
    Aii_dense = np.diag(A)
    if bPrint:
        print("check_diagonal()\n#i Sparse_Aii vs Dense_Aii:")
        for i in range(len(Aii)):
            print(f"{i}: {Aii[i]:.6f} vs {Aii_dense[i]:.6f}")
    error = np.linalg.norm(Aii - Aii_dense)
    print(f"check_diagonal() Total error: {error:.6f}")

if __name__ == "__main__":
    from projective_dynamics import make_pd_matrix, make_pd_rhs, make_pd_Aii0
    from truss               import Truss

    bWheel = False
    dt = 1.0
    #dt = 0.5
    #dt = 0.2
    gravity = np.array([0., -9.81, 0.0 ])

    fixed_points = [0]

    truss = Truss()
    if bWheel:
        p0 = np.array([0., 0., 0.])
        p1 = np.array([1., 0., 0.])
        ax = np.array([0., 0., 1.])
        k_dict = [10000.0, 5000.0, 2000.0, 2000.0]  # [long, perp, zigIn, zigOut]
        truss.wheel(p0, p1, ax, width=0.2, n=8, k_scale=1.0, k_dict=k_dict)
    else:
        truss.build_grid_2d(nx=5, ny=5, m=1.0, m_end=1000.0, l=1.0, k=10000.0, k_diag=1000.0)

    bonds, points, masses, ks, fixed, l0s, neighbs = truss.get_pd_quantities()
    neighbs = build_neighbor_list(bonds, len(points))
    A = make_pd_matrix(neighbs, bonds, masses, dt, ks)
    velocity = points*0 + gravity[None,:] * dt
    pos_pred = points + velocity * dt
    if fixed_points is not None:
        pos_pred[fixed_points] = points[fixed_points]
    b = make_pd_rhs(neighbs, bonds, masses, dt, ks, points, l0s, pos_pred)

    neighs, kngs, n_max = neigh_stiffness(neighbs, bonds, ks)
    print("n_max = ",  n_max)
    #print("neighs = ", neighs)
    #print("kngs = ",   kngs)

    Aii0 = make_pd_Aii0(  masses, dt )
    bx = b[:,1]
    x0 = pos_pred[:,1]


    # check diagonal
    Aii_sparse = make_Aii(neighs, kngs, Aii0)
    check_diagonal(A, Aii_sparse, bPrint=False)

    # check Matrix-vector product
    Ax_ref    = A@x0
    Ax_sparse = dot_sparse(x0, neighs, kngs, Aii_sparse)
    #print("Ax_ref = ", Ax_ref)
    #print("Ax_sparse = ", Ax_sparse)
    print("| Ax_sparse - Ax_ref | = ", np.linalg.norm(Ax_sparse - Ax_ref))  

    #exit()

    x_ref = np.linalg.solve(A, bx)



    #x_jacobi = linsolve_Jacobi       ( bx, A,            x0=x0,            bPrint=True )
    #x_sparse = linsolve_Jacobi_sparse( bx, neighs, kngs, x0=x0, Aii0=Aii0, bPrint=True )

    niter = 50
    err_jacobi = []
    err_JacobMix = []
    err_sparse = []
    x_jacobi   = linsolve_iterative( lambda A, b, x: jacobi_iteration(A, b, x),                                bx, A,    x0=x0, niter=niter, tol=1e-8, errs=err_jacobi )
    x_JacobMix = linsolve_iterative( lambda A, b, x: jacobi_iteration(A, b, x),                                bx, A,    x0=x0, niter=niter, tol=1e-8, errs=err_JacobMix, niter_mix=2, bmix=0.75, bPrint=True )
    x_sparse = linsolve_iterative( lambda A, b, x: jacobi_iteration_sparse(x, b, neighs, kngs, Aii_sparse ), bx, None, x0=x0, niter=niter, tol=1e-8, bPrint=True, errs=err_sparse )

    #print("x_ref = ", x_ref)
    print("| x_jacobi - x_ref | = ", np.linalg.norm(x_jacobi - x_ref))
    print("| x_JacobMix - x_ref | = ", np.linalg.norm(x_JacobMix - x_ref))
    print("| x_sparse - x_ref | = ", np.linalg.norm(x_sparse - x_ref))


    import matplotlib.pyplot as plt
    plt.plot(err_jacobi, label="Jacobi")
    plt.plot(err_JacobMix, label="JacobMix")
    plt.yscale('log')
    plt.legend()
    plt.grid()
    plt.show()




