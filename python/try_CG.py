#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse.linalg as spla

# =========== Functions
itr=0

def stepCG( K, x, d, f, f2old ):
    global itr
    #print( "stepCG ", itr, "=========================="  )
    Kd  = np.dot(K,d)
    dt  = np.dot(d,f) / np.dot(d,Kd)    # step length such that f is orthogonal to d 
    #dt  = np.dot(f,f) / np.dot(d,Kd)
    #print( "d", d )
    #print( "ox", x )
    #print( "of", f )
    #print( "dt", dt )
    x   = x + d *dt
    f   = f - Kd*dt
    print( "x:", x, " f:", f )
    #print( "f", f )
    f2  = np.dot(f,f)
    #print( itr, "dt ", dt, "f2", f2, "f2old", f2old )
    d   = f + d*(f2/f2old)
    itr = itr+1
    # NOTE: we can always normalize d to 1, but it is not necessary
    return x, f, d, f2

def SolveCG( K, x0, f0, niter=10, eps=1e-6 ):
    print( "===== SolveCG" )
    x  = x0.copy()
    f  = f0 - np.dot(K,x)
    d  = f.copy()
    f2 = np.dot(f,f)
    print( "x:", x )
    print( "f:", f )
    print( "d:", d )
    print( "---- CG loop:", d )
    eps2 = eps**2
    for i in range(niter):
        x, f, d, f2 = stepCG( K, x, d, f, f2 )
        if( f2 < eps2 ): break
    return x

def makeMat( couplings, n ):
    A = np.zeros((n,n))
    A += np.diag( np.ones(n) )
    for (i,j,k) in couplings:
        A[i,j] = k
        A[j,i] = k
    return A

itr=0
def myPrint( x ):
    global itr
    print( itr, "x", x )
    itr = itr+1

# =========== Main

k0 = -1000.0
couplings =[
 ( 0,1, k0*3 ),
 ( 1,2, k0*2 ),
 ( 2,3, k0 ),
]

f0 = np.array([ 
    -1.0, 
     0.0,
     0.0, 
    +1.0,
])
x0 = np.zeros(len(f0))

K = makeMat( couplings, len(f0) )    ;print( "K \n", K )
x = SolveCG( K, x0, f0 )
# print( "x \n", x )



#x_cg = spla.cg(K, f0, x0, tol=1e-05, maxiter=5, callback=myPrint )
x_cg = spla.cg(K, f0, x0, tol=1e-05, maxiter=5 )
print( "x_cg: ", x_cg )
print( "x_ref: ", np.linalg.solve( K, f0 ) )

