#!/usr/bin/python

import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse.linalg as spla
from   matplotlib import collections  as mc

# =========== Functions
iCGstep=0

def stepCG( K, x, d, f, f2old ):
    global iCGstep
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
    #print( "x:", x, " f:", f )
    #print( "f", f )
    f2  = np.dot(f,f)
    print( "CG[%i]" %iCGstep, "dt ", dt, "f2", f2, "f2old", f2old )
    d   = f + d*(f2/f2old)
    iCGstep = iCGstep+1
    # NOTE: we can always normalize d to 1, but it is not necessary
    return x, f, d, f2

def SolveCG( K, f0, x0=None, niter=10, eps=1e-6 ):
    if x0 is None: x0 = np.zeros(len(f0))
    #print( "===== SolveCG" )
    x  = x0.copy()
    f  = f0 - np.dot(K,x)
    d  = f.copy()
    f2 = np.dot(f,f)
    #print( "x:", x )
    #print( "f:", f )
    #print( "d:", d )
    #print( "---- CG loop:", d )
    eps2 = eps**2
    for i in range(niter):
        x, f, d, f2 = stepCG( K, x, d, f, f2 )
        if( f2 < eps2 ): break
    return x, f2

def Jacobi( K, f0, x0=None, niter=10, eps=1e-6 ):
    n = len(f0)
    if x0 is None: 
        x = np.zeros(n)
    else:
        x = x0.copy()
    for i in range(niter):
        f2 = 0.0
        for i in range(n):
            fi = f0[i] - np.dot( K[i,:], x )
            x[i] =  fi / K[i,i]
            f2 += fi*fi
        if( f2 < eps**2 ): break
    return x, f2

def SOR( K, f0, x0=None, niter=10, eps=1e-6, w=0.9 ):
    n = len(f0)
    if x0 is None: 
        x = np.zeros(n)
    else:
        x = x0.copy()
    for i in range(niter):
        f2 = 0.0
        for i in range(n):
            fi = f0[i] - np.dot( K[i,:], x )
            x[i] = x[i]*(1.-w) + fi * (w/K[i,i])
        f2 += fi*fi
        if( f2 < eps**2 ): break
    return x, f2

def makeMat_stick_2d_( sticks, ps, l0s=None, constrKs=None, kReg=1e-2 ):
    '''
    sticks:   list of (i,j,k) where k is the spring constant
    ps:       list of (x,y) coordinates of the points
    constrKs: list of spring constants constraining the points in place (i.e. fixed points), this is important to ensure that the matrix is well conditioned
    kReg:     regularization constant to ensure that the matrix is well conditioned (e.g. if there are no constrKs)
    '''
    n   = len(ps)
    fdl = np.zeros((n*3))
    A   = np.zeros((n*3,n*3))
    A  += np.diag( np.ones(3*n)*kReg )   # regularization, so that the matrix is well conditioned and points does not move too much from their original position
    if constrKs is None: constrKs = np.zeros(n)
    bIsRelaxed = False
    if l0s  is None: bIsRelaxed = True
    ls = np.zeros(len(sticks))

    for i in range(2*n):
        A[i,i] += constrKs[i//3]
    
    for ib,( i,j,k) in enumerate(sticks):

        # --- strick vector
        x = ps[j,0] - ps[i,0]
        y = ps[j,1] - ps[i,1]
        z = ps[j,2] - ps[i,2]
        l  = np.sqrt( x*x + y*y + z*z )

        # --- stick length and normalized stick direction
        ls[ib] = l
        il = 1./l
        x*=il
        y*=il
        z*=il

        # --- force due to change of stick length
        dl = 0.0
        if not bIsRelaxed: 
            dl  = ls[ib] - l0s[ib]
        fdlij = k*dl

        # construct the system in the standard form Ax=b
        i3 = i*3
        j3 = j*3
        
        # force vector ( apply sticks stress to the points )
        fdl[i3+0] += x*fdlij
        fdl[i3+1] += y*fdlij
        fdl[i3+2] += z*fdlij
        fdl[j3+0] -= x*fdlij
        fdl[j3+1] -= y*fdlij
        fdl[j3+2] -= z*fdlij
        
        # --- stiffness matrix
        # diagonal i,i
        A[i3+0,i3+0] += k*x*x
        A[i3+1,i3+0] += k*y*x
        A[i3+2,i3+0] += k*z*x
        A[i3+0,i3+1] += k*x*y
        A[i3+1,i3+1] += k*y*y
        A[i3+2,i3+1] += k*z*y
        A[i3+0,i3+2] += k*x*z
        A[i3+1,i3+2] += k*y*z
        A[i3+2,i3+2] += k*z*z
        # diagonal j,j
        A[j3+0,j3+0] += k*x*x
        A[j3+1,j3+0] += k*y*x
        A[j3+2,j3+0] += k*z*x
        A[j3+0,j3+1] += k*x*y
        A[j3+1,j3+1] += k*y*y
        A[j3+2,j3+1] += k*z*y
        A[j3+0,j3+2] += k*x*z
        A[j3+1,j3+2] += k*y*z
        A[j3+2,j3+2] += k*z*z
        # off-diagonal i,j
        A[i3+0,j3+0] -= k*x*x
        A[i3+1,j3+0] -= k*y*x
        A[i3+2,j3+0] -= k*z*x
        A[i3+0,j3+1] -= k*x*y
        A[i3+1,j3+1] -= k*y*y
        A[i3+2,j3+1] -= k*z*y
        A[i3+0,j3+2] -= k*x*z
        A[i3+1,j3+2] -= k*y*z
        A[i3+2,j3+2] -= k*z*z
        # off-diagonal j,i
        A[j3+0,i3+0] -= k*x*x
        A[j3+1,i3+0] -= k*y*x
        A[j3+2,i3+0] -= k*z*x
        A[j3+0,i3+1] -= k*x*y
        A[j3+1,i3+1] -= k*y*y
        A[j3+2,i3+1] -= k*z*y
        A[j3+0,i3+2] -= k*x*z
        A[j3+1,i3+2] -= k*y*z
        A[j3+2,i3+2] -= k*z*z
    return A, fdl, ls

# =========== Main

ps = np.array([     
    [-2.0, -0.0, 0.0],
    [-1.0, -0.1, 0.0],
    [ 0.0, -0.2, 0.0],
    [+1.0, -0.1, 0.0],
    [+2.0, -0.0, 0.0],
])

f0 = np.array([ 
[ 0.0, 0.0, 0.0],
[ 0.0, 0.0, 0.0],
[ 0.0,-5.0, 0.0],
[ 0.0, 0.0, 0.0],
[ 0.0, 0.0, 0.0],
])

k0 = 50.0
sticks =[
 ( 0,1, k0 ),
 ( 1,2, k0 ),
 ( 2,3, k0 ),
 ( 3,4, k0 ),
]

constrKs=[50.0, 0.0, 0.0, 0.0,50.0]

os.system('mode con: cols=100 lines=50')
np.set_printoptions(linewidth=np.inf)

dp = np.zeros(ps.shape )
dp[0,1] = 1.0
dp[2,1] = 1.0

K, fdl, ls = makeMat_stick_2d_( sticks, ps, constrKs=constrKs )      # fixed end points

print( "K: \n", K )
x = dp.flat.copy()
f = np.dot( K, x )
print( "x: ", x )
print( "f: ", f )
exit(0)

K, fdl, ls = makeMat_stick_2d_( sticks, ps, constrKs=constrKs )      # fixed end points
x, err2 = SolveCG( K, f0.flat+fdl, x0=None, niter=10, eps=1e-6 )
print( "x_CG ", x )
print( "x_def", np.linalg.solve( K, f0.flat+fdl ) )