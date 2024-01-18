#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse.linalg as spla
from   matplotlib import collections  as mc

# =========== Functions

def makeMat_stick_2d_( sticks, ps, l0s=None, constrKs=None, kReg=1e-2 ):
    '''
    sticks:   list of (i,j,k) where k is the spring constant
    ps:       list of (x,y) coordinates of the points
    constrKs: list of spring constants constraining the points in place (i.e. fixed points), this is important to ensure that the matrix is well conditioned
    kReg:     regularization constant to ensure that the matrix is well conditioned (e.g. if there are no constrKs)
    '''
    n   = len(ps)
    fdl = np.zeros((n*2))
    A   = np.zeros((2*n,2*n))
    A  += np.diag( np.ones(2*n)*kReg )   # regularization, so that the matrix is well conditioned and points does not move too much from their original position
    if constrKs is None: constrKs = np.zeros(n)
    bIsRelaxed = False
    if l0s  is None: bIsRelaxed = True
    ls = np.zeros(len(sticks))

    for i in range(2*n):
        A[i,i] += constrKs[i//2]
    
    for ib,( i,j,k) in enumerate(sticks):

        # --- strick vector
        x = ps[j,0] - ps[i,0]
        y = ps[j,1] - ps[i,1]
        l  = np.sqrt( x*x + y*y )

        # --- stick length and normalized stick direction
        ls[ib] = l
        il = 1./l
        x*=il
        y*=il

        # --- force due to change of stick length
        dl = 0.0
        if not bIsRelaxed: 
            dl  = ls[ib] - l0s[ib]
        fdlij = k*dl

        # construct the system in the standard form Ax=b
        i2 = i*2
        j2 = j*2

        # force vector ( apply sticks stress to the points )
        fdl[i2+0] += x*fdlij
        fdl[i2+1] += y*fdlij
        fdl[j2+0] -= x*fdlij
        fdl[j2+1] -= y*fdlij
        
        # --- stiffness matrix
        # diagonal i,i
        A[i2+0,i2+0] += k*x*x
        A[i2+1,i2+0] += k*y*x
        A[i2+0,i2+1] += k*x*y
        A[i2+1,i2+1] += k*y*y
        # diagonal j,j
        A[j2+0,j2+0] += k*x*x
        A[j2+1,j2+0] += k*y*x
        A[j2+0,j2+1] += k*x*y
        A[j2+1,j2+1] += k*y*y
        # off-diagonal i,j
        A[i2+0,j2+0] -= k*x*x
        A[i2+1,j2+0] -= k*y*x
        A[i2+0,j2+1] -= k*x*y
        A[i2+1,j2+1] -= k*y*y
        # off-diagonal j,i
        A[j2+0,i2+0] -= k*x*x
        A[j2+1,i2+0] -= k*y*x
        A[j2+0,i2+1] -= k*x*y
        A[j2+1,i2+1] -= k*y*y
    
    return A, fdl, ls

def dynamics( ps, f0, niter = 10, dt=0.05 ):
    #cmap   = plt.get_cmap('rainbow')
    cmap   = plt.get_cmap('gist_rainbow')
    #cmap   = plt.get_cmap('jet')
    #cmap   = plt.get_cmap('turbo')
    colors = [cmap(i/float(niter)) for i in range(niter)]
    global iCGstep
    n = len(ps)
    iCGstep = 0
    constrKs=np.array([50.0, 0.0, 0.0, 0.0,50.0])
    ps0 = ps.copy()
    _, _, l0s = makeMat_stick_2d_( sticks, ps, constrKs=constrKs )  
    
    plt.figure(figsize=(3*niter,3))
    plt.subplot(1,niter,1)
    plt.plot( ps[:,0], ps[:,1], 'o-k' )
    plt.quiver( ps[:,0], ps[:,1], f0[:,0], f0[:,1] )
    #plt.plot( ps[:,0]+x, ps[:,1]+y, 'o-' )

    v = np.zeros((n,2))

    for i in range(niter):
        plt.subplot(1,niter,i+1)
        clr = colors[i%len(colors)]

        # ---- Predictor step ( move mass points by external forces )
        # Here we do normal dynamical move v+=(f/m)*dt, p+=v*dt
        f   = f0[:,:] #- ps[:,:]*constrKs[:,None]
        v  += f*dt    # move by external forces (ignoring constraints)   # NOTE: now we use steep descent, but we could verlet or other integrator of equations of motion
        ps += v*dt    # move by external forces (ignoring constraints)   # NOTE: now we use steep descent, but we could verlet or other integrator of equations of motion
        

        # ---- Corrector step ( to satisfy constraints e.g. stick length )
        # ---- Linearize the force around the current position to be able to use CG or other linear solver
        K, fdl, ls = makeMat_stick_2d_( sticks, ps, l0s=l0s, constrKs=constrKs, kReg=5.0 )      # fixed end points
        
        #print( "K\n", K )

        f = f0.flatten() + fdl

        fx = f[0::2]
        fy = f[1::2]
        plt.plot  ( ps[:,0], ps[:,1], 'o:', label=("predicted[%i]" % i), color='k' )
        plt.quiver( ps[:,0], ps[:,1], fx , fy, color='r', scale=1000.0 )

        # ---- Solve for the correction to the position using linear solver
        dp = np.linalg.solve( K, f )    # solve (f0x+fdlx) = Kx*dx    aka b=A*x ( A=Kx, x=dx, b=f0x+fdlx )

        #print( x.shape, y.shape, ps.shape )
        print( "move[%i,%i] " %(i,iCGstep)," |d|=",  np.linalg.norm(dp) )
        ps[:,0] += dp[0::2]
        ps[:,1] += dp[1::2]
        mask = constrKs>1; ps[ mask,:] = ps0[ mask,:]   # return the constrained points to their original position
        
        #plt.plot( ps[:,0], ps[:,1], 'o-', label=("step[%i]" % i) )
        plt.plot( ps[:,0], ps[:,1], 'o-', label=("corected[%i]" % i), color='k' )
        #plt.xlim(-5,5); 
        plt.ylim(-3,1)

# =========== Main

# 5 point in line along x-axis
ps = np.array([     
    [-2.0, 0.0],
    [-1.0, 0.0],
    [ 0.0, 0.0],
    [+1.0, 0.0],
    [+2.0, 0.0],
])

# basically a rope with 5 sticks between end points
k0 = 5000.0
sticks =[
 ( 0,1, k0 ),
 ( 1,2, k0 ),
 ( 2,3, k0 ),
 ( 3,4, k0 ),
]

# pull down the middle point
f0 = np.array([ 
[ 0.0, 0.0],
[ 0.0, 0.0],
[ 0.0,-5.0],   # pull down the middle point
[ 0.0, 0.0],
[ 0.0, 0.0],
])


dynamics( ps, f0 )

#plt.legend( loc='lower left' )
#plt.xlim(-5,5); plt.ylim(-5,5)
# plt.grid()
plt.savefig( "try_Linearized_Move.png", bbox_inches='tight' )
plt.show()
