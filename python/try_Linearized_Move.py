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
    n    = len(ps)
    Ax   = np.zeros((n,n))
    Ay   = np.zeros((n,n))
    fdlx = np.zeros((n))
    fdly = np.zeros((n))
    Ax  += np.diag( np.ones(n)*kReg )   # regularization, so that the matrix is well conditioned and points does not move too much from their original position
    Ay  += np.diag( np.ones(n)*kReg )
    if constrKs is None: constrKs = np.zeros(n)
    bIsRelaxed = False
    if l0s  is None: bIsRelaxed = True
    ls = np.zeros(len(sticks))
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
        fdl = k*dl
        fdlx[i] += x*fdl
        fdly[i] += y*fdl
        fdlx[j] -= x*fdl
        fdly[j] -= y*fdl
        
        # --- stiffness matrix
        kx = k* np.abs(x)
        ky = k* np.abs(y)
        #kx = k* x
        #ky = k* y
        Ax[i,i] += kx + constrKs[i]
        Ay[i,i] += ky + constrKs[i]
        Ax[j,j] += kx + constrKs[j]
        Ay[j,j] += ky + constrKs[j]
        Ax[i,j] = -kx
        Ay[i,j] = -ky
        Ax[j,i] = -kx
        Ay[j,i] = -ky
    return Ax, Ay, fdlx, fdly, ls

def dynamics( ps, f0, niter = 20, dt=0.1 ):
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
    _, _, _, _, l0s = makeMat_stick_2d_( sticks, ps, constrKs=constrKs )  
    
    plt.figure(figsize=(3*niter,3))
    plt.subplot(1,niter,1)
    plt.plot( ps[:,0], ps[:,1], 'o-k' )
    plt.quiver( ps[:,0], ps[:,1], f0[:,0], f0[:,1] )
    #plt.plot( ps[:,0]+x, ps[:,1]+y, 'o-' )


    for i in range(niter):
        plt.subplot(1,niter,i+1)
        clr = colors[i%len(colors)]

        # ---- Predictor step ( move mass points by external forces )
        # Here we do normal dynamical move v+=(f/m)*dt, p+=v*dt
        f   = f0[:,:] #- ps[:,:]*constrKs[:,None]
        ps += f*dt    # move by external forces (ignoring constraints)   # NOTE: now we use steep descent, but we could verlet or other integrator of equations of motion
        

        # ---- Corrector step ( to satisfy constraints e.g. stick length )
        # ---- Linearize the force around the current position to be able to use CG or other linear solver
        Kx, Ky, fdlx, fdly, ls = makeMat_stick_2d_( sticks, ps, l0s=l0s, constrKs=constrKs, kReg=5.0 )      # fixed end points
        
        print( fdlx, fdly )

        plt.plot( ps[:,0], ps[:,1], 'o:', label=("predicted[%i]" % i), color='k' )
        plt.quiver( ps[:,0], ps[:,1], fdlx , fdly, color='r', scale=1000.0 )

        # ---- Solve for the correction to the position using linear solver
        dx = np.linalg.solve( Kx, f0[:,0]+fdlx )    # solve (f0x+fdlx) = Kx*dx    aka b=A*x ( A=Kx, x=dx, b=f0x+fdlx )
        dy = np.linalg.solve( Ky, f0[:,1]+fdly )    # solve (f0y+fdly) = Ky*dy

        #print( x.shape, y.shape, ps.shape )
        print( "move[%i,%i] " %(i,iCGstep)," |dx|=",  np.linalg.norm(dx), " |dy|=", np.linalg.norm(dy) )
        ps[:,0] += dx
        ps[:,1] += dy
        mask = constrKs>1; ps[ mask,:] = ps0[ mask,:]   # return the constrained points to their original position
        
        #plt.plot( ps[:,0], ps[:,1], 'o-', label=("step[%i]" % i) )
        plt.plot( ps[:,0], ps[:,1], 'o-', label=("corected[%i]" % i), color='k' )

# =========== Main

# ps = np.array([     
#     [-2.0, 0.0],
#     [-1.0,-0.1],
#     [ 0.0,-0.5],
#     [+1.0,-0.1],
#     [+2.0, 0.0],
# ])

# 5 point in line along x-axis
ps = np.array([     
    [-2.0, 0.0],
    [-1.0, 0.0],
    [ 0.0, 0.0],
    [+1.0, 0.0],
    [+2.0, 0.0],
])

# basically a rope with 5 sticks between end points
k0 = 500.0
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


dynamics( ps, f0, niter=10 )

#plt.legend( loc='lower left' )
# plt.xlim(-5,3)
# plt.ylim(-5,3)
# plt.grid()
plt.savefig( "try_Linearized_Move.png", bbox_inches='tight' )
plt.show()
