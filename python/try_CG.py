#!/usr/bin/python

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

def makeMat( couplings, n ):
    A = np.zeros((n,n))
    #A += np.diag( np.ones(n) )
    for (i,j,k) in couplings:
        A[i,j]  = -k
        A[j,i]  = -k
        A[i,i] +=  k
        A[j,j] +=  k
    return A

def makeMat_stick_2d( sticks, ps ):
    n = len(ps)
    A = np.zeros((n*2,n*2))
    #A += np.diag( np.ones(n*2) )
    for ( i,j,k) in sticks:
        x = ps[j,0] - ps[i,0]
        y = ps[j,1] - ps[i,1]
        il = 1./np.sqrt( x*x + y*y )
        kx = k*x*il
        ky = k*y*il
        i2 = i*2
        j2 = j*2
        A[i2+0,i2+0] += kx
        A[i2+1,i2+1] += ky

        A[j2+0,j2+0] += kx
        A[j2+1,j2+1] += ky

        A[i2+0,j2+0] = -kx
        A[i2+1,j2+1] = -ky

        A[j2+0,i2+0] = -kx
        A[j2+1,i2+1] = -ky

    return A

def makeMat_stick_2d_( sticks, ps, l0s=None, constrKs=None, kReg=1e-2 ):
    '''
    sticks: list of (i,j,k) where k is the spring constant
    ps:     list of (x,y) coordinates of the points
    constrKs: list of spring constants constraining the points in place (i.e. fixed points), this is important to ensure that the matrix is well conditioned
    kReg:   regularization constant to ensure that the matrix is well conditioned (e.g. if there are no constrKs)
    '''
    n    = len(ps)
    Ax   = np.zeros((n,n))
    Ay   = np.zeros((n,n))
    fdlx = np.zeros((n))
    fdly = np.zeros((n))
    Ax  += np.diag( np.ones(n)*kReg )
    Ay  += np.diag( np.ones(n)*kReg )
    if constrKs is None: constrKs = np.zeros(n)
    bIsRelaxed = False
    if l0s  is None: 
        bIsRelaxed = True
    ls = np.zeros(len(sticks))
    for ib,( i,j,k) in enumerate(sticks):
        x = ps[j,0] - ps[i,0]
        y = ps[j,1] - ps[i,1]
        l  = np.sqrt( x*x + y*y )
        ls[ib] = l
        il = 1./l
        x*=il
        y*=il
        dl = 0.0
        if not bIsRelaxed: 
            dl  = ls[ib] - l0s[ib]
        fdl = k*dl
        fdlx[i] += x*fdl
        fdly[i] += y*fdl
        fdlx[j] -= x*fdl
        fdly[j] -= y*fdl
        
        kx = k* np.abs(x)
        ky = k* np.abs(y)
        i2 = i*2
        j2 = j*2
        Ax[i,i] += kx + constrKs[i]
        Ay[i,i] += ky + constrKs[i]

        Ax[j,j] += kx + constrKs[j]
        Ay[j,j] += ky + constrKs[j]

        Ax[i,j] = -kx
        Ay[i,j] = -ky

        Ax[j,i] = -kx
        Ay[j,i] = -ky
    return Ax, Ay, fdlx, fdly, ls

def move_CG( ps, f0, l0s, nitr = 20, dt=0.01, nCGmax=5, fCGconv=1e-3  ):
    global iCGstep
    n = len(ps)
    iCGstep = 0
    constrKs=np.array([50.0, 0.0, 0.0, 0.0,50.0])
    ps0 = ps.copy()
    _, _, _, _, ls = makeMat_stick_2d_( sticks, ps, l0s=l0s, constrKs=constrKs )  
    for i in range(nitr):
        # ---- Here we do normal dynamical move v+=(f/m)*dt, p+=v*dt
        f = f0[:,:] #- ps[:,:]*constrKs[:,None]
        ps += f*dt
        mask = constrKs>1
        ps[ mask,:] = ps0[ mask,:]
        #plt.plot( ps[:,0], ps[:,1], 'o:', label=("step[%i]" % i) )

        # ---- Linearize the force around the current position to be able to use CG
        Kx, Ky, fdlx, fdly, ls = makeMat_stick_2d_( sticks, ps, l0s=l0s, constrKs=constrKs )      # fixed end points
        # ---- CG is like corrector step (to ensure Truss contrains)
        # x, f2x = SolveCG( Kx, f0[:,0]+fdlx, niter=nCGmax, eps=fCGconv )
        # y, f2y = SolveCG( Ky, f0[:,1]+fdly, niter=nCGmax, eps=fCGconv )

        # # ---- Jacobi is like corrector step (to ensure Truss contrains)
        # x, f2x = Jacobi( Kx, f0[:,0]+fdlx, niter=nCGmax, eps=fCGconv )
        # y, f2y = Jacobi( Ky, f0[:,1]+fdly, niter=nCGmax, eps=fCGconv )

        # ---- SOR is like corrector step (to ensure Truss contrains)
        x, f2x = SOR( Kx, f0[:,0]+fdlx, niter=nCGmax, eps=fCGconv )
        y, f2y = SOR( Ky, f0[:,1]+fdly, niter=nCGmax, eps=fCGconv )

        #print( x.shape, y.shape, ps.shape )
        print( "move[%i,%i] |x|" %(i,iCGstep),  np.linalg.norm(x), "|y|", np.linalg.norm(y), "fCGx:", np.sqrt(f2x),"fCGx:", np.sqrt(f2y) )
        ps[:,0] += x
        ps[:,1] += y
        
        plt.plot( ps[:,0], ps[:,1], 'o-', label=("step[%i]" % i) )

# =========== Main
'''
k0 = 1000.0
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
'''

ps = np.array([     
    [-2.0, 0.0],
    [-1.0,-0.1],
    [ 0.0,-0.5],
    [+1.0,-0.1],
    [+2.0, 0.0],
])

f0 = np.array([ 
[-0.5*0,+0.5*0],
[ 0.0, 0.0],
[ 0.0,-5.0],
[ 0.0, 0.0],
[ 0.5*0,+0.5*0],
])

k0 = 50.0
sticks =[
 ( 0,1, k0 ),
 ( 1,2, k0 ),
 ( 2,3, k0 ),
 ( 3,4, k0 ),
]
'''
f0 = f0.flatten()
x0 = np.zeros(len(f0))
K = makeMat_stick_2d(sticks , ps )    ;print( "K \n", K )
print( K.shape, f0.shape, x0.shape )
# x = SolveCG( K, x0, f0 )
# print( "x \n", x )
'''



#Kx, Ky = makeMat_stick_2d_( sticks, ps )     # regularization by homogeneous constrain

Kx, Ky, fdlx, fdly, ls = makeMat_stick_2d_( sticks, ps, constrKs=[50.0, 0.0, 0.0, 0.0,50.0] )      # fixed end points

print( "Kx \n", Kx )
print( "Ky \n", Ky )

print( " ====== solve Kx" )
x = SolveCG( Kx,  f0[:,0] )
print( "f0_x : ", f0[:,0] )
print( "x_CG : ", x  )
print( "x_ref: ", np.linalg.solve( Kx, f0[:,0] ) )

print( " ====== solve Ky" )
y = SolveCG( Ky,  f0[:,1] )
print( "f0_y : ", f0[:,1] )
print( "y_CG : ", y  )
print( "y_ref: ", np.linalg.solve( Ky, f0[:,1] ) )


# plot arrow for each point in direction of force

   

plt.figure()
plt.plot( ps[:,0], ps[:,1], 'o-k' )
plt.quiver( ps[:,0], ps[:,1], f0[:,0], f0[:,1] )
#plt.plot( ps[:,0]+x, ps[:,1]+y, 'o-' )

move_CG( ps, f0, ls )

plt.legend( loc='lower left' )
plt.xlim(-5,3)
plt.ylim(-5,3)
plt.grid()
plt.show()
