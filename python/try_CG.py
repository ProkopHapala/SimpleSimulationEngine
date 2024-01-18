#!/usr/bin/python

import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse.linalg as spla
from   matplotlib import collections  as mc

# =========== Functions
iCGstep=0
bPrint = False

def stepCG( x, d, f, f2old, K=None, dotFunc=None ):
    global iCGstep
    if dotFunc is not None:
        Kd  = dotFunc(d)
    else:
        Kd = np.dot(K,d)
    if bPrint:
        print("[%i]r :" %(iCGstep-1), f  );
        print("[%i]p :" %(iCGstep-1), d  );
        print("[%i]Ap:" %(iCGstep-1), Kd );
    dt  = np.dot(d,f) / np.dot(d,Kd)    # step length such that f is orthogonal to d 
    #dt  = np.dot(f,f) / np.dot(d,Kd)
    if bPrint:
        print( "### CG_step ", iCGstep, " dt=", dt, " rho=", f2old )
        print( "[%i]Kd:" %iCGstep, Kd )
        #print( "ox", x )
        #print( "of", f )
        print( "dt", dt )
    x   = x + d *dt
    f   = f - Kd*dt
    if bPrint:
        print( "[%i]x :"%iCGstep, x )
        print( "[%i]f :"%iCGstep, f )
    f2  = np.dot(f,f)
    #print( "CG[%i]" %iCGstep, "dt ", dt, "f2", f2, "f2old", f2old )
    d   = f + d*(f2/f2old)
    iCGstep = iCGstep+1
    # NOTE: we can always normalize d to 1, but it is not necessary
    return x, f, d, f2

def SolveCG( f0, K=None, dotFunc=None, x0=None, niter=10, eps=1e-6, bPrint=False ):
    if x0 is None: x0 = np.zeros(len(f0))
    #print( "===== SolveCG" )
    x  = x0.copy()
    if dotFunc is not None:
        f = f0 - dotFunc(x)
    else:
        f  = f0 - np.dot(K,x)
    d  = f.copy()
    f2 = np.dot(f,f)
    if bPrint:
        print( "x:", x )
        print( "f:", f )
        print( "d:", d )
        print( "---- CG loop: f2= ", f2 )
    eps2 = eps**2
    for i in range(niter):
        x, f, d, f2 = stepCG( x, d, f, f2, K=K, dotFunc=dotFunc )
        if bPrint: print( "CG[%i] x: ", x )
        if( f2 < eps2 ): 
            if bPrint: print ("### CG converged at step", i )
            break
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

def apply_sticks( dx, sticks, hdirs, constrKs=None, kReg=1e-2 ):

    n = len(dx)//3
    f = np.zeros((n*3))
    #print( "kReg", kReg )
    #print( "constrKs", constrKs )
    #print( "constrKs_glob ", constrKs_glob )
    for ia in range(n):
        i3 = ia*3
        k  = constrKs[ia] + kReg
        #k *= 0.0 # DEBUG
        f[i3+0] = dx[i3+0]*k
        f[i3+1] = dx[i3+1]*k
        f[i3+2] = dx[i3+2]*k
    
    for ib,( i,j,k) in enumerate(sticks):
        i3 = i*3
        j3 = j*3
        # --- strick vector
        h    = hdirs[ib]; #print( "hdir", hdir )
        dxij = dx[j3:j3+3] - dx[i3:i3+3]
        #k*= 0 # DEBUG
        dl   =  np.dot( h, dxij )
        fb   =  h * ( dl * -k )    # scalar force
        #print( "dot_bond[%i] k %g fl %g h(%g,%g,%g)" %( ib, k, dl, h[0], h[1], h[2]) );
        # apply force to the points
        # point i
        f[i3+0] += fb[0]
        f[i3+1] += fb[1]
        f[i3+2] += fb[2]
        # point j
        f[j3+0] -= fb[0]
        f[j3+1] -= fb[1]
        f[j3+2] -= fb[2]
    return f

def globalize( sticks, hdirs ):
    global sticks_glob, hdirs_glob #, constrKs_glob, kReg_glob
    sticks_glob    = sticks
    hdirs_glob     = hdirs
    #constrKs_glob  = np.zeros(len(hdirs)//2)
    #kReg_glob      = 1e-2

def dot_func( dx ):
    return apply_sticks( dx, sticks_glob, hdirs_glob, constrKs=constrKs_glob, kReg=kReg_glob )

def linearize( sticks, ps, l0s=None, constrKs=None, kReg=1e-2 ):
    n   = len(ps)
    fdl = np.zeros((n*3))
    if constrKs is None: constrKs = np.zeros(n)
    bIsRelaxed = False
    if l0s  is None: bIsRelaxed = True
    ls = np.zeros(len(sticks))

    hdirs = np.zeros((len(sticks),3))
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
        hdirs[ib,:] = [x,y,z]

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
        
    return fdl, ls, hdirs


def makeMat( sticks, ps, l0s=None, constrKs=None, kReg=1e-2 ):
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

    for i in range(3*n):
        A[i,i] += constrKs[i//3]

    #A[:,:] = 0.0 # DEBUG

    #print( "kReg", kReg )
    #print( "constrKs", constrKs )

    hdirs = np.zeros((len(sticks),3))
    
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
        hdirs[ib,:] = [x,y,z]

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
        

        #k*=0.0 # DEBUG
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
    return A, fdl, ls, hdirs

def dynamics( ps, f0, niter = 5, dt=0.05, constrKs=None, kReg=5.0, bUseCG=True ):
    #cmap   = plt.get_cmap('rainbow')
    cmap   = plt.get_cmap('gist_rainbow')
    #cmap   = plt.get_cmap('jet')
    #cmap   = plt.get_cmap('turbo')
    colors = [cmap(i/float(niter)) for i in range(niter)]
    global iCGstep
    n = len(ps)
    iCGstep = 0
    if constrKs is None: constrKs = np.zeros(n)
    #constrKs=np.array([50.0, 0.0, 0.0, 0.0,50.0])
    ps0 = ps.copy()
    _, _, l0s, hdirs = makeMat( sticks, ps, constrKs=constrKs, kReg=kReg )  
    
    global constrKs_glob, kReg_glob, sticks_glob, hdirs_glob
    constrKs_glob  = constrKs
    kReg_glob      = kReg
    sticks_glob    = sticks

    plt.figure(figsize=(3*niter,3))
    plt.subplot(1,niter,1)
    plt.plot(   ps[:,0], ps[:,1], 'o-k' )
    plt.quiver( ps[:,0], ps[:,1], f0[:,0], f0[:,1] )
    #plt.plot( ps[:,0]+x, ps[:,1]+y, 'o-' )

    v = np.zeros((n,3))

    for i in range(niter):
        plt.subplot(1,niter,i+1)
        clr = colors[i%len(colors)]

        # ---- Predictor step ( move mass points by external forces )
        # Here we do normal dynamical move v+=(f/m)*dt, p+=v*dt
        f   = f0[:,:] #- ps[:,:]*constrKs[:,None]
        v  += f*dt    # move by external forces (ignoring constraints)   # NOTE: now we use steep descent, but we could verlet or other integrator of equations of motion
        ps += v*dt    # move by external forces (ignoring constraints)   # NOTE: now we use steep descent, but we could verlet or other integrator of equations of motion
        
        if bUseCG:
            # K, fdl, ls, hdirs = makeMat( sticks, ps, l0s=l0s, constrKs=constrKs, kReg=kReg )      # fixed end points
            # f  = f0.flatten() + fdl
            # #dp = np.linalg.solve( K, f )    # solve (f0x+fdlx) = Kx*dx    aka b=A*x ( A=Kx, x=dx, b=f0x+fdlx )
            # dp, err2 = SolveCG( f, niter=10, eps=1e-6, bPrint=True, K=K )

            fdl, ls, hdirs = linearize( sticks, ps, l0s=l0s, constrKs=constrKs, kReg=kReg )
            f = f0.flatten() + fdl
            hdirs_glob = hdirs
            dp, err2 = SolveCG( f, niter=10, eps=1e-6, bPrint=False, dotFunc=dot_func )

            #exit()
        else:
            K, fdl, ls, hdirs = makeMat( sticks, ps, l0s=l0s, constrKs=constrKs, kReg=kReg )      # fixed end points
            f  = f0.flatten() + fdl
            dp = np.linalg.solve( K, f )    # solve (f0x+fdlx) = Kx*dx    aka b=A*x ( A=Kx, x=dx, b=f0x+fdlx )
            #dp, err2 = SolveCG( f, niter=10, eps=1e-6, bPrint=False, K=K )
            
            #print( "K\n", K )

        fx = f[0::3]
        fy = f[1::3]
        plt.plot  ( ps[:,0], ps[:,1], 'o:', label=("predicted[%i]" % i), color='k' )
        plt.quiver( ps[:,0], ps[:,1], fx , fy, color='r', scale=1000.0 )

        #print( x.shape, y.shape, ps.shape )
        print( "move[%i,%i] " %(i,iCGstep)," |d|=",  np.linalg.norm(dp) )
        ps[:,0] += dp[0::3]
        ps[:,1] += dp[1::3]
        mask = constrKs>1; ps[ mask,:] = ps0[ mask,:]   # return the constrained points to their original position
        
        #plt.plot( ps[:,0], ps[:,1], 'o-', label=("step[%i]" % i) )
        plt.plot( ps[:,0], ps[:,1], 'o-', label=("corected[%i]" % i), color='k' )
        #plt.xlim(-5,5); 
        plt.ylim(-3,1)

# =========== Main

ps = np.array([     
    [-2.0, -0.0, 0.0],
    [-1.0, -0.1, 0.0],
    [ 0.0, -0.2, 0.0],
    [+1.0, -0.1, 0.0],
    [+2.0, -0.0, 0.0],
])

k0 = 50.0
sticks =[
 ( 0,1, k0 ),
 ( 1,2, k0 ),
 ( 2,3, k0 ),
 ( 3,4, k0 ),
]

f0 = np.array([ 
[ 0.0, 0.0, 0.0],
[ 0.0, 0.0, 0.0],
[ 0.0,-5.0, 0.0],
[ 0.0, 0.0, 0.0],
[ 0.0, 0.0, 0.0],
])

constrKs = np.array([50.0, 0.0, 0.0, 0.0,50.0])
kReg     = 5.0


dynamics( ps, f0, constrKs=constrKs, kReg=kReg )

'''
#constrKs = [0.0, 0.0, 0.0, 0.0,0.0]
#kReg     = 0.0

os.system('mode con: cols=100 lines=50')
np.set_printoptions(linewidth=np.inf)

dp = np.zeros(ps.shape )
dp[0,1] = 1.0
dp[2,1] = 1.0
dp[4,1] = 0.5

K, fdl, ls, hdirs = makeMat( sticks, ps, constrKs=constrKs, kReg=kReg )      # fixed end points

x = dp.flat.copy()
f  = np.dot( K, x )
ff = apply_sticks( x, sticks, hdirs, constrKs=constrKs, kReg=kReg )

#print( "K: \n", K )
print( "x    : ", x )
print( "f_mat: ", f )
print( "f_fun: ", ff )
#exit(0)

print("====================set initial conditions");
print( "x0", np.zeros(len(fdl))  )
print( "f0", f0.flat+fdl         )
print("==================== SOLVE");
x, err2 = SolveCG( K, f0.flat+fdl, x0=None, niter=10, eps=1e-6 )
print( "x_CG ", x )
print( "x_ref", np.linalg.solve( K, f0.flat+fdl ) )
'''


plt.savefig( "try_Linearized_Move.png", bbox_inches='tight' )
plt.show()
