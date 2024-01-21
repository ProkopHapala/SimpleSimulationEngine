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
    global iCGstep
    iCGstep = 0
    if x0 is None: x0 = np.zeros(len(f0))
    #print( "===== SolveCG" )
    x  = x0.copy()
    
    if dotFunc is not None:
        f = f0 - dotFunc(x)
    else:
        f  = f0 - np.dot(K,x)
    d  = f.copy()
    f2 = np.dot(f,f)
    #if bPrint:
    #print( "[_]x :", x )
    #print( "[_]f0:", f0)
    #print( "[_]f :", f )
    #print( "[_]d :", d )
    #    print( "---- CG loop: f2= ", f2 )
    eps2 = eps**2
    for i in range(niter):
        x, f, d, f2 = stepCG( x, d, f, f2, K=K, dotFunc=dotFunc )
        #print( "CG[%i] x: " %i, x )
        if bPrint: print( "CG[%i] x: " %i, x )
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

def makeNeights( sticks, npoint ):
    neighBs = [ [] for i in range(npoint) ]
    for i,(a,b,k) in enumerate(sticks):
        neighBs[a].append( (b, i) )
        neighBs[b].append( (a, i) )
    nNeighMax = max( [len(neighBs[i]) for i in range(len(neighBs))] )
    return neighBs, nNeighMax

def constr_jacobi_neighs2_absolute(ps, neighBs, ms, l0s ):
    '''
     This is constrain-solver
       we know how long should be each bond
       we need to solve for dpos (i.e. how much we should shift each point to make the bond correct length)
       we shift more the points which are lighter (i.e. we use inverse mass as a weight)
       it should correspond to to position based dynamics (PBD) developed by Muller et al. 2007 (https://matthias-research.github.io/pages/) https://doi.org/10.1016/j.jvcir.2007.01.005  https://matthias-research.github.io/pages/publications/posBasedDyn.pdf
    '''
    dlmax = 0.0
    nPoint = len(ps)
    dpos = np.zeros( ps.shape )
    for iG in range(nPoint):
        pi = ps[iG]
        mi = ms[iG]
        dp = np.zeros(3)
        wsum = 0
        for ja,ib in neighBs[iG]:
            if ib == -1:
                break
            d  = ps[ja,:] - pi
            l  = np.linalg.norm(d)
            dl = l - l0s[ ib ]
            mj = ms[ja]
            w  = mj / ( mj + mi )  # (1/mi)/(1/mi+1/mj) = mi*mj/( mi*(mi + mj)) = mj/(mi+mj)

            #print( "w[%i,%i]=%g (mi=%g,mj=%g)" %( iG, ja, w, mi,mj ) )  
            dp   += d * dl * w * w
            wsum += w
            dl    = np.abs(dl)
            dlmax = max(dl, dlmax)
        dp *= 1.0/wsum
        dpos[iG,:] = dp
    return dpos, dlmax

def apply_sticks( dx, sticks, hdirs, constrKs=None, kReg=1e-2 ):

    n = len(dx)//3
    f = np.zeros((n*3))
    #print( "kReg", kReg )
    #print( "constrKs", constrKs )
    #print( "constrKs_glob ", constrKs_glob )
    for ia in range(n):
        i3 = ia*3
        k  = constrKs[ia] + kReg
        #print( "Kdp[%i] k %g dp(%g,%g,%g)" %( ia, k, dx[i3+0], dx[i3+1], dx[i3+2] ) );
        #k *= 0.0 # DEBUG
        f[i3+0] = dx[i3+0]*k
        f[i3+1] = dx[i3+1]*k
        f[i3+2] = dx[i3+2]*k
    
    #print( "fdpos:", f )
    for ib,( i,j,k) in enumerate(sticks):
        i3 = i*3
        j3 = j*3
        # --- strick vector
        h    = hdirs[ib]; #print( "hdir", hdir )
        dxij = dx[j3:j3+3] - dx[i3:i3+3]
        #k*= 0 # DEBUG
        dl   =  np.dot( h, dxij )
        fb   =  h * ( dl * -k )    # scalar force
        #print( "dot_bond[%i] k %g dl %g h(%g,%g,%g) f(%g,%g,%g)" %( ib, k, dl, h[0],h[1],h[2],   fb[0],fb[1],fb[2]    ) );
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

def stickLenghs( sticks, ps ):
    ls = np.zeros(len(sticks))
    for ib,( i,j,k) in enumerate(sticks):
        x = ps[j,0] - ps[i,0]
        y = ps[j,1] - ps[i,1]
        z = ps[j,2] - ps[i,2]
        ls[ib]  = np.sqrt( x*x + y*y + z*z )
    return ls

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

        #print( "prepare_lin[%i|%i,%i] f0=%g l=%g hdir(%g,%g,%g)"  %( ib,i,j, fdlij, l, hdirs[ib,0], hdirs[ib,1], hdirs[ib,2]) );

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
