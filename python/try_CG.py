#!/usr/bin/python

import os
import numpy as np
import matplotlib.pyplot as plt
import scipy.sparse.linalg as spla
from   matplotlib import collections  as mc

import constrain_solver as cs

def drawStricks( ps, sticks, ax=None, color='k', lw=1, ls='-' ):
    if ax is None: ax = plt.gca()
    lines = []
    for i1,i2,k in sticks:
        p1 = ps[i1][:2]
        p2 = ps[i2][:2]
        lines.append( [p1,p2] )
    lc = mc.LineCollection(lines, colors=color, linewidths=lw, linestyle=ls  )
    ax.add_collection(lc)
    #ax.axis('equal')


def dynamics_PBD_GS( ps, ms, f0, sticks, niter=1, dt=0.2, nsolvmax=20, dlconv=1e-3 ):
    cmap   = plt.get_cmap('rainbow')
    colors = [cmap(i/float(nsolvmax)) for i in range(nsolvmax)]
    n = len(ps)
    plt.figure(figsize=(3*niter,3))
    plt.subplot(1,niter,1)
    plt.plot(   ps[:,0], ps[:,1], 'o-k' )

    l0s                = cs.stickLenghs( sticks, ps       )
    neighBs, nNeighMax = cs.makeNeights( sticks, npoint=n )
    v = np.zeros((n,3))

    vdp = np.zeros((n,3))

    #print( "l0s", l0s )
    #print( "ms",  ms  )
    #for i,b in enumerate(neighBs): 
    #    print("neighBs[%i]" %i, b)
    #exit()

    # f0[2,1] = -5.0
    # f0[2,0] =  5.0

    dp_sc = 1.2
    sc_1  = 2.0
    #bmix  = 0.7
    bmix_anti = 0.0
    bmix_syn  = 1.0

    for i in range(niter):
        print( "--------- move[%i]" %i );
        plt.subplot(1,niter,i+1)
        
        v[:,:]  += f0[:,:]*dt   
        ps_bak = ps.copy()
        #plt.plot( ps_bak[:,0], ps_bak[:,1], 'o--k', label=("start[%i]" % i) )
        drawStricks( ps_bak, sticks, color='k' )
        ps[:,:] += v[:,:]*dt   
        ps_mov = ps.copy()
        #plt.plot( ps_mov[:,0], ps_mov[:,1], 'o:k',  label=("moved[%i]" % i) )
        drawStricks( ps_mov, sticks, color='k', ls = ':' )
        
        vdp[:,:] = 0.0
        
        for isolv in range(nsolvmax):
            dps, dlmax = cs.constr_jacobi_neighs2_absolute( ps, neighBs, ms, l0s )
            
            #vdp[:,:] = dps[:,:]*1.5
            #vdp[:,:] = dps[:,:]
            cvf = 0  

            # if isolv==0:
            #     vdp[:,:]= dps*sc_1
            # else:
            #     dps *= dp_sc
            #     vdp  = dps

            if isolv==0:
                vdp[:,:]= dps*sc_1
            else:
                # vv  = np.dot( dps.flat, dps.flat )
                # dd  = np.dot( vdp.flat, vdp.flat )               
                cvf = np.dot( dps.flat, vdp.flat ) #/np.sqrt(vv*dd)
                # #vdp += dps*(-cvf/dd)
                # vdp *=0.5
                # vdp += dps*dp_sc
                if cvf<0:
                    vdp *= bmix_anti
                    vdp += dps*dp_sc
                else:
                    vdp *= bmix_syn
                    vdp += dps*sc_1

            print( "solve[%i] slmax=%g cvf=%g" %(isolv, dlmax, cvf) )

            plt.quiver( ps[:,0], ps[:,1], vdp[:,0], vdp[:,1], scale=1., scale_units='xy', width=0.002 )
            #plt.quiver( ps[:,0], ps[:,1], dps[:,0], dps[:,1] )
            ps[:,:] += vdp
            drawStricks( ps, sticks, color=colors[isolv] )
            #plt.plot( ps    [:,0], ps[:,1], 'o-', c=colors[isolv], label=("final[%i|%i]" %(i,isolv)) )
            if dlmax<dlconv: break
        v = (ps-ps_bak)/dt

        plt.legend( loc='upper left' )
        plt.axis('equal')

def dynamics( ps, f0, niter = 3, dt=0.05, constrKs=None, kReg=5.0, bUseCG=True ):
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
    _, _, l0s, hdirs = cs.makeMat( sticks, ps, constrKs=constrKs, kReg=kReg )  
    
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
        print( "--------- move[%i]" %i );
        plt.subplot(1,niter,i+1)
        clr = colors[i%len(colors)]

        # ---- Predictor step ( move mass points by external forces )
        # Here we do normal dynamical move v+=(f/m)*dt, p+=v*dt
        f   = f0[:,:] #- ps[:,:]*constrKs[:,None]

        print( "pre-move ps", ps.flat.copy() )
        #print( "pre-move fs", f .flat.copy() )

        v  += f*dt    # move by external forces (ignoring constraints)   # NOTE: now we use steep descent, but we could verlet or other integrator of equations of motion
        ps += v*dt    # move by external forces (ignoring constraints)   # NOTE: now we use steep descent, but we could verlet or other integrator of equations of motion
        print( "post-move ps", ps.flat.copy() )
        #print( "pre-move vs", v .flat.copy() )

        if bUseCG:
            # K, fdl, ls, hdirs = makeMat( sticks, ps, l0s=l0s, constrKs=constrKs, kReg=kReg )      # fixed end points
            # f  = f0.flatten() + fdl
            # #dp = np.linalg.solve( K, f )    # solve (f0x+fdlx) = Kx*dx    aka b=A*x ( A=Kx, x=dx, b=f0x+fdlx )
            # dp, err2 = SolveCG( f, niter=10, eps=1e-6, bPrint=True, K=K )

            #print( "linearize ps", ps.flat.copy() )
            fdl, ls, hdirs = cs.linearize( sticks, ps, l0s=l0s, constrKs=constrKs, kReg=kReg )
            #print( "fdl", fdl )
            f = f0.flatten() + fdl
            print( "f(Ax=f)", f )
            hdirs_glob = hdirs
            dp, err2 = cs.SolveCG( f, niter=10, eps=1e-6, bPrint=False, dotFunc=dot_func )
            print( "dp ", dp )

            #exit()
        else:
            K, fdl, ls, hdirs = cs.makeMat( sticks, ps, l0s=l0s, constrKs=constrKs, kReg=kReg )      # fixed end points
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
os.system('mode con: cols=100 lines=50')
np.set_printoptions(linewidth=np.inf)


# =====  ROPE
''''
ps = np.array([     
    [-2.0, -0.0, 0.0],
    [-1.0, -0.0, 0.0],
    [ 0.0, -0.0, 0.0],
    [+1.0, -0.0, 0.0],
    [+2.0, -0.0, 0.0],
])
ms = np.array([ 1.0e+300, 1.0, 1.0, 1.0, 1.0e+300 ])
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
'''


# =====  Minimal-Truss
ps = np.array([     
    [ 0.0,  0.0, 0.0],
    [-1.0, -1.0, 0.0],
    [+1.0, -1.0, 0.0],
    [ 0.0, -2.0, 0.0],
])
ms = np.array([ 1.0e+300, 1.0, 1.0, 1.0, ])
k0 = 50.0
sticks =[
 ( 0,1, k0 ),
 ( 0,2, k0 ),
 ( 3,1, k0 ),
 ( 3,2, k0 ),
 ( 1,2, k0 ),
]
f0 = np.array([ 
[ 0.0, 0.0, 0.0],
[ 0.0, 0.0, 0.0],
[ 0.0, 0.0, 0.0],
[ 0.0,-5.0, 0.0],
])


# dynamics( ps, f0, constrKs=constrKs, kReg=kReg )

dynamics_PBD_GS( ps, ms, f0, sticks )



'''
#constrKs = [0.0, 0.0, 0.0, 0.0,0.0]
#kReg     = 0.0

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

#plt.legend()
plt.savefig( "try_Linearized_Move.png", bbox_inches='tight' )
plt.show()
