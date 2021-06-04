import os
import sys
import numpy as np

#sys.path.append("../")

sys.path.append('../')
from pyMeta import cpp_utils 
cpp_utils.clean_build    = False   # Recompile only if changed

import eFF
import CLCFGO as effmc
import eFF_terms as effpy

# ========== Globals

bDebug = 0
natom=0; norb=2; perOrb=1; nqOrb=1
bPrintInfo = False
label=""
plt  =None

rnd_pos  = 0
rnd_size = 0
rnd_coef = 0

# ========= Functions

def init_effmc( natom_=0, norb_=1, perOrb_=1, sz=0.5, dist=1.0 ):
    global natom,norb,perOrb,nqOrb
    natom=natom_; norb=norb_; perOrb=perOrb_; nqOrb=perOrb*(perOrb+1)/2
    effmc.init(natom,norb,perOrb,1)  #  natom, nOrb, perOrb, natypes
    ecoef = effmc.getBuff( "ecoef",(norb,perOrb)   )
    esize = effmc.getBuff( "esize",(norb,perOrb)   )
    epos  = effmc.getBuff( "epos" ,(norb,perOrb,3) )
    rhoQ  = effmc.getBuff( "rhoQ" ,(norb,nqOrb) )
    rhoS  = effmc.getBuff( "rhoS" ,(norb,nqOrb) )
    rhoP  = effmc.getBuff( "rhoP" ,(norb,nqOrb,3) )
    epos [:,:,:]= 0              + (np.random.rand(norb,perOrb,3)-0.5)*rnd_pos
    esize[:,:  ]=sz              + (np.random.rand(norb,perOrb  )-0.5)*rnd_size
    ecoef[:,:  ]=1               + (np.random.rand(norb,perOrb  )-0.5)*rnd_coef
    rhoP [:,:,:]=0               + (np.random.rand(norb,nqOrb,3 )-0.5)*rnd_pos
    rhoS [:,:  ]=sz*np.sqrt(0.5) + (np.random.rand(norb,nqOrb   )-0.5)*rnd_size
    rhoQ [:,:  ]=1               + (np.random.rand(norb,nqOrb   )-0.5)*rnd_coef
    if norb_>1:
        epos [1,0,0] += dist
        #ecoef[1,1]=0   # psi2=(1,0)
        #ecoef[1,0]=1   # psi2=(1,0)
    if perOrb_>1:    
        epos [:,1,0] = epos [:,0,0] + dist
    if(bPrintInfo): effmc.printAtomsAndElectrons()

def test_ProjectWf( Etoll=1e-5 ):
    print "\n ============ test_ProjectWf ( rho = |wf|^2 )"
    init_effmc( norb_=2, perOrb_=2, sz=0.5, dist=1.0 )
    # ---- test
    nps = 400+1
    ps  = np.zeros((nps,3))
    ps[:,0] = np.linspace( -5.0, 5.0, nps )
    effmc.eval() # we have to run it to project wavefuction to aux density
    wf  = effmc.orbAtPoints(ps)
    rho = effmc.rhoAtPoints(ps)
    wf2 = wf**2
    err = rho - wf2
    Err = np.sqrt( (err**2).sum()/len(err) )
    print " Error ", Err
    if plt is not None:
        plt.figure(figsize=(5,5))
        plt.plot( ps[:,0],    wf, label='wf'  )
        plt.plot( ps[:,0],   rho, label='rho' )
        plt.plot( ps[:,0],   wf2,':', label='wf^2')
        plt.legend(); plt.grid()
        plt.title( " test_ProjectWf( rho = |wf|^2 )" )
    #print "DEBUG 4 "
    return Err

def test_Poisson( Etoll=1e-5 ):
    print " ===== test_Poisson ( rho = dd_xyz V )"
    init_effmc( norb_=2, perOrb_=2, sz=0.5, dist=1.0 )
    effmc.eval() # we have to run it to project wavefuction to aux density
    dx=0.03; R=3.0
    err2, rho, rho_ =  effmc.test_Poisson( dx=dx, Rmax=R )
    Err = np.sqrt( err2/len(rho) )
    if(plt):
        xs=np.arange(0,R*2,dx)
        plt.figure(figsize=(5,5))
        plt.plot( xs, rho ,      label=('rho ' ) ); 
        plt.plot( xs, rho_, ":", label=('rho_' ) ); 
        plt.title( "test_Poisson( rho = dd_xyz V )" )
        plt.legend(); plt.grid()
    print " Error ", Err
    return Err

def test_OrbInteraction( Etoll=1e-5, iMODE=1 ):
    labels=[ "NONE", "Overlap Sij", "Kinetic Tij", "Coulomb Kij", ]
    label=labels[iMODE]
    print " ===== test_OrbInteraction "+label
    init_effmc( norb_=2, perOrb_=1, sz=0.5, dist=1.0 )

    effmc.eval() # we have to run it to project wavefuction to aux density
    dx=0.2; 
    #nint=30;
    nint=50
    err2, Ek, Ek_, f1, f2 =  effmc.test_OrbInteraction( iMODE=iMODE, io=0,jo=1, nint=nint, dx=dx, Rmax=5.0, bPrint=0, bSave=0  )
    #err2, rho, rho_ =  effmc.test_Poisson( dx=dx, Rmax=R, useWf=False )
    Err = np.sqrt( err2/len(Ek) )
    print "Error ", Err
    if(plt):
        xs=np.arange(0,dx*nint,dx)
        plt.figure(figsize=(10,5))
        plt.subplot(1,2,1)
        plt.plot( f1,      label=('f1' ) ); 
        plt.plot( f2, ":", label=('f2' ) );
        plt.title( "test_OrbInteraction funcs "+ label )
        plt.legend(); plt.grid()
        plt.subplot(1,2,2)
        plt.plot( xs, Ek ,    label=('Ek_ana' ) ); 
        plt.plot( xs, Ek_, ":", label=('Ek_num' ) ); 
        plt.title( "test_OrbInteraction "+label )
        #plt.ylim( 0, 20.0 )
        plt.legend(); plt.grid()
    return Err
def test_Overlap_Sij():
    return test_OrbInteraction(iMODE=1)
def test_Kinetic_Tij():
    return test_OrbInteraction(iMODE=2)
def test_Coulomb_Kij():
    effmc.setSwitches_( normalize=1, kinetic=-1, coulomb=1, exchange=-1, pauli=-1, AA=-1, AE=-1, AECoulomb=-1, AEPauli=-1 )
    return test_OrbInteraction(iMODE=3)

# ========= Check Forces

def processForces( xs,Es,fs ):
    n=len(xs)
    dx=xs[1]-xs[0]
    fs_num=(Es[2:]-Es[:-2])/(-2*dx)
    normF2 = (fs**2).sum()
    Err = np.sqrt( ( (fs_num-fs[1:-1])**2/normF2 ).sum()/(n-1) )
    print "Error ", Err
    if(plt):
        plt.figure(figsize=(5,5))
        plt.plot( xs,      Es    ,      label="E" )
        plt.plot( xs,      fs    ,      label="f_ana" )
        plt.plot( xs[1:-1],fs_num, ":", label="f_num" )
        #plt.plot( xs[1:-1],(fs_num-fs[1:-1])*10.0, label="(f_ana-f_num)*10.0" )
        plt.grid();plt.legend();
        plt.title(label)
    return Err

def checkForces( xname="ecoef", fname="efcoef", inds=(0,0), x0=0 ):
    nind = len(inds)
    szs  = (norb,perOrb,3)[:nind]
    xbuf = effmc.getBuff( xname,szs )
    fbuf = effmc.getBuff( fname,szs )
    x0  += x0_glob
    xs = np.arange(x0,x0+nx*dx,dx)
    Es = np.zeros(len(xs))
    fs = np.zeros(len(xs))
    #effmc.eval()  # ---- This is to make sure initial normalization is not a problem
    for i in range(nx):
        if(nind>2):
            xbuf[inds[0],inds[1],inds[2]] = xs[i]
        else:
            xbuf[inds[0],inds[1]]         = xs[i]
        Es[i] = effmc.eval()
        if(nind>2):
            fs[i] = fbuf[inds[0],inds[1],inds[2]]
        else:
            fs[i] = fbuf[inds[0],inds[1]]
    #print "Es ", Es
    #print "fs ", fs
    return processForces( xs,Es,fs )

def checkForces_Kinetic_epos( ):
    init_effmc( norb_=1, perOrb_=2, sz=0.75, dist=-0.1 )
    effmc.setSwitches_( normalize=-1, kinetic=1 )
    return checkForces( xname="epos",fname="efpos",inds=(0,0,0) )

def checkForces_Kinetic_esize( ):
    init_effmc( norb_=2, perOrb_=2, sz=0.5, dist=10.0 )
    effmc.setSwitches_( normalize=-1, kinetic=1 )
    return checkForces( xname="esize",fname="efsize",inds=(0,0), x0=0.25 )

def checkForces_Kinetic_ecoef( ):
    init_effmc( norb_=1, perOrb_=2, sz=0.75, dist=-0.1 )
    effmc.setSwitches_( normalize=-1, kinetic=1 )
    return checkForces( xname="ecoef",fname="efcoef",inds=(0,0) )

def checkForces_Hartree_epos( ):
    init_effmc( norb_=2, perOrb_=2, sz=0.2, dist=1.0 )
    effmc.setSwitches_( normalize=-1, coulomb=1 )
    return checkForces( xname="epos",fname="efpos",inds=(0,0,0) )

def checkForces_Hartree_esize( ):
    init_effmc( norb_=2, perOrb_=2, sz=0.75, dist=0.25 )
    effmc.setSwitches_( normalize=-1, coulomb=1 )
    return checkForces( xname="esize",fname="efsize",inds=(0,0) )

def checkForces_Hartree_ecoef( ):
    init_effmc( norb_=2, perOrb_=2, sz=0.75, dist=0.25 )
    effmc.setSwitches_( normalize=-1, coulomb=1 )
    return checkForces( xname="ecoef",fname="efcoef",inds=(0,0) )

# ========= Check Normalization derivatives

def check_dS( xname="ecoef", fname="enfcoef", inds=(0,0), x0=0 ):
    nind = len(inds)
    szs  = (norb,perOrb,3)[:nind]
    xbuf = effmc.getBuff( xname,szs )
    fbuf = effmc.getBuff( fname,szs )
    oQs  = effmc.getBuff( "oQs",(norb) )
    x0  += x0_glob
    xs   = np.arange(x0,x0+nx*dx,dx)
    Es   = np.zeros(len(xs))
    fs   = np.zeros(len(xs))
    effmc.eval()  # ---- This is to make sure initial normalization is not a problem
    for i in range(nx):
        if(nind>2):
            xbuf[inds[0],inds[1],inds[2]] = xs[i]
        else:
            xbuf[inds[0],inds[1]]         = xs[i]
        E     = effmc.eval()
        Es[i] = oQs[0]
        if(nind>2):
            fs[i] = fbuf[inds[0],inds[1],inds[2]]
        else:
            fs[i] = fbuf[inds[0],inds[1]]
    return processForces( xs,Es,fs )

def check_dS_epos( ):
    init_effmc( norb_=1, perOrb_=2, sz=1.0, dist=0.5 )
    effmc.setSwitches_( normalize=-1, normForce=1 )
    return check_dS( xname="epos", fname="enfpos", inds=(0,0,0) )

def check_dS_esize(  ):
    init_effmc( norb_=1, perOrb_=2, sz=1.0, dist=0.5 )
    effmc.setSwitches_( normalize=-1, normForce=1 )
    return check_dS( xname="esize", fname="enfsize", inds=(0,0), x0=0.25 )

def check_dS_ecoef( ):
    init_effmc( norb_=1, perOrb_=2, sz=1.0, dist=0.5 )
    effmc.setSwitches_( normalize=-1, normForce=1 )
    return check_dS( xname="ecoef", fname="enfcoef", inds=(0,0) )

# ========= Check Normalization derivatives

def check_Coulomb( xname="rhoQ", fname="rhofQ", inds=(0,0), x0=0 ):
    nind = len(inds)
    szs  = (norb,perOrb,3)[:nind]
    xbuf = effmc.getBuff( xname,szs )
    fbuf = effmc.getBuff( fname,szs )
    x0  += x0_glob
    xs = np.arange(x0,x0+nx*dx,dx)
    Es = np.zeros(len(xs))
    fs = np.zeros(len(xs))
    for i in range(nx):
        if(nind>2):
            xbuf[inds[0],inds[1],inds[2]] = xs[i]
        else:
            xbuf[inds[0],inds[1]]         = xs[i]
        Es[i] = effmc.coulombOrbPair( 0, 1 )
        if(nind>2):
            fs[i] = fbuf[inds[0],inds[1],inds[2]]
        else:
            fs[i] = fbuf[inds[0],inds[1]]
    return processForces( xs,Es,fs )

def check_Coulomb_rhoP_( ):
    init_effmc( norb_=2, perOrb_=1, sz=0.75, dist=-0.5 )
    return check_Coulomb( xname="rhoP", fname="rhofP", inds=(0,0,0) )

def check_Coulomb_rhoS_( ):
    init_effmc( norb_=2, perOrb_=1, sz=0.75, dist=-0.5 )
    return check_Coulomb( xname="rhoS", fname="rhofS", inds=(0,0), x0=0.5 )

def check_Coulomb_rhoQ_( ):
    init_effmc( norb_=2, perOrb_=1, sz=0.5, dist=-0.1 )
    return check_Coulomb( xname="rhoQ", fname="rhofQ", inds=(0,0) )

if __name__ == "__main__":
    import matplotlib.pyplot as plt_
    global plt,label
    global dx,nx,x0_glob
    global rnd_pos, rnd_size, rnd_coef
    x0_glob = 0.0001
    dx=0.05
    nx=40
    #nx=10
    #nx=2
    plt=plt_
    #bPrintInfo = True
    rnd_pos  = 0.2
    rnd_size = 0.2
    rnd_coef = 0.2
    '''
    xs = np.arange(0.0,1.0,0.25)
    ys = np.zeros(len(xs))
    for i,x in enumerate(xs):
        x = xs[i]
        s = 0.2
        E = effmc.evalFunc(xs[i],s)
        print(  x, s, "-> ", E )
    #exit(0)
    '''


    tests_funcs = []
    tests_funcs += [ test_ProjectWf, test_Poisson ]
    #tests_funcs += [ check_dS_epos,            check_dS_esize,              check_dS_ecoef             ]
    #tests_funcs += [ checkForces_Kinetic_epos, checkForces_Kinetic_esize ,  checkForces_Kinetic_ecoef  ]
    tests_funcs += [ checkForces_Hartree_epos, checkForces_Hartree_esize ,  checkForces_Hartree_ecoef  ]
    tests_funcs += [ check_Coulomb_rhoP_, check_Coulomb_rhoS_, check_Coulomb_rhoQ_ ]
    #tests_funcs += [ test_Overlap_Sij, test_Kinetic_Tij, test_Coulomb_Kij ]
    tests_results = []

    for test_func in tests_funcs:
        label = test_func.__name__
        effmc.setSwitches_( normalize=1, normForce=-1, kinetic=-1, coulomb=-1, exchange=-1, pauli=-1, AA=-1, AE=-1, AECoulomb=-1, AEPauli=-1 )
        tests_results.append( test_func() )
    print ""
    print " ##### Test Result Summary ##### "
    for i,test_func in enumerate(tests_funcs):
        print test_func.__name__," Error: ", tests_results[i]
    plt.show()
    print( " ===== ALL DONE !!! \n" )

