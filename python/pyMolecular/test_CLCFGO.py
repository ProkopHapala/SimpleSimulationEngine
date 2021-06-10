import os
import sys
import numpy as np

#sys.path.append("../")

sys.path.append('../')
from pyMeta import cpp_utils 
cpp_utils.clean_build    = False   # Recompile only if changed

#import eFF
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

def init_effmc_2x2( sz=0.5 ):
    global natom,norb,perOrb,nqOrb
    natom=0; norb=2; perOrb=2; nqOrb=perOrb*(perOrb+1)/2
    effmc.init(natom,norb,perOrb,1)  #  natom, nOrb, perOrb, natypes
    ecoef = effmc.getBuff( "ecoef",(norb,perOrb)   )
    esize = effmc.getBuff( "esize",(norb,perOrb)   )
    epos  = effmc.getBuff( "epos" ,(norb,perOrb,3) )
    epos [:,:,:]= 0
    #epos [0,0,0]=  0.0; epos [0,1,0]= -0.25; 
    #epos [1,0,1]= -0.5; epos [1,1,1]= -0.25     
    epos [0,0,0]= 0.0; epos [0,1,0]= 1.0; 
    epos [1,0,1]= 1.0;  epos [1,1,1]= 2.0   
    #esize[:,:]= 1
    #esize[:,:]= 0.5
    esize[0,0]= 1.8;    esize[0,1]= 1.1;   
    esize[1,0]= 1.5;    esize[1,1]= 0.3;     
    ecoef[:,:]= 1
    ecoef[0,0  ]= 0.4;  ecoef[0,1]= 0.3; 
    epos [:,:,:] += (np.random.rand(norb,perOrb,3)-0.5)*0.2 
    esize[:,:  ] += (np.random.rand(norb,perOrb  )-0.5)*0.2
    ecoef[:,:  ] += (np.random.rand(norb,perOrb  )-0.5)*0.2
    #ecoef[1,0  ]= 0.7;  ecoef[1,1]= 1.6               
    if(bPrintInfo): effmc.printAtomsAndElectrons()

def init_effmc( natom_=0, norb_=1, perOrb_=1, sz=0.5, dist=1.0, aQ=-1, aP=1.0, aPsz=0.1, aQsz=0.25 ):
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
    if natom_>0:
        apos   = effmc.getBuff( "apos"  ,(natom,3)  )
        aQs    = effmc.getBuff( "aQs"   ,(natom,) )
        aQsize = effmc.getBuff( "aQsize",(natom,) )
        aPcoef = effmc.getBuff( "aPcoef",(natom,) )
        aPsize = effmc.getBuff( "aPsize",(natom,) )
        apos[:,:] = 0
        apos[:,0] = np.arange(0,dist*natom_,dist)
        aQs   [:] = aQ
        aQsize[:] = aQsz
        aPcoef[:] = aP
        aPsize[:] = aPsz
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

def test_ETerms( xname="epos", inds=(0,0), x0=0 ):
    init_effmc( natom_=1, norb_=2, perOrb_=2, sz=0.25, dist=0.2 )
    effmc.setSwitches_( normalize=-1, normForce=-1, kinetic=1, coulomb=1, pauli=-1,    AA=1, AE=1, AECoulomb=1, AEPauli=1 )
    nind = len(inds)
    szs  = (norb,perOrb,3)[:nind]
    xbuf = effmc.getBuff( xname, szs )
    Ebuf = effmc.getEnergyTerms( )
    nterm = len(Ebuf)
    x0  += x0_glob
    xs = np.arange(x0,x0+nx*dx,dx)
    Etot = np.zeros(len(xs))
    Es   = np.zeros((len(xs),nterm))
    fs = np.zeros(len(xs))
    for i in range(nx):
        if(nind>2):
            xbuf[inds[0],inds[1],inds[2]] = xs[i]
        else:
            xbuf[inds[0],inds[1]]         = xs[i]
        Etot[i] = effmc.eval()
        #print "Ebuf", Ebuf
        Es[i,:] = Ebuf[:]
    if(plt):
        #print "test_ETerms PLOT"
        term_names = ["Ek","Eee","EeePaul","EeeExch","Eae","EaePaul","Eaa" ]
        plt.figure(figsize=(5,5))
        for i in range(nterm):
            plt.plot( xs, Es[:,i], label=term_names[i] )
        plt.plot    ( xs, Etot, "-k", lw=2, label="Etot"  )
        plt.grid();plt.legend();
        plt.title(label)
    return xs, Etot, Es

# ========= Check Forces

def processForces( xs,Es,fs ):
    n=len(xs)
    dx=xs[1]-xs[0]
    fs_num=(Es[2:]-Es[:-2])/(-2*dx)
    normF2 = (fs**2).sum()
    Err = np.sqrt( ( (fs_num-fs[1:-1])**2/normF2 ).sum()/(n-1) )
    print label, "Error ", Err
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
    if(nind==1):
        szs  = (natom,)
    else:    
        szs  = (norb,perOrb,3)[:nind]
    xbuf = effmc.getBuff( xname,szs )
    fbuf = effmc.getBuff( fname,szs )
    x0  += x0_glob
    xs = np.arange(x0,x0+nx*dx,dx)
    Es = np.zeros(len(xs))
    fs = np.zeros(len(xs))
    #effmc.eval()  # ---- This is to make sure initial normalization is not a problem
    for i in range(nx):
        xbuf[inds]= xs[i]
        Es[i] = effmc.eval()
        fs[i] = fbuf[inds]
    #print "Es ", Es
    #print "fs ", fs
    return processForces( xs,Es,fs )

# ========= Kinetic

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

# ========= Hartree

def checkForces_Hartree_epos( ):
    #init_effmc( norb_=2, perOrb_=2, sz=0.2, dist=1.0 )
    init_effmc_2x2( sz=0.5 )
    effmc.setSwitches_( normalize=-1, coulomb=1 )
    return checkForces( xname="epos",fname="efpos",inds=(0,0,0) )

def checkForces_Hartree_esize( ):
    #init_effmc( norb_=2, perOrb_=2, sz=0.75, dist=0.25 )
    init_effmc_2x2( sz=0.5 )
    effmc.setSwitches_( normalize=-1, coulomb=1 )
    return checkForces( xname="esize",fname="efsize",inds=(0,0) )

def checkForces_Hartree_ecoef( ):
    #init_effmc( norb_=2, perOrb_=2, sz=0.75, dist=0.25 )
    init_effmc_2x2( sz=0.5 )
    effmc.setSwitches_( normalize=-1, coulomb=1 )
    return checkForces( xname="ecoef",fname="efcoef",inds=(0,0) )

# ========= Pauli

def checkForces_Pauli_epos( ):
    init_effmc( norb_=2, perOrb_=2, sz=0.5, dist=1.0 )
    effmc.setSwitches_( normalize=-1, pauli=1 )
    return checkForces( xname="epos",fname="efpos",inds=(0,0,0) )

def checkForces_Pauli_esize( ):
    init_effmc( norb_=2, perOrb_=2, sz=0.5, dist=1.0 )
    effmc.setSwitches_( normalize=-1, pauli=1 )
    return checkForces( xname="esize",fname="efsize",inds=(0,0) )

def checkForces_Pauli_ecoef( ):
    init_effmc( norb_=2, perOrb_=2, sz=0.5, dist=1.0 )
    effmc.setSwitches_( normalize=-1, pauli=1 )
    return checkForces( xname="ecoef",fname="efcoef",inds=(0,0) )

# ========= Atom-Electron Coulomb

def checkForces_AQ_epos( ):
    init_effmc( natom_=1, norb_=1, perOrb_=1, sz=0.5, dist=1.0 )
    effmc.setSwitches_( normalize=-1,   AE=1, AECoulomb=1, AEPauli=-1 )
    return checkForces( xname="epos",fname="efpos",inds=(0,0,0) )

def checkForces_AQ_esize( ):
    init_effmc( natom_=1, norb_=1, perOrb_=1, sz=0.5, dist=1.0 )
    effmc.setSwitches_( normalize=-1,  AE=1, AECoulomb=1, AEPauli=-1 )
    return checkForces( xname="esize",fname="efsize",inds=(0,0) )

def checkForces_AQ_ecoef( ):
    init_effmc( natom_=1, norb_=1, perOrb_=1, sz=0.5, dist=1.0 )
    effmc.setSwitches_( normalize=-1,  AE=1, AECoulomb=1, AEPauli=-1 )
    return checkForces( xname="ecoef",fname="efcoef",inds=(0,0) )

# ========= Atom-Electron Pauli

def checkForces_AP_epos( ):
    init_effmc( natom_=1, norb_=1, perOrb_=1, sz=0.5, dist=1.0 )
    effmc.setSwitches_( normalize=-1,   AE=1, AECoulomb=-1, AEPauli=+1 )
    return checkForces( xname="epos",fname="efpos",inds=(0,0,0) )

def checkForces_AP_esize( ):
    init_effmc( natom_=1, norb_=1, perOrb_=1, sz=0.5, dist=1.0 )
    effmc.setSwitches_( normalize=-1,  AE=1, AECoulomb=-1, AEPauli=+1 )
    return checkForces( xname="esize",fname="efsize",inds=(0,0) )

def checkForces_AP_ecoef( ):
    init_effmc( natom_=1, norb_=1, perOrb_=1, sz=0.5, dist=1.0 )
    effmc.setSwitches_( normalize=-1,  AE=1, AECoulomb=-1, AEPauli=+1 )
    return checkForces( xname="ecoef",fname="efcoef",inds=(0,0) )

# ========= Atom-Atom Electrostatics

def checkForces_AA_pos( ):
    init_effmc( natom_=1, norb_=0, perOrb_=0, sz=0.5, dist=1.0 )
    effmc.setSwitches_( normalize=-1,   AE=1, AECoulomb=+1, AEPauli=+1 )
    return checkForces( xname="apos",fname="aforce",inds=(0,) )

# ========= Total

def checkForces_Tot_epos( ):
    init_effmc( natom_=1, norb_=2, perOrb_=2, sz=0.5, dist=1.0 )
    #effmc.setSwitches_( normalize=-1, normForce=-1, kinetic=1, coulomb=1, pauli=-1     )
    effmc.setSwitches_( normalize=-1, normForce=-1, kinetic=1, coulomb=1, pauli=-1,    AA=1, AE=1, AECoulomb=1, AEPauli=1 )
    return checkForces( xname="epos",fname="efpos",inds=(0,0,0) )

def checkForces_Tot_esize( ):
    init_effmc( natom_=1, norb_=2, perOrb_=2, sz=0.5, dist=1.0 )
    #effmc.setSwitches_( normalize=-1, normForce=-1, kinetic=1, coulomb=1, pauli=-1     )
    effmc.setSwitches_( normalize=-1, normForce=-1, kinetic=1, coulomb=1, pauli=-1,    AA=1, AE=1, AECoulomb=1, AEPauli=1 )
    return checkForces( xname="esize",fname="efsize",inds=(0,0) )

def checkForces_Tot_ecoef( ):
    init_effmc( natom_=1, norb_=2, perOrb_=2, sz=0.5, dist=1.0 )
    #effmc.setSwitches_( normalize=-1, normForce=-1, kinetic=1, coulomb=1, pauli=-1     )
    effmc.setSwitches_( normalize=-1, normForce=-1, kinetic=1, coulomb=1, pauli=-1,    AA=1, AE=1, AECoulomb=1, AEPauli=1 )
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
        xbuf[inds] = xs[i]
        E     = effmc.eval()
        Es[i] = oQs[0]
        fs[i] = fbuf[inds]
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
        xbuf[inds] = xs[i]
        Es[i] = effmc.coulombOrbPair( 0, 1 )
        fs[i] = fbuf[inds]
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

    #np.random.seed( 451)
    #np.random.seed( 1)   # good
    np.random.seed( 3)
    #np.random.seed( 5)    # exact
    #np.random.seed( 8)
    import matplotlib.pyplot as plt_
    global plt,label
    global dx,nx,x0_glob
    global rnd_pos, rnd_size, rnd_coef
    x0_glob = 0.25001
    dx=0.0125
    nx=150
    #nx=50
    #nx=10
    #nx=2
    plt=plt_
    bPrintInfo = True
    #rnd_pos  = 0.0; rnd_size = 0.0; rnd_coef = 0.0
    #rnd_pos  = 0.2; rnd_size = 0.2; rnd_coef = 0.2
    rnd_pos  = 0.2; rnd_size = 0.2; rnd_coef = 0.5
    #rnd_pos  = 1.0; rnd_size = 0.2; rnd_coef = 0.2
    
    #xs = np.arange(0.0,6.0,0.025)
    xs = np.arange(0.1,3.0,0.005)
    ys = np.zeros(len(xs))

    #Es,Fs = effmc.evalFuncDerivs(xs,0.6)
    Es,Fs = effmc.evalFuncDerivs(0.165,xs)
    plt.plot(xs,Es,label="E")
    plt.plot(xs,Fs,label="F")
    plt.plot(xs[1:-1],(Es[2:]-Es[:-2])/(-2*(xs[1]-xs[0])),':',label="F_num")
    plt.legend(); plt.grid(); plt.show(); exit(0)

    '''
    for i,x in enumerate(xs):
        x = xs[i]
        s = 0.2
        E = effmc.evalFunc(xs[i],s)
        print(  x, s, "-> ", E )
    exit(0)
    '''
    #effmc.setPauliMode(0)  # E = K*S^2
    effmc.setPauliMode(2)  # E = Sij^2/(1-Sij^2) * ( Tii + Tjj - 2Tij/Sij )
    #effmc.setPauliMode(3)  # E=T
    #effmc.setPauliMode(4)  # E=S
    #effmc.setPauliMode(5)   # Ep = ( Sij/(1-Sij^2) )* Tij 
    #effmc.setPauliMode(6)   # Ep = Sij*Tij

    '''
    test_ETerms( xname="epos", inds=(0,0), x0=0 )
    plt.show()    
    exit(0)
    '''

    tests_results = []
    tests_funcs = []
    #tests_funcs += [ test_ProjectWf, test_Poisson ]
    #tests_funcs += [ check_dS_epos,            check_dS_esize,              check_dS_ecoef             ]
    #tests_funcs += [ checkForces_Kinetic_epos, checkForces_Kinetic_esize ,  checkForces_Kinetic_ecoef  ]
    #tests_funcs += [ checkForces_Pauli_epos ]
    #tests_funcs += [ checkForces_Pauli_epos, checkForces_Pauli_esize ,  checkForces_Pauli_ecoef  ]
    #tests_funcs += [ checkForces_Hartree_epos, checkForces_Hartree_esize ,  checkForces_Hartree_ecoef  ]
    #tests_funcs += [ checkForces_AQ_epos, checkForces_AQ_esize ,  checkForces_AQ_ecoef  ]
    #tests_funcs += [ checkForces_AP_epos, checkForces_AP_esize ,  checkForces_AP_ecoef  ]
    #tests_funcs += [ checkForces_AA_pos ]
    tests_funcs += [ checkForces_Tot_epos, checkForces_Tot_esize ,  checkForces_Tot_ecoef  ]
    #tests_funcs += [ check_Coulomb_rhoP_, check_Coulomb_rhoS_, check_Coulomb_rhoQ_ ]
    #tests_funcs += [ test_Overlap_Sij, test_Kinetic_Tij, test_Coulomb_Kij ]
    
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

