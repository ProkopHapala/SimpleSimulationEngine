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

iNorm = -1
bDebug = 0
natom=0; norb=2; perOrb=1; nqOrb=1
bPrintInfo = False
label=""
plt  =None

ss_glob=None; xs_glob=None; cs_glob=None

rnd_pos  = 0
rnd_size = 0
rnd_coef = 0

fnnScale = 10

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

def init_effmc_H2mol_1g( adist=1.0, edist=0.1, sz=0.9, aQ=-1,aQsz=0.0, aP=0.0,aPsz=0.1 ):
    global natom,norb,perOrb,nqOrb
    natom=2; norb=2; perOrb=1; nqOrb=perOrb*(perOrb+1)/2
    effmc.init(natom,norb,perOrb,1)  #  natom, nOrb, perOrb, natypes
    ecoef = effmc.getBuff( "ecoef",(norb,perOrb)   )
    esize = effmc.getBuff( "esize",(norb,perOrb)   )
    epos  = effmc.getBuff( "epos" ,(norb,perOrb,3) )
    ospin = effmc.getIBuff( "ospin",(norb)  )
    apos  = effmc.getBuff( "apos" ,(natom,3)  )
    aPars = effmc.getBuff( "aPars",(natom,4)  )
    ospin[0]=1; ospin[1]=-1
    epos [:,:,:]= 0;  epos[0,:,0]=-edist/2;  epos[1,:,0]=+edist/2; 
    esize[:,:]= sz   
    ecoef[:,:]= 1
    #epos [:,:,:] += (np.random.rand(norb,perOrb,3)-0.5)*0.2 
    #esize[:,:  ] += (np.random.rand(norb,perOrb  )-0.5)*0.2
    #ecoef[:,:  ] += (np.random.rand(norb,perOrb  )-0.5)*0.2    
    apos    [:,:] = 0;  apos[0,0] = -adist/2; apos[1,0] = +adist/2   
    aPars   [:,0] = aQ
    aPars   [:,1] = aQsz
    aPars   [:,2] = aPsz
    aPars   [:,3] = aP      
    if(bPrintInfo): effmc.printAtomsAndElectrons()

def init_effmc( natom_=0, norb_=1, perOrb_=1, sz=0.5, dist=1.0, aQ=-1,aQsz=0.0, aP=0.0,aPsz=0.1  ):
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
        aPars  = effmc.getBuff( "aPars",(natom,4)  )
        #aQs    = effmc.getBuff( "aQs"   ,(natom,) )
        #aQsize = effmc.getBuff( "aQsize",(natom,) )
        #aPcoef = effmc.getBuff( "aPcoef",(natom,) )
        #aPsize = effmc.getBuff( "aPsize",(natom,) )
        apos[:,:] = 0
        apos[:,0] = np.arange(0,dist*natom_,dist)
        #aQs   [:] = aQ
        #aQsize[:] = aQsz
        #aPcoef[:] = aP
        #aPsize[:] = aPsz
        aPars   [:,0] = aQ
        aPars   [:,1] = aQsz
        aPars   [:,2] = aPsz
        aPars   [:,3] = aP
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
        plt.plot( xs, Ek ,      label=('Ek_ana' ) ); 
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

def plot_Terms( xs=None, xname="epos", inds=(0,0) ):
    nind = len(inds)
    szs  = (norb,perOrb,3)[:nind]
    xbuf = effmc.getBuff( xname, szs )
    Ebuf = effmc.getEnergyTerms( )
    nterm = len(Ebuf)
    Etot  = np.zeros(len(xs))
    Es    = np.zeros((len(xs),nterm))
    #fs    = np.zeros(len(xs))
    for i in range(len(xs)):
        if(nind>2):
            xbuf[inds[0],inds[1],inds[2]] = xs[i]
        else:
            xbuf[inds[0],inds[1]]         = xs[i]
        Etot[i] = effmc.eval()
        Es[i,:] = Ebuf[:]
    if(plt):
        #print "test_ETerms PLOT"
        term_names = ["Ek","Eee","EeePaul","EeeExch","Eae","EaePaul","Eaa" ]
        plt.figure(figsize=(5,5))
        for i in range(nterm):
            plt.plot( xs, Es[:,i], label=term_names[i] )
        plt.plot    ( xs, Etot, "k", lw=4, label="Etot"  )
        plt.grid();plt.legend();
        plt.title(label)    
    return xs, Etot, Es

def test_ETerms( xname="epos", inds=(0,0), xs=xs_glob ):
    init_effmc( natom_=1, norb_=2, perOrb_=2, sz=0.25, dist=0.2 )
    effmc.setSwitches_( normalize=-1, normForce=-1, kinetic=1, coulomb=1, pauli=-1,    AA=1, AE=1, AECoulomb=1, AEPauli=1 )
    return plot_Terms( xname=xname, inds=inds, xs=xs ) 

# ========= Check Forces

def processForces( xs,Es,fs, fns=None, es=None ):
    n=len(xs)
    dx=xs[1]-xs[0]
    fs_num=(Es[2:]-Es[:-2])/(2*dx)
    normF2 = (fs**2).sum()
    Err = np.sqrt( ( (fs_num-fs[1:-1])**2/normF2 ).sum()/(n-1) )
    print label, "Error ", Err
    if(plt):
        plt.figure(figsize=(5,5))
        plt.plot( xs,      Es    ,      label="E" )
        plt.plot( xs,     -fs    ,      label="f_ana" )
        plt.plot( xs[1:-1],fs_num, ":", label="f_num" )
        if fns is not None:
            plt.plot( xs,fns*fnnScale, ":", label="fns" )
        if fns is not None:
            plt.plot( xs,es, ":", label="es" )
        #plt.plot( xs[1:-1],(fs_num-fs[1:-1])*10.0, label="(f_ana-f_num)*10.0" )
        plt.grid();plt.legend();
        plt.title(label)
    return Err

def checkForces( xname="ecoef", fname="efcoef", inds=(0,0), xs=None ):
    nind = len(inds)
    if(nind==1):
        szs  = (natom,)
    else:    
        szs  = (norb,perOrb,3)[:nind]
    xbuf = effmc.getBuff( xname,szs )
    fbuf = effmc.getBuff( fname,szs )
    #if xs is None:
    #    x0  += x0_glob
    #    xs = np.arange(x0,x0+nx*dx,dx)
    Es = np.zeros(len(xs))
    fs = np.zeros(len(xs))
    #effmc.eval()  # ---- This is to make sure initial normalization is not a problem
    for i in range(len(xs)):
        xbuf[inds]= xs[i]
        Es[i] = effmc.eval()
        fs[i] = fbuf[inds]
    #print "Es ", Es
    #print "fs ", fs
    return processForces( xs,Es,fs )

def checkForces_norm( xname="ecoef", fname="efcoef", fnname="enfcoef", inds=(0,0), xs=None ):
    nind = len(inds)
    print norb,perOrb, " in checkForces_norm "
    if(nind==1):
        szs  = (natom,)
    else:    
        szs  = (norb,perOrb,3)[:nind]
    xbuf   = effmc.getBuff( xname,szs )
    fbuf   = effmc.getBuff( fname,szs )
    fnbuf  = effmc.getBuff( fnname,szs )
    epsbuf = effmc.getBuff( "oEs",szs[0] )
    #if xs is None:
    #    x0  += x0_glob
    #    xs = np.arange(x0,x0+nx*dx,dx)
    Es  = np.zeros(len(xs))
    fs  = np.zeros(len(xs))
    fns = np.zeros(len(xs))
    es  = np.zeros(len(xs))
    #effmc.eval()  # ---- This is to make sure initial normalization is not a problem
    print "xbuf.shape", xbuf.shape, szs, norb, perOrb
    for i in range(len(xs)):
        if xname=="ecoef":
            xbuf[(0,1)] = 1
        xbuf[inds]= xs[i]
        Es[i]  = effmc.eval()
        fs[i]  = fbuf[inds]
        fns[i] = fnbuf[inds]
        es[i]  = epsbuf[inds[0]]
    #print "Es ", Es
    #print "fs ", fs
    return processForces( xs,Es,fs, fns=fns, es=es )

# ========= Kinetic

def checkForces_Kinetic_epos( ):
    init_effmc( norb_=1, perOrb_=2, sz=0.75, dist=-0.1 )
    effmc.setSwitches_( normalize=iNorm, normForce=iNorm, kinetic=1 )
    return checkForces( xname="epos",fname="efpos",inds=(0,0,0), xs=xs_glob )

def checkForces_Kinetic_esize( ):
    init_effmc( norb_=2, perOrb_=2, sz=0.5, dist=10.0 )
    effmc.setSwitches_( normalize=iNorm, normForce=iNorm, kinetic=1 )
    return checkForces( xname="esize",fname="efsize",inds=(0,0), xs=ss_glob )

def checkForces_Kinetic_ecoef( ):
    init_effmc( norb_=1, perOrb_=2, sz=0.75, dist=-0.1 )
    effmc.setSwitches_( normalize=iNorm, normForce=iNorm, kinetic=1 )
    return checkForces( xname="ecoef",fname="efcoef",inds=(0,0), xs=cs_glob )

# ========= Hartree

def checkForces_Hartree_epos( ):
    #init_effmc( norb_=2, perOrb_=2, sz=0.2, dist=1.0 )
    init_effmc_2x2( sz=0.5 )
    effmc.setSwitches_( normalize=iNorm, normForce=iNorm, coulomb=1 )
    return checkForces( xname="epos",fname="efpos",inds=(0,0,0), xs=xs_glob )

def checkForces_Hartree_esize( ):
    #init_effmc( norb_=2, perOrb_=2, sz=0.75, dist=0.25 )
    init_effmc_2x2( sz=0.5 )
    effmc.setSwitches_( normalize=iNorm, normForce=iNorm, coulomb=1 )
    return checkForces( xname="esize",fname="efsize",inds=(0,0), xs=ss_glob )

def checkForces_Hartree_ecoef( ):
    #init_effmc( norb_=2, perOrb_=2, sz=0.75, dist=0.25 )
    init_effmc_2x2( sz=0.5 )
    effmc.setSwitches_( normalize=iNorm, normForce=iNorm, coulomb=1 )
    return checkForces( xname="ecoef",fname="efcoef",inds=(0,0), xs=cs_glob )

# ========= Pauli

def checkForces_Pauli_epos( ):
    init_effmc( norb_=2, perOrb_=2, sz=0.5, dist=1.0 )
    effmc.setSwitches_( normalize=iNorm, normForce=iNorm, pauli=1 )
    return checkForces( xname="epos",fname="efpos",inds=(0,0,0), xs=xs_glob )

def checkForces_Pauli_esize( ):
    init_effmc( norb_=2, perOrb_=2, sz=0.5, dist=1.0 )
    effmc.setSwitches_( normalize=iNorm, normForce=iNorm, pauli=1 )
    return checkForces( xname="esize",fname="efsize",inds=(0,0), xs=ss_glob )

def checkForces_Pauli_ecoef( ):
    init_effmc( norb_=2, perOrb_=2, sz=0.5, dist=1.0 )
    effmc.setSwitches_( normalize=iNorm, normForce=iNorm, pauli=1 )
    return checkForces( xname="ecoef",fname="efcoef",inds=(0,0), xs=cs_glob )

# ========= Atom-Electron Coulomb

def checkForces_AQ_epos( ):
    init_effmc( natom_=1, norb_=1, perOrb_=2, sz=0.5, dist=1.0 )
    effmc.setSwitches_( normalize=iNorm, normForce=iNorm,   AE=1, AECoulomb=1, AEPauli=-1 )
    return checkForces( xname="epos",fname="efpos",inds=(0,0,0), xs=xs_glob )

def checkForces_AQ_esize( ):
    init_effmc( natom_=1, norb_=1, perOrb_=2, sz=0.5, dist=1.0 )
    effmc.setSwitches_( normalize=iNorm, normForce=iNorm,  AE=1, AECoulomb=1, AEPauli=-1 )
    return checkForces( xname="esize",fname="efsize",inds=(0,0), xs=ss_glob )

def checkForces_AQ_ecoef( ):
    init_effmc( natom_=1, norb_=1, perOrb_=2, sz=0.5, dist=1.0 )
    effmc.setSwitches_( normalize=iNorm, normForce=iNorm,  AE=1, AECoulomb=1, AEPauli=-1 )
    return checkForces( xname="ecoef",fname="efcoef",inds=(0,0), xs=cs_glob )

# ========= Atom-Electron Pauli

def checkForces_AP_epos( ):
    init_effmc( natom_=1, norb_=1, perOrb_=2, sz=0.5, dist=1.0 )
    effmc.setSwitches_( normalize=iNorm, normForce=iNorm,   AE=1, AECoulomb=-1, AEPauli=+1 )
    return checkForces( xname="epos",fname="efpos",inds=(0,0,0), xs=xs_glob )

def checkForces_AP_esize( ):
    init_effmc( natom_=1, norb_=1, perOrb_=2, sz=0.5, dist=1.0 )
    effmc.setSwitches_( normalize=iNorm, normForce=iNorm,  AE=1, AECoulomb=-1, AEPauli=+1 )
    return checkForces( xname="esize",fname="efsize",inds=(0,0), xs=ss_glob )

def checkForces_AP_ecoef( ):
    init_effmc( natom_=1, norb_=1, perOrb_=2, sz=0.5, dist=1.0 )
    effmc.setSwitches_( normalize=iNorm, normForce=iNorm,  AE=1, AECoulomb=-1, AEPauli=+1 )
    return checkForces( xname="ecoef",fname="efcoef",inds=(0,0), xs=cs_glob )

# ========= Atom-Atom Electrostatics

def checkForces_AA_pos( ):
    init_effmc( natom_=1, norb_=0, perOrb_=0, sz=0.5, dist=1.0 )
    effmc.setSwitches_( normalize=-1,   AE=1, AECoulomb=+1, AEPauli=+1 )
    return checkForces( xname="apos",fname="aforce",inds=(0,), xs=xs_glob )

# ========= Total

def checkForces_Tot_epos( ):
    init_effmc( natom_=1, norb_=2, perOrb_=2, sz=0.5, dist=1.0 )
    effmc.setSwitches_( normalize=iNorm, normForce=iNorm, kinetic=1, coulomb=1, pauli=-1,    AA=1, AE=1, AECoulomb=1, AEPauli=1 )
    return checkForces( xname="epos",fname="efpos",inds=(0,0,0), xs=xs_glob )

def checkForces_Tot_esize( ):
    init_effmc( natom_=1, norb_=2, perOrb_=2, sz=0.5, dist=1.0 )
    effmc.setSwitches_( normalize=iNorm, normForce=iNorm, kinetic=1, coulomb=1, pauli=-1,    AA=1, AE=1, AECoulomb=1, AEPauli=1 )
    return checkForces( xname="esize",fname="efsize",inds=(0,0), xs=ss_glob )

def checkForces_Tot_ecoef( ):
    init_effmc( natom_=1, norb_=2, perOrb_=2, sz=0.5, dist=1.0 )
    effmc.setSwitches_( normalize=iNorm, normForce=iNorm, kinetic=1, coulomb=1, pauli=-1,    AA=1, AE=1, AECoulomb=1, AEPauli=1 )
    return checkForces( xname="ecoef",fname="efcoef",inds=(0,0), xs=cs_glob )

# ========= Check Normalization derivatives

def check_dS( xname="ecoef", fname="enfcoef", inds=(0,0), xs=None ):
    nind = len(inds)
    szs  = (norb,perOrb,3)[:nind]
    xbuf = effmc.getBuff( xname,szs )
    fbuf = effmc.getBuff( fname,szs )
    oQs  = effmc.getBuff( "oQs",(norb) )
    #x0  += x0_glob
    #xs   = np.arange(x0,x0+nx*dx,dx)
    #xs=xs_glob
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
    return check_dS( xname="epos", fname="enfpos", inds=(0,0,0), xs=xs_glob )

def check_dS_esize(  ):
    init_effmc( norb_=1, perOrb_=2, sz=1.0, dist=0.5 )
    effmc.setSwitches_( normalize=-1, normForce=1 )
    return check_dS( xname="esize", fname="enfsize", inds=(0,0), xs=ss_glob )

def check_dS_ecoef( ):
    init_effmc( norb_=1, perOrb_=2, sz=1.0, dist=0.5 )
    effmc.setSwitches_( normalize=-1, normForce=1 )
    return check_dS( xname="ecoef", fname="enfcoef", inds=(0,0), xs=cs_glob )

# ========= Check Normalization derivatives

def check_Coulomb( xname="rhoQ", fname="rhofQ", inds=(0,0), xs=None ):
    nind = len(inds)
    szs  = (norb,perOrb,3)[:nind]
    xbuf = effmc.getBuff( xname,szs )
    fbuf = effmc.getBuff( fname,szs )
    #x0  += x0_glob
    #xs = np.arange(x0,x0+nx*dx,dx)
    #print " x0  ", x0
    Es = np.zeros(len(xs))
    fs = np.zeros(len(xs))
    for i in range(nx):
        xbuf[inds] = xs[i]
        Es[i] = effmc.coulombOrbPair( 0, 1 )
        fs[i] = fbuf[inds]
    return processForces( xs,Es,fs )

def check_Coulomb_rhoP_( ):
    init_effmc( norb_=2, perOrb_=1, sz=0.8, dist=-0.0 )
    return check_Coulomb( xname="rhoP", fname="rhofP", inds=(0,0,0), xs=xs_glob )

def check_Coulomb_rhoS_( ):
    init_effmc( norb_=2, perOrb_=1, sz=0.8, dist=-0.5 )
    return check_Coulomb( xname="rhoS", fname="rhofS", inds=(0,0), xs=ss_glob )

def check_Coulomb_rhoQ_( ):
    init_effmc( norb_=2, perOrb_=1, sz=0.8, dist=-0.1 )
    return check_Coulomb( xname="rhoQ", fname="rhofQ", inds=(0,0), xs=cs_glob )

def test_EvalFuncDerivs( r0=None, s0=None, label="test_EvalFuncDerivs"):
    # - Note check firts what is currently in evalFuncDerivs()
    # ====== test Gauss::Coulomb()
    if r0 is None:
        xs = np.arange( 0.001, 5.0, 0.05 )
        Es,Fs = effmc.evalFuncDerivs(xs,s0)
    else:
        xs = np.arange( 0.3, 3.0, 0.005 )
        Es,Fs = effmc.evalFuncDerivs(r0,xs)
    plt.plot(xs,Es,label="E")
    plt.plot(xs,Fs,label="F")
    plt.plot(xs[1:-1],(Es[2:]-Es[:-2])/(-2*(xs[1]-xs[0])),':',label="F_num")
    plt.legend(); plt.grid(); 
    plt.title(label)
    #plt.show(); 

def test_Hatom():
    # ====== test H-atom
    # - Note check first what is currently in evalFuncDerivs()
    init_effmc( natom_=1, norb_=1, perOrb_=1, sz=0.5, dist=1.0, aQ=-1,aQsz=0.0, aP=0.0,aPsz=0.0 )
    effmc.setSwitches_( normalize=-1, normForce=-1, kinetic=1, coulomb=1, pauli=1,    AA=1, AE=1, AECoulomb=1, AEPauli=-1 )
    xs = np.arange(0.3,3.0,0.05)
    plot_Terms( xs=xs, xname="esize" )

    import eFF_terms as eff
    Ek,Eae = eff.Hatom(xs)
    plt.plot(xs, Ek, ':',     label="Ek_eFF")
    plt.plot(xs, Eae,':',     label="Eae_eFF")
    plt.plot(xs, Ek+Eae, ':', label="Etot_eFF")
    plt.legend()

def checkForces_H_2g( inds=(0,0) ):
    global label
    #effmc.setSwitches_( normalize=-1, normForce=-1, kinetic=1,  AE=1, AECoulomb=1 )
    #effmc.setSwitches_( normalize=1, normForce=-1, kinetic=1, AE=-1, AECoulomb=-1 )
    #effmc.setSwitches_( normalize=1, normForce=1, kinetic=1,  AE=1, AECoulomb=1 )
    effmc.setSwitches_( normalize=1, normForce=+1, kinetic=1, AE=-1, AECoulomb=-1 )
    #effmc.loadFromFile( "../../cpp/sketches_SDL/Molecular/data/H_2g_problem_sym.fgo"  )
    effmc.loadFromFile( "../../cpp/sketches_SDL/Molecular/data/H_2g_compare.fgo"  )
    label="epos"; checkForces( xname="epos",  fname="efpos",  inds=inds, xs=np.arange(-1.0,1.0,0.01) )
    #label="epos"; checkForces( xname="epos",  fname="efpos",  inds=inds, xs=np.arange( 0.1,1.0,0.01) )
    label="esize"; checkForces( xname="esize", fname="efsize", inds=inds, xs=np.arange(1.0,2.0,0.01) )
    label="ecoef"; checkForces( xname="ecoef", fname="efcoef", inds=inds, xs=np.arange(5.0,7.0,0.01) )

def compareForces_H_2g( inds=(0,0), bNormalize=1 ):
    import CLCFGO_normalization_derivs_2 as effpy
    global label
    
    #bNormalize = -1
    bNormalize = +1

    #effmc.setSwitches_( normalize=bNormalize, normForce=bNormalize, kinetic=+1,  AE=-1, AECoulomb=+1, AEPauli=-1 )    # Kinetic
    effmc.setSwitches_( normalize=bNormalize, normForce=bNormalize, kinetic=-1,  AE=+1, AECoulomb=+1, AEPauli=-1 )    # AECoulomb
    #effmc.setSwitches_( normalize=bNormalize, normForce=bNormalize, kinetic=-1,  AE=+1, AECoulomb=-1, AEPauli=-1 )    # AEPauli
    #effmc.setSwitches_( normalize=bNormalize, normForce=bNormalize, kinetic=+1,  AE=+1, AECoulomb=+1, AEPauli=+1 )    # Total

    #effmc.loadFromFile( "../../cpp/sketches_SDL/Molecular/data/H_2g-.fgo"  )
    effmc.loadFromFile( "../../cpp/sketches_SDL/Molecular/data/H_2g_compare.fgo"  )
    effmc.printSetup()
    effmc.printAtomsAndElectrons()
    global norb,perOrb
    # natypes,natom,nOrb,nBas,perOrb,perOrb2,nqOrb,nQtot
    arr = effmc.getDimPointer();    
    norb   = arr[2]
    perOrb = arr[4]
    print norb,perOrb, arr

    #Q,E, (dEdxa,dEdsa,dEdca),(dQdxa,dQdsa,dQdca),xs = effpy.evalTest( what="xa", bNormalize=bNormalize,    xa=-0.4,sa=0.35,ca=1.6,     xb=+0.5,sb=0.55,cb=-0.4 )
    Q,E, (dEdxa,dEdsa,dEdca),(dQdxa,dQdsa,dQdca),xs = effpy.evalTest( what="ca", bNormalize=(bNormalize>0),     xa=-0.4,sa=0.35,ca=1.6,     xb=+0.5,sb=0.55,cb=-0.4 )

    #label="epos";   checkForces_norm( xname="epos",  fname="efpos", fnname="enfpos", inds=inds, xs=np.arange(-2.0,3.0,0.1) )
    label="ecoef";  checkForces_norm( xname="ecoef",  fname="efcoef", fnname="enfcoef", inds=inds, xs=np.arange( 0.0, 2.0, 0.01 ) )
    ##label="epos";  checkForces( xname="epos",  fname="enfpos",  inds=inds, xs=np.arange(-2.0,3.0,0.1) )
    #plt.plot( xs, dQdxa*-fnnScale, label="dQdxa_py" )
    plt.plot( xs, dQdca*-fnnScale, label="dQdca_py" )
    #effpy.plotNumDeriv( xs, E, dEdxa, F_=None, title="", bNewFig=False )
    effpy.plotNumDeriv( xs, E, dEdca, F_=None, title="", bNewFig=False )
    
    plt.grid()

def checkForces_He_2g( inds=(0,0), bNormalize=1 ):
    import CLCFGO_normalization_derivs_2 as effpy
    global label
    effmc.setPauliMode(0)
    #effmc.setPauliMode(2)
    #effmc.setSwitches_( normalize=bNormalize, normForce=bNormalize, kinetic=1, pauli=-1, coulomb=-1, AE=-1, AECoulomb=-1, AEPauli=-1 )   # Kinetic -  OK
    #effmc.setSwitches_( normalize=bNormalize, normForce=bNormalize, kinetic=-1, pauli=-1, coulomb=-1, AE=1, AECoulomb=1, AEPauli=-1 )   # AECoulomb -  OK
    #effmc.setSwitches_( normalize=bNormalize, normForce=bNormalize, kinetic=-1, pauli=-1, coulomb=-1, AE=1, AECoulomb=-1, AEPauli=1 )    # AEPauli   - OK
    #effmc.setSwitches_( normalize=bNormalize, normForce=bNormalize, kinetic=-1, pauli=-1, coulomb=1, AE=-1, AECoulomb=-1, AEPauli=-1 )    # ee Coulomb - OK 
    #effmc.setSwitches_( normalize=bNormalize, normForce=bNormalize, kinetic=-1, pauli=1, coulomb=-1, AE=-1, AECoulomb=-1, AEPauli=-1 )    # ee Pauli - probably OK (some small difference when energy is high for iPaulModel=2; iPauliModel=0 work well) 
    #effmc.setSwitches_( normalize=bNormalize, normForce=bNormalize, kinetic=1, pauli=1, coulomb=1, AE=1, AECoulomb=1, AEPauli=1 )    # Total - OK
    effmc.loadFromFile( "../../cpp/sketches_SDL/Molecular/data/He_2g_triplet_asym.fgo"  )
    #effmc.loadFromFile( "../../cpp/sketches_SDL/Molecular/data/H_2g-.fgo"  )
    effmc.printSetup()
    effmc.printAtomsAndElectrons()
    #label="epos";  checkForces_norm( xname="epos",  fname="efpos", fnname="enfpos", inds=inds, xs=np.arange(-2.0,3.0,0.1) )    
    #label="epos";  checkForces_norm( xname="epos",  fname="efpos", fnname="enfpos", inds=inds, xs=np.arange(-2.0,3.0,0.02) )  
    label="ecoef";  checkForces_norm( xname="ecoef",  fname="efcoef", fnname="enfcoef", inds=inds, xs=np.arange(-1.0,1.0,0.02) )  
    plt.grid()


def test_H2molecule():
    '''
    according to http://aip.scitation.org/doi/10.1063/1.3272671
    s_opt  = 0.94 A (1.77 bohr)
    lHHopt = 1.47 A     ()
    E_opt  = -2.91 [eV] (67kcal/mol)
    '''
    # ====== test H-atom
    init_effmc_H2mol_1g()
    #init_effmc( natom_=1, norb_=1, perOrb_=1, sz=0.5, dist=1.0, aQ=-1,aQsz=0.0, aP=0.0,aPsz=0.0 )
    #effmc.setSwitches_( normalize=-1, normForce=-1, kinetic=1, coulomb=1, pauli=1,    AA=1, AE=1, AECoulomb=1, AEPauli=-1 )
    #effmc.setSwitches_( normalize=-1, normForce=-1, kinetic=1, coulomb=1, pauli=1, AA=1, AE=1, AECoulomb=1, AEPauli=-1  )
    xs = np.arange(-3.0,3.0,0.05)
    plot_Terms( xs=xs, xname="epos", inds=(0,0) )

def opt_He_Triplet( xs=None, s0=0.5, s1=0.5 ):
    #effmc.init(1,2,2,1)
    #effmc.loadFromFile( "../../cpp/sketches_SDL/Molecular/data/He_2g_triplet_asym.fgo" )
    effmc.init(1,2,2,1)
    effmc.getBuffs( 1, 2, 2 )
    effmc.setSwitches_( normalize=1, normForce=1, kinetic=1, coulomb=1, pauli=1,    AA=-1, AE=1, AECoulomb=1, AEPauli=1 )
    #epos  = effmc.epos
    #esize = effmc.esize
    #ecoef = effmc.ecoef
    #Ebuf  = effmc.Ebuf
    #apos  = effmc.apos
    #aPars = effmc.aPars
    ff = effmc
    ff.apos [:,:]=0.0
    ff.aPars[:,0]=-2;ff.aPars[:,1]=0.2;ff.aPars[:,2]=0.2;ff.aPars[:,3]=0.0;
    nterm = len(ff.Ebuf)
    Es  = np.zeros((len(xs),nterm+1))
    Es_ = np.zeros((len(xs),nterm+1))
    # -- Orb1
    ff.epos [0,:,:] = 0.0
    ff.esize[0,:]   = s0
    ff.ecoef[0,:]   = 1.0
    # -- Orb2
    ff.esize[1,:] = s1
    ff.ecoef[1,0] =  1.0
    ff.ecoef[1,1] = -1.0
    for i in range(len(xs)):
        ff.epos[1,0,0] = +xs[i]
        ff.epos[1,1,0] = -xs[i]
        Es[i,0 ] = ff.eval()
        Es[i,1:] = ff.Ebuf[:]
    ff.ecoef[1,0] =  1.0
    ff.ecoef[1,1] =  1.0
    for i in range(len(xs)):
        ff.epos[1,0,0] = +xs[i]
        ff.epos[1,1,0] = -xs[i]
        Es_[i,0 ] = ff.eval()
        Es_[i,1:] = ff.Ebuf[:]
    if(plt):
        #print "test_ETerms PLOT"
        term_clr   = ['k', 'r', 'b',    'm',       '',      'g',   '',        ''  ]
        term_mask  = [1,     1,   1,      1,        0,        1,    0,        0   ]
        term_names = ["Etot","Ek","Eee","EeePaul","EeeExch","Eae","EaePaul","Eaa" ]
        plt.figure(figsize=(5,5))
        for i in range(nterm):
            if term_mask[i]==1: plt.plot( xs, Es [:,i],'-', c=term_clr[i], label=term_names[i] )
        for i in range(nterm):
            if term_mask[i]==1: plt.plot( xs, Es_[:,i],':', c=term_clr[i] )
        plt.grid();plt.legend();
        plt.title(label)  
    return xs, Es, Es_


if __name__ == "__main__":

    #np.random.seed( 451)
    #np.random.seed( 1)   # good
    np.random.seed( 3)
    #np.random.seed( 5)    # exact
    #np.random.seed( 8)
    import matplotlib.pyplot as plt_
    global plt,label
    #global dx,nx,x0_glob
    global ss_glob,xs_glob,cs_glob
    global rnd_pos, rnd_size, rnd_coef
    global iNorm
    x0_glob = 0.00001
    dx=0.05
    nx=100
    #nx=50
    #nx=10
    #nx=2
    plt=plt_
    bPrintInfo = True
    rnd_pos  = 0.0; rnd_size = 0.0; rnd_coef = 0.0
    #rnd_pos  = 0.2; rnd_size = 0.2; rnd_coef = 0.2
    #rnd_pos  = 0.2; rnd_size = 0.2; rnd_coef = 0.5
    #rnd_pos  = 1.0; rnd_size = 0.2; rnd_coef = 0.2

    ss_glob = np.arange( 0.1,3.0,0.01)
    xs_glob = np.arange(-2.0,3.0,0.01)
    cs_glob = np.arange(-1.0,2.0,0.01)
    
    #iNorm = +1

    #test_Hatom()     #; plt.show(); exit(0)
    compareForces_H_2g( inds=(0,0), bNormalize=True ); plt.show(); exit(0)
    #checkForces_He_2g( inds=(0,0), bNormalize=True ); plt.show(); exit(0)
    #compareForces_H_2g( inds=(0,0), bNormalize=False ); plt.show(); exit(0)
    #checkForces_H_2g( inds=(0,0) ); plt.show(); exit(0)

    #dx=0.02
    #opt_He_Triplet( np.arange(x0_glob,x0_glob+dx*100,dx), s0=0.5, s1=0.7 )
    
    #test_EvalFuncDerivs()  # plt.show; exit()
    #test_EvalFuncDerivs( s0=0.8 )

    effmc.setPauliMode(0)  # E = K*S^2
    #effmc.setPauliMode(2)  # E = Sij^2/(1-Sij^2) * ( Tii + Tjj - 2Tij/Sij )
    #effmc.setPauliMode(3)  # E=T
    #effmc.setPauliMode(4)  # E=S
    #effmc.setPauliMode(5)   # Ep = ( Sij/(1-Sij^2) )* Tij 
    #effmc.setPauliMode(6)   # Ep = Sij*Tij
    #test_ETerms( xname="epos", inds=(0,0), x0=0 );   # plt.show(); exit(0)
    #test_H2molecule()

    tests_results = []
    tests_funcs = []

    tests_funcs += [ test_ProjectWf, test_Poisson ]
    #tests_funcs += [ check_dS_epos           , check_dS_esize           , check_dS_ecoef             ]
    #tests_funcs += [ checkForces_Kinetic_epos, checkForces_Kinetic_esize, checkForces_Kinetic_ecoef  ]
    #tests_funcs += [ checkForces_Pauli_epos  , checkForces_Pauli_esize  , checkForces_Pauli_ecoef    ]
    tests_funcs += [ checkForces_Hartree_epos, checkForces_Hartree_esize, checkForces_Hartree_ecoef  ]
    #tests_funcs += [ checkForces_AQ_epos     , checkForces_AQ_esize     ,  checkForces_AQ_ecoef      ]
    #tests_funcs += [ checkForces_AP_epos     , checkForces_AP_esize     ,  checkForces_AP_ecoef      ]
    #tests_funcs += [ checkForces_AA_pos ]
    #tests_funcs += [ checkForces_Tot_epos, checkForces_Tot_esize ,  checkForces_Tot_ecoef  ]
    #tests_funcs += [ check_Coulomb_rhoP_, check_Coulomb_rhoS_, check_Coulomb_rhoQ_ ]
    #tests_funcs += [ check_Coulomb_rhoP_ ]
    #tests_funcs += [ check_Coulomb_rhoP_, check_Coulomb_rhoS_ ]
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

