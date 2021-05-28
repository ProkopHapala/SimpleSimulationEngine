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

import matplotlib.pyplot as plt

bDebug = 0

# ========= Functions

def init_2x1_electrons( sz = 0.5, dist=1.0 ):
    # ---- setup
    natom=0; norb=2; perOrb=1
    effmc.init(natom,norb,perOrb,1)  #  natom, nOrb, perOrb, natypes
    ecoef = effmc.getBuff( "ecoef",(norb,perOrb)   )
    esize = effmc.getBuff( "esize",(norb,perOrb)   )
    epos  = effmc.getBuff( "epos" ,(norb,perOrb,3) )
    ecoef[:,:  ]=1
    esize[:,:  ]=sz
    epos [1,0,0]=dist

def init_2x2_electrons( sz = 0.5, dist=1.0 ):
    # ---- setup
    natom=0; norb=2; perOrb=2
    effmc.init(natom,norb,perOrb,1)  #  natom, nOrb, perOrb, natypes
    ecoef = effmc.getBuff( "ecoef",(norb,perOrb)   )
    esize = effmc.getBuff( "esize",(norb,perOrb)   )
    epos  = effmc.getBuff( "epos" ,(norb,perOrb,3) )
    ecoef[:,:  ]=1
    esize[:,:  ]=sz
    epos [:,1,0]=dist

def test_ProjectWf( plt = None, Etoll=1e-5 ):
    print "\n ============ test_ProjectWf ( rho = |wf|^2 )"
    init_2x2_electrons( sz = 0.5, dist=1.0 )
    #print "DEBUG 1 "
    # ---- test
    nps = 400+1
    ps  = np.zeros((nps,3))
    ps[:,0] = np.linspace( -5.0, 5.0, nps )
    effmc.eval() # we have to run it to project wavefuction to aux density
    wf  = effmc.orbAtPoints(ps)
    rho = effmc.rhoAtPoints(ps)
    #print "DEBUG 2 "
    wf2 = wf**2
    err = rho - wf2
    Err = np.sqrt( (err**2).sum()/len(err) )
    print " Error ", Err
    #Vh  = effmc.hartreeAtPoints(ps)
    #print "DEBUG 3 "
    if plt is not None:
        plt.figure(figsize=(5,5))
        #plt.subplot(2,1,1);
        plt.plot( ps[:,0],    wf, label='wf'  )
        plt.plot( ps[:,0],   rho, label='rho' )
        plt.plot( ps[:,0],   wf2,':', label='wf^2')
        plt.legend(); plt.grid()
        plt.title( " test_ProjectWf( rho = |wf|^2 )" )
    #print "DEBUG 4 "
    return Err

def test_Poisson( plt = None, Etoll=1e-5 ):
    print " ===== test_Poisson ( rho = dd_xyz V )"
    init_2x2_electrons( sz = 0.5, dist=1.0 )
    effmc.eval() # we have to run it to project wavefuction to aux density
    #print "DEBUG 1 "
    dx=0.03; R=3.0
    err2, rho, rho_ =  effmc.test_Poisson( dx=dx, Rmax=R )
    #err2, rho, rho_ =  effmc.test_Poisson( dx=dx, Rmax=R, useWf=False )
    Err = np.sqrt( err2/len(rho) )
    #print "DEBUG 2 "
    if(plt):
        xs=np.arange(0,R*2,dx)
        plt.figure(figsize=(5,5))
        plt.plot( xs, rho ,      label=('rho ' ) ); 
        plt.plot( xs, rho_, ":", label=('rho_' ) ); 
        plt.title( "test_Poisson( rho = dd_xyz V )" )
        plt.legend(); plt.grid()
    #print "DEBUG 3 "
    print " Error ", Err
    return Err

def test_OrbInteraction( plt = None, Etoll=1e-5, iMODE=1 ):
    labels=[ "NONE", "Overlap Sij", "Kinetic Tij", "Coulomb Kij", ]
    label=labels[iMODE]
    print " ===== test_OrbInteraction "+label
    #init_2x1_electrons( sz = 0.1, dist=0.0 )
    #init_2x1_electrons ( sz = 0.25, dist=0.0 )
    #init_2x1_electrons( sz = 0.5, dist=0.0 )
    init_2x1_electrons( sz = 0.75, dist=0.0 )
    #init_2x1_electrons( sz = 1.0, dist=0.0 )
    #init_2x1_electrons( sz = 1.5, dist=0.0 )
    #print "DEBUG 1 "

    effmc.eval() # we have to run it to project wavefuction to aux density
    #print "DEBUG 2 "
    dx=0.2; 
    #nint=30;
    nint=50
    err2, Ek, Ek_, f1, f2 =  effmc.test_OrbInteraction( iMODE=iMODE, io=0,jo=1, nint=nint, dx=dx, Rmax=5.0, bPrint=0, bSave=0  )
    #err2, rho, rho_ =  effmc.test_Poisson( dx=dx, Rmax=R, useWf=False )
    Err = np.sqrt( err2/len(Ek) )
    print "Error ", Err
    #print "DEBUG 3 "
    if(plt):
        #print "DEBUG 3.0 "
        xs=np.arange(0,dx*nint,dx)
        #print "DEBUG 3.0.1 "
        plt.figure(figsize=(10,5))
        plt.subplot(1,2,1)
        #print "DEBUG 3.1 "
        plt.plot( f1,      label=('f1' ) ); 
        plt.plot( f2, ":", label=('f2' ) );
        plt.title( "test_OrbInteraction funcs "+ label )
        plt.legend(); plt.grid()
        plt.subplot(1,2,2)
        #print "DEBUG 3.2 "
        #plt.figure(figsize=(5,5))
        plt.plot( xs, Ek ,    label=('Ek_ana' ) ); 
        plt.plot( xs, Ek_, ":", label=('Ek_num' ) ); 
        plt.title( "test_OrbInteraction "+label )
        #plt.ylim( 0, 20.0 )
        plt.legend(); plt.grid()
    #print "DEBUG 4 "
    return Err
def test_Overlap_Sij(plt=None):
    return test_OrbInteraction(iMODE=1,plt=plt)
def test_Kinetic_Tij(plt=None):
    return test_OrbInteraction(iMODE=2,plt=plt)
def test_Coulomb_Kij(plt=None):
    return test_OrbInteraction(iMODE=3,plt=plt)


#def test_Kinetic( plt = None, Etoll=1e-5 ):
#

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    tests_funcs = [ test_ProjectWf, test_Poisson, test_Overlap_Sij, test_Kinetic_Tij, test_Coulomb_Kij ]
    #tests_funcs = [ test_ProjectWf ,  test_Poisson  ]
    #tests_funcs = [ test_Overlap_Sij, test_Kinetic_Tij, test_Coulomb_Kij  ]
    #tests_funcs = [ test_Overlap_Sij  ]
    #tests_funcs = [ test_Coulomb_Kij  ]
    tests_results = []
    for test_func in tests_funcs:
        tests_results.append( test_func(plt) )
    print ""
    print " ##### Test Result Summary ##### "
    for i,test_func in enumerate(tests_funcs):
        print test_func.__name__," Error: ", tests_results[i]
    plt.show()
    print( " ===== ALL DONE !!! \n" )

