import os
import sys
import numpy as np

#sys.path.append("../")

sys.path.append('../')
from pyMeta import cpp_utils 
cpp_utils.clean_build    = False   # Recompile only if changed

import eFF as eff
import eFF_terms as effpy

# ========== Globals

bDebug = 0
natom=0; nelec=0
bPrintInfo = False
label=""
plt  =None

rnd_pos  = 0
rnd_size = 0
rnd_coef = 0

# ========= Functions

def init_eff( natom_=0, nelec_=1, s=0.5,  aQ=1.0,aQs=0.1,aP=0.0,aPs=0.1 ):
    global natom,nelec
    natom=natom_; nelec=nelec_; 
    eff.init( natom, nelec )
    aPar  = eff.getBuff( "aPars",(natom,4) )
    apos  = eff.getBuff( "apos",(natom,3) )
    epos  = eff.getBuff( "epos",(nelec,3) )
    esize = eff.getBuff( "esize",(nelec)  )
    aPar [:,0]=aQ;aPar[:,1]=aQs;aPar[:,2]=aPs;aPar[:,3]=aP;
    apos [:,:] = 0
    epos [:,:] = 0
    esize[:]   = s
    '''
    epos [:,:,:]= 0              + (np.random.rand(norb,perOrb,3)-0.5)*rnd_pos
    esize[:,:  ]=sz              + (np.random.rand(norb,perOrb  )-0.5)*rnd_size
    ecoef[:,:  ]=1               + (np.random.rand(norb,perOrb  )-0.5)*rnd_coef
    rhoP [:,:,:]=0               + (np.random.rand(norb,nqOrb,3 )-0.5)*rnd_pos
    rhoS [:,:  ]=sz*np.sqrt(0.5) + (np.random.rand(norb,nqOrb   )-0.5)*rnd_size
    rhoQ [:,:  ]=1               + (np.random.rand(norb,nqOrb   )-0.5)*rnd_coef
    '''

if __name__ == "__main__":
    import matplotlib.pyplot as plt

    init_eff( natom_=1, nelec_=1, s=0.5 )
    ss = np.arange( 0.4,3.0, 0.05 )
    Ek,Eae = effpy.Hatom( ss );  E_ref=Ek+Eae 
    E,dE   = eff.evalFuncDerivs(1,ss)
    #E_lim  = np.sqrt(8/np.pi)* (1/ss)    ; plt.plot( ss, E_lim,         label="Elim"  )
    #dE_lim = np.sqrt(8/np.pi)* (1/ss**2) ;     plt.plot( ss, dE_lim,        label="dElim" )
    #print "ss    ", ss
    #print "E     ", E
    #print "E_ref ", E_ref
    plt.plot( ss, E_ref,':',lw=3, label="E_ref" )
    plt.plot( ss, Eae,  ':',     label="Eae_ref" )
    plt.plot( ss, Ek,   ':',     label="Ek_ref" )
    plt.plot( ss, E,    'k',lw=3,label="E"     )
    
    xs=ss
    plt.plot(xs      ,dE                               ,'-',label="F_ana")
    plt.plot(xs[1:-1],(E[2:]-E[:-2])/(-2*(xs[1]-xs[0])),':',label="F_num")

    plt.xlabel('size [A]')
    plt.legend()
    plt.grid()
    plt.show()