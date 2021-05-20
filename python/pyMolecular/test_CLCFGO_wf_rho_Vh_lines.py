import os
import sys
import numpy as np

#sys.path.append("../")

sys.path.append('../')
from pyMeta import cpp_utils 
cpp_utils.clean_build    = False   # Recompile only if changed
#cpp_utils.recompile_glob = False  # don't recompile

import eFF
import CLCFGO as effmc
import eFF_terms as effpy

import matplotlib.pyplot as plt

print( "\n ============ CLCFGO \n" );
#effmc.loadFromFile( "../../cpp/sketches_SDL/Molecular/data/e2_1g_2o.fgo" )
#effmc.setPauli(2)
#effmc.printSetup()
#effmc.printAtomsAndElectrons()

#natom=0; norb=2; perOrb=1;
natom=0; norb=2; perOrb=2;

effmc.init(natom,norb,perOrb,1)  #  natom, nOrb, perOrb, natypes
ecoef = effmc.getBuff( "ecoef",(norb,perOrb)   )
esize = effmc.getBuff( "esize",(norb,perOrb)   )
epos  = effmc.getBuff( "epos" ,(norb,perOrb,3) )

ecoef[:,:  ]=1
esize[:,:  ]=0.5
#epos [:,:,:]=0.5
epos [:,1,0]=1.0
#epos [:,1,0]=0.7
#epos [:,1,0]=0.0
#epos [:,1,0]=1.4

nps = 400+1
ps  = np.zeros((nps,3))
ps[:,0] = np.linspace( -5.0, 5.0, nps )


effmc.eval() # we have to run it to project wavefuction to aux density
wf  = effmc.orbAtPoints(ps)
rho = effmc.rhoAtPoints(ps)
Vh  = effmc.hartreeAtPoints(ps)

plt.figure(figsize=(10,5))
plt.subplot(2,1,1);
plt.plot( ps[:,0],    wf, label='wf'  )
plt.plot( ps[:,0], wf**2, label='wf^2')
plt.plot( ps[:,0],   rho, label='rho' )
plt.legend(); plt.grid()
plt.subplot(2,1,2);
plt.plot( ps[:,0], Vh, label='Vh' )
plt.legend(); plt.grid()


'''
print( " ===== Test wf wf^2 rho Vhartree for different width !!! \n" )
plt.figure(figsize=(10,5))
#for sz in [0.1,0.2,0.4,0.8,1.6]:
for sz in [0.2,0.4,0.8]:
    esize[:,:  ]=sz
    effmc.eval() # we have to run it to project wavefuction to aux density
    wf  = effmc.orbAtPoints(ps)
    rho = effmc.rhoAtPoints(ps)
    Vh  = effmc.hartreeAtPoints(ps)
    plt.subplot(2,1,1); 
    plt.plot( ps[:,0], rho, label=('rho s=%1.1f' %sz) ); 
    #plt.plot( ps[:,0], wf**2, ":", label=('wf^2 s=%1.1f' %sz) ); 
    plt.legend(); plt.grid()
    plt.subplot(2,1,2); 
    plt.plot( ps[:,0], Vh,  label=('Vh s=%1.1f' %sz) ); 
    plt.legend(); plt.grid()
'''


print( " ===== Test Poisson !!! \n" )
plt.figure()
effmc.eval()
#esize[:,:  ]=1.0
#esize[:,:  ]=0.5
esize[:,:  ]=0.35
#esize[:,:  ]=0.25
dx=0.03; R=3.0
xs=np.arange(0,R*2,dx)
err2, rho, rho_ =  effmc.test_Poisson( dx=dx, Rmax=R )
#err2, rho, rho_ =  effmc.test_Poisson( dx=dx, Rmax=R, useWf=False )
plt.plot( xs, rho ,      label=('rho ' ) ); 
plt.plot( xs, rho_, ":", label=('rho_' ) ); 
plt.legend(); plt.grid()


plt.show()

print( " ===== ALL DONE !!! \n" )