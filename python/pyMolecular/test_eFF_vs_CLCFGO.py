import os
import sys
import numpy as np

#sys.path.append("../")

sys.path.append('../')
from pyMeta import cpp_utils 
cpp_utils.clean_build    = False  # Recompile only if changed
#cpp_utils.recompile_glob = False  # don't recompile

import eFF
import CLCFGO as effmc
import eFF_terms as effpy

import matplotlib.pyplot as plt

iPauli = 1

def sampleEnergy( func=None, n=30, dx=0.1, buff=None ):
    Es = np.zeros(n)
    xs = np.arange(0,n*dx-1e-8,dx) 
    ind = (0,)*len(buff.shape)
    for i in range(n):
        buff[ind] += xs[i]
        Es[i] = func()
    plt.plot( xs, Es  )
    return Es

#==================================================
#==================================================
print( "\n ============ CLCFGO \n" );

effmc.loadFromFile( "../../cpp/sketches_SDL/Molecular/data/e2_1g_2o_singlet.fgo" )
#effmc.loadFromFile( "../../cpp/sketches_SDL/Molecular/data/e2_1g_2o_triplet.fgo" )
#effmc.loadFromFile( "../../cpp/sketches_SDL/Molecular/data/H2O_1g_8o.fgo" )
effmc.setPauliMode(iPauli)
effmc.printSetup()
effmc.printAtomsAndElectrons()

print "effmc dims ", effmc.getDimPointer()
effmc_esize = effmc.getBuff("esize",2)    
effmc_pos = effmc.getBuff("epos" ,3)  
#effmc_esize[:] *= np.sqrt(2)   # ToDo : :If esize multiplied by sqrt(2)   Ek and Eee agrees (eFF-vs-CLCDGO) but EeePaul, if without factor sqrt(2) EeePaul agree but Ek and Eee does not !!!
print " effmc_esize ", effmc_esize


E = effmc.eval()
print " # Ek Eee EeePaul EeeExch Eae EaePaul Eaa"
print "E_terms ", effmc.getEnergyTerms()
print "Etot ", E

sampleEnergy( effmc.eval, buff=effmc_pos )


#==================================================
#==================================================
print( "\n ============ eFF \n" );

eFF.load_fgo( "../../cpp/sketches_SDL/Molecular/data/e2_1g_2o_singlet.fgo" )
#eFF.load_fgo( "../../cpp/sketches_SDL/Molecular/data/e2_1g_2o_triplet.fgo" )
#eFF.load_xyz("../../cpp/sketches_SDL/Molecular/data/e2_eFF_singlet.xyz")
#eFF.load_xyz("../../cpp/sketches_SDL/Molecular/data/e2_eFF_triplet.xyz")
#eff.load_xyz("../../cpp/sketches_SDL/Molecular/data/H2O_eFF.xyz")
eFF.setPauliModel(iPauli)
eFF.info()

print "eFF dims ", eFF.getDimPointer()
eFF_esize = eFF.getBuff("esize",1) 
eFF_epos = eFF.getBuff("epos" ,2)
print " eFF_esize ", eFF_esize

E = eFF.eval()
print " # Ek Eee EeePaul EeeExch Eae EaePaul Eaa"
print "E_terms ", eFF.getEnergyTerms()
print "Etot ", E 


sampleEnergy( eFF.eval,buff=eFF_epos )





print "effpy.Kinetic(0.5)*2 ", effpy.Kinetic(0.5)*2 

plt.grid()
plt.show()

print( " ===== ALL DONE !!! \n" )