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

print( "\n ============ CLCFGO \n" );
effmc.loadFromFile( "../../cpp/sketches_SDL/Molecular/data/e2_1g_2o.fgo" )
#effmc.loadFromFile( "../../cpp/sketches_SDL/Molecular/data/H2O_1g_8o.fgo" )
effmc.setPauli(2)
effmc.printSetup()
effmc.printAtomsAndElectrons()

print "effmc dims ", effmc.getDimPointer()
effmc_esize = effmc.getBuff("esize",2)    
#effmc_esize[:] *= np.sqrt(2)   # ToDo : :If esize multiplied by sqrt(2)   Ek and Eee agrees (eFF-vs-CLCDGO) but EeePaul, if without factor sqrt(2) EeePaul agree but Ek and Eee does not !!!
print " effmc_esize ", effmc_esize

E = effmc.eval()
print " # Ek Eee EeePaul EeeExch Eae EaePaul Eaa"
print "E_terms ", effmc.getEnergyTerms()
print "Etot ", E


print( "\n ============ eFF \n" );
eFF.load_xyz("../../cpp/sketches_SDL/Molecular/data/e2_eFF.xyz")
#eff.load_xyz("../../cpp/sketches_SDL/Molecular/data/H2O_eFF.xyz")
eFF.setPauliModel(2)
eFF.info()

print "eFF dims ", eFF.getDimPointer()
eFF_esize = eFF.getBuff("esize",2)      
print " eFF_esize ", eFF_esize

E = eFF.eval()
print " # Ek Eee EeePaul EeeExch Eae EaePaul Eaa"
print "E_terms ", eFF.getEnergyTerms()
print "Etot ", E 

print "effpy.Kinetic(0.5)*2 ", effpy.Kinetic(0.5)*2 


print( " ===== ALL DONE !!! \n" )