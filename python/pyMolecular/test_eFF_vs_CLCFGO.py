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

print( "\n ============ CLCFGO \n" );
effmc.loadFromFile( "../../cpp/sketches_SDL/Molecular/data/e2_1g_2o.fgo" )
#effmc.loadFromFile( "../../cpp/sketches_SDL/Molecular/data/H2O_1g_8o.fgo" )
effmc.printSetup()
effmc.printAtomsAndElectrons()
effmc.eval()

print( "\n ============ eFF \n" );
eFF.load_xyz("../../cpp/sketches_SDL/Molecular/data/e2_eFF.xyz")
#eff.load_xyz("../../cpp/sketches_SDL/Molecular/data/H2O_eFF.xyz")
eFF.info()
eFF.eval()

print( " ===== ALL DONE !!! \n" )