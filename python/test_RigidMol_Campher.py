#!/usr/bin/python

import sys
import os
import numpy as np
import pyMolecular.RigidMol as rmol

os.chdir("/u/25/prokoph1/unix/git/SimpleSimulationEngine/cpp/Build/apps/MolecularEditor2")

# ========= Molecules

rmol.initParams( "common_resources/AtomTypes.dat", "common_resources/BondTypes.dat" )
water = rmol.loadMolType( "inputs/Campher.xyz" );

rot = np.array([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]])

print rmol.insertMolecule( water, np.array([ 5.78, 6.7, 12.24 ]), rot, True );

rmol.save2xyz( "world_debug_00.xyz" )

# ========= RigidSurface

rmol.initRigidSubstrate ( "inputs/Cu111_6x6_2L.xyz", np.array([60,60,100],dtype=np.int32), np.array([0.0,0.0,0.0]), np.array([[15.31593,0.0,0.0],[0.0,13.26399,0.0],[0.0,0.0,20.0]]) )
#rmol.initRigidSubstrate( "inputs/NaCl_wo4.xyz", np.array([60,60,100],dtype=np.int32), np.array([0.0,0.0,0.0]), np.array([[12.0173,0.0,0.0],[0.0,12.0173,0.0],[0.0,0.0,20.0]]) )


if os.path.isfile("data/FFPauli.bin"):
    print "gridFF found on disk => loading "
    rmol.loadGridFF()
else:
    print "gridFF not found on disk => recalc "
    rmol.recalcGridFF( np.array([1,1,1],dtype=np.int32) )
    rmol.saveGridFF()

rmol.debugSaveGridFF( "FFtot_z_Na.xsf", np.array([1.3,0.0447214,0.0]) )

# ========= Relaxation

rmol.bakeMMFF()

rmol.prepareOpt()

fout = rmol.openf( "movie.xyz", -1, "w" )

for i in range(550):
    print ">>> i ", i
    F2 = rmol.relaxNsteps( 1, 0.0 ); 
    print "|F| ", np.sqrt(F2)
    rmol.write2xyz( fout )
    #rmol.save2xyz( "world_debug_%03i.xyz" %i )

