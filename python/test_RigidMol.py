#!/usr/bin/python

import sys
import os
import numpy as np
import pyMolecular.RigidMol as rmol

os.chdir("/u/25/prokoph1/unix/git/SimpleSimulationEngine/cpp/Build/apps/MolecularEditor2")

rmol.initParams( "common_resources/AtomTypes.dat", "common_resources/BondTypes.dat" )
water = rmol.loadMolType( "inputs/water_T5_ax.xyz" );
Na    = rmol.loadMolType( "inputs/NaIon.xyz" );
Cl    = rmol.loadMolType( "inputs/ClIon.xyz" );

rot = np.array([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]])

print rmol.insertMolecule( water, np.array([0.0,0.0,4.0]), rot, True );
print rmol.insertMolecule( water, np.array([4.0,0.0,4.0]), rot, True );

print rmol.insertMolecule( Na, np.array([4.0,6.0,5.0]), rot, True );
print rmol.insertMolecule( Na, np.array([4.0,4.0,2.0]), rot, True );
print rmol.insertMolecule( Na, np.array([4.0,8.0,2.0]), rot, True );

print rmol.insertMolecule( Cl, np.array([2.0,6.0,2.0]), rot, True );
print rmol.insertMolecule( Cl, np.array([6.0,6.0,2.0]), rot, True );

rmol.bakeMMFF()

rmol.save2xyz( "world_debug_1.xyz" )



