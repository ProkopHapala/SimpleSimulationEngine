#!/usr/bin/python

import numpy as np
import pyMolecular as mol
import pyMolecular.testing as moltest
import matplotlib.pyplot as plt

# ==================== Compare two point distributions (permutation inveriant)

'''

points_ref = np.array([
 [1.0,0.0,0.0],   [-1.0, 0.0, 0.0], 
 [0.0,1.0,0.0],   [ 0.0,-1.0, 0.0],
 [0.0,0.0,1.0],   [ 0.0, 0.0,-1.0]
], dtype=np.float64 )

mol.initComparator( points_ref )

points = points_ref.copy()
dist   = mol.compDistance( points );   print( "dist (identical)", dist )

drnd    = np.random.rand( points.shape[0], points.shape[1] )
points += drnd*0.01
print "========== drandom"
print points
dist   = mol.compDistance( points.copy() );   print( "dist (drandom)", dist )

np.random.shuffle(points)
print points
dist = mol.compDistance( points.copy() );   print( "dist (shuffled)", dist )
'''

# ==================== Compare two TypePoint distributions (like atoms in molecule with different atom types)  (permutation inveriant)

'''
atoms=np.genfromtxt( "/home/prokop/git/SimpleSimulationEngine/cpp/apps/MolecularEditor/inputs/PTCDA/PTCDA.bas", skip_header=1 )

#print "atoms=", atoms

points_ref = atoms[:,1:4].copy();
types_ref  = atoms[:,0  ].astype(np.int32).copy();   

print points_ref

mol.initComparatorT     ( points_ref, types_ref )


print "========= identical"
points = atoms[:,1:4].copy();
types  = atoms[:,0  ].astype(np.int32).copy();    print( "types = ", types)
dist   = mol.compDistanceT( points_ref, types_ref );  print " >>> dist = ", dist 

print "========= shuffled"
np.random.shuffle(atoms);   
points = atoms[:,1:4].copy();
types  = atoms[:,0  ].astype(np.int32).copy();    print( "types = ", types)
dist   = mol.compDistanceT( points, types );     print " >>> dist = ", dist 

print "========= drandom"
drnd    = np.random.rand( points.shape[0], points.shape[1] )
points += drnd*0.01
dist    = mol.compDistanceT( points, types );  print " >>> dist = ", dist 
'''

# ==================== Compute fast Hash by plane waves projection of atomic coordinets (permutation inveriant)

'''
atoms=np.genfromtxt( "/home/prokop/git/SimpleSimulationEngine/cpp/apps/MolecularEditor/inputs/PTCDA/PTCDA.bas", skip_header=1 )

ks = np.array([
  [1.0,0.0,0.0],
  [0.0,1.0,0.0],
  [0.0,0.0,1.0]
])

points_ref = atoms[:,1:4].copy();
coefs_ref  = mol.getPlaneWaveDescriptor( points_ref, ks );   print "coefs (ref)       ", coefs_ref

points_1 = points_ref.copy()
coefs = mol.getPlaneWaveDescriptor( points_1, ks );          print "coefs (identical) ", coefs

np.random.shuffle(points_1); points_1 = points_1.copy() 
coefs = mol.getPlaneWaveDescriptor( points_1, ks );          print "coefs (shufled)   ", coefs

points_3 = points_ref.copy() + np.random.rand( len(atoms), 3 ) * 0.25
coefs = mol.getPlaneWaveDescriptor( points_3, ks );          print "coefs (drand)     ", coefs
'''

# ==================== Testing of statistical poperties of plane-wave hash

nrep   = 10
natoms = 100
Ns = range( 1, natoms )
dx = 0.5
k  = 3.0 

'''
xs = np.linspace(-10.0,10.0,1000)
ys = moltest.saw_sine( xs+100 )
plt.plot( xs, ys )
'''

xs_ref   = np.random.rand( natoms ); #print "xs_ref = ", xs_ref

#xs   = moltest.mutateN( xs_ref.copy(), 3, 0.1 ); print "xs = ", xs
coef_ref = moltest.hash_saw( xs_ref, k )

result = np.zeros((len(Ns)*nrep,2))

ires = 0
for N in Ns:
    dx_ = dx/float(N)
    for i in range(nrep):
        xs   = moltest.mutateN( xs_ref.copy(), N, dx_ )
        coef = moltest.hash_saw( xs, k )
        result[ ires, 0 ] = N; result[ ires, 1 ] = coef; 
        ires+=1

plt.axhline(coef_ref)    
plt.plot( result[:,0], result[:,1], '.' )

plt.show()


