#!/usr/bin/python

import numpy as np
import pyVis3D as vis
from time import sleep

vis.lib.initWindow()
vis.lib.asyncLoop(1000000000)  # asynchronous window - does not wait for exit

'''
# --- spheres (atoms)
poss_sph = (np.random.random((100,3)) - 0.5 )*10.0
colors   = np.array( [[0.9,0.1,0.1] for xyz in poss_sph] )
radius   = np.ones(len(poss_sph))*0.5
print poss_sph.shape, colors.shape, radius.shape
glint_spheres = vis.spheres( poss_sph, colors, radius )
'''

'''
# --- polyline (curve)
ts   = np.linspace(0,2*np.pi,100)
poss = np.transpose( np.stack([ np.sin(ts), np.cos(ts), np.sin(ts*3) ]), (1,0) ).copy()
#print poss 
glint_curve = vis.polyline( poss )
'''

'''
# --- Tetrahedron ( tetrahedron wireframe )
verts =np.array( [ [-1.0,0.0,0.0], [+1.0,0.0,0.0], [0.0,-1.0,0.0], [0.0,+1.0,0.0], [0.0,0.0,-1.0], [0.0,0.0,+1.0] ])
edges =np.array( [ [0,2], [0,3], [0,4], [0,5],   [1,2], [1,3], [1,4], [1,5],   [2,4], [2,5], [3,4], [3,5]  ], dtype=np.int32)
faces =np.array( [ [0,4,2], [0,2,5], [0,3,4], [0,5,3],   [1,2,4], [1,5,2], [1,4,3], [1,3,5] ],                dtype=np.int32)
glint_lines = vis.lines( edges, verts, icolor = int("0xFF10FF10", 0) )
glint_trinagles = vis.triangles( faces, verts )
'''

# --- vector field
def linspace3D( rmin=(-1.0,-1.0,-1.0), rmax=(1.0,1.0,1.0), n=(100,100,100) ):
    XYZs = np.mgrid[0:n[0],0:n[1],0:n[2]]
    X = XYZs[0]*(rmax[0]-rmin[0])/float(n[0]) + rmin[0]
    Y = XYZs[1]*(rmax[1]-rmin[1])/float(n[1]) + rmin[1]
    Z = XYZs[2]*(rmax[2]-rmin[2])/float(n[2]) + rmin[2]
    return X,Y,Z
    
X,Y,Z = linspace3D( n=(30,30,30) )

#poss = np.stack([X,Y*3,Z*5])
#print poss
#print poss.shape
#print poss

R  = np.sqrt(X**2+Y**2+Z**2) 
dX = X/(1+R*R)
dY = Y/(1+R*R)
dZ = Z/(1+R*R)

poss = np.transpose( np.stack([X,Y,Z]),    (1,2,3,0) ).reshape(-1,3).copy() 
vecs = np.transpose( np.stack([dX,dY,dZ]), (1,2,3,0) ).reshape(-1,3).copy()

vis.vectors( vecs, poss, icolor=int("0xFF10FF10", 0) )

#vis.lines( edges, verts, icolor = int("0xFF10FF10", 0) )

#vis.lib.loop(1000000)  # synchronous window - wait for exit

# --- asynchronous real time modification


'''
for i in range(20):
    sleep(0.5)
    print( glint_curve, glint_spheres );
    # --- change curve shape
    vis.lib.erase( glint_curve )
    poss        = np.transpose( np.stack([ np.sin(ts), np.cos(ts), np.sin(ts*i) ]), (1,0) ).copy()
    glint_curve = vis.polyline( poss )
    # --- change sphere colors
    vis.lib.erase( glint_spheres );                                    print "DEBUG 1 "
    colors = np.array( [[0.05*i,0.5,1-0.05*i] for xyz in poss_sph] );  print "DEBUG 2 "
    glint_spheres = vis.spheres( poss_sph, colors, radius );           print "DEBUG 3 "
'''

for i in range(20):
    sleep(10)

'''
vis.lib.initWindow()
vis.lib.asyncLoop(10000)
for i in range(10):
    vis.lib.setGobVar( i*100.0001 );
    sleep(0.05)
'''




