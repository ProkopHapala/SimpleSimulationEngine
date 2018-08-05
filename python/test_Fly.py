#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pyFlight as flight

flight.loadFromFile( "/home/prokop/git/SimpleSimulationEngine/cpp/apps/AeroCombat/data/AeroCraftMainWing1.ini" )

def getRotXZ( angle ):
    ca=np.cos(angle)
    sa=np.sin(angle)
    return np.array( [[ca,0.0,-sa],[0.0,1.0,0.0],[sa,0.0,ca]] )

rot = np.array( [[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]] )
flight.setTilt( 2, 0.01 )
flight.setTilt( 0,  0.1 )
flight.setTilt( 1, -0.1 )
#flight.setTilt( 3, 0.1 )

n = 10
poss=np.zeros((n,3  ))
vels=np.zeros((n,3  ))
rots=np.zeros((n,3,3))

'''
plt.figure(figsize=(15,5));
plt.subplot(1,3,1); plt.plot( poss[:,0], poss[:,2] ); plt.title("X-Z")
plt.subplot(1,3,2); plt.plot( poss[:,0], poss[:,1] ); plt.title("X-Y")
plt.subplot(1,3,3); plt.plot( poss[:,2], poss[:,1] ); plt.title("Z-Y")
plt.show()
'''

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

for angle in np.linspace(0.0,1.0,16,False)*np.pi*4.0:
    rot = getRotXZ( angle )
    print "angle", angle
    flight.setPose( np.array([0.0,200.0,0.0]), rot[2]*100.0, rot )
    flight.fly( poss, vels, rots, nsub=10, dt=0.01 )
    ps = poss.copy()
    ax.plot(ps[:,0], ps[:,2], ps[:,1] )

plt.axis('equal')

ax.set_xlabel('X'); ax.set_ylabel('Z'); ax.set_zlabel('Y');

plt.show()
