#!/usr/bin/python

import numpy as np
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
import pyFlight as flight
import time

work_dir = "/home/prokop/git/SimpleSimulationEngine/cpp/Build/apps/AeroCombat"
flight.loadFromFile( work_dir+"/data/AeroCraftMainWing1.ini" )

fview = flight.FlightView( work_dir )

flight.setPose( np.array([0.0,200.0,0.0]), np.array([0.0,0.0,100.0]) , np.array([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]) )
flight.setTilt( 2, 0.02 )
#flight.setTilt( 0,  0.1 )
#flight.setTilt( 1, -0.1 )
#flight.setTilt( 3, 0.1 )

targets = np.random.random((30,4));
targets[:,:3] += -0.5; 
targets[:, 0] *= 1000;
targets[:, 1] *= 200;
targets[:, 2] *= 1000; 
targets[:, 3] += 2.0
#print targets

flight.setTargets( targets )

n    = 10
poss=np.zeros((n,3  ))
vels=np.zeros((n,3  ))
rots=np.zeros((n,3,3))

for i in range(100000):
    flight.fly( poss, vels, rots, nsub=10, dt=0.001 )
    fview .draw()
    time.sleep( 0.05)
