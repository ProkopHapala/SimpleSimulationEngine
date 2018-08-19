#!/usr/bin/python

import numpy as np
#import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D
import pyFlight as flight
from pynput import keyboard
import time

# --- set workdir with configuration of aircraft, textures etc.
work_dir = "/home/prokop/git/SimpleSimulationEngine/cpp/Build/apps/AeroCombat"
flight.loadFromFile( work_dir+"/data/AeroCraftMainWing1.ini" )

# --- initialize visualization (optional, just for debugging)
fview = flight.FlightView( work_dir )

# --- set initial state of aircraft   setPose( position, velocity, rotMat==[left,up,forward] )
flight.setPose( np.array([0.0,200.0,0.0]), np.array([0.0,0.0,100.0]) , np.array([[1.0,0.0,0.0],[0.0,1.0,0.0],[0.0,0.0,1.0]]) )
flight.setTilt( 2, 0.02 ) # set Angle-of-incidence for wing #2 ( = elevator)
#flight.setTilt( 0,  0.1 )
#flight.setTilt( 1, -0.1 )
#flight.setTilt( 3, 0.1 )

# --- initialize targets Spheres (x,y,z,Radius)
targets = np.random.random((30,4));
targets[:, 0] += -0.5; targets[:, 0] *= 10000;  # x
targets[:, 2] += -0.5; targets[:, 2] *= 10000;  # z
targets[:, 1] *= 200;   # y=height
targets[:, 3] += 50.0   # radius
#print targets
flight.setTargets( targets )

# ---  initialize C-interface buffers
iPitch=0; iYaw=1; iRoll=2; iThrottle=3; iTrigger=4
controlBuff  = np.zeros(5)                   # python->C [pich,yaw,roll,throttle,trigger] differeces (rate-of-change)
stateBuff    = np.zeros(4*3 + 4)             # C->python [ pos.xyz, vel.xyz, forward.xyz, up.xyz       pich,yaw,roll,throttle  ]
targetHits   = np.zeros( (len(targets),2) )  # C->python if hit [ageMin,ageMax] of projectiles which hit it, if no-hit [+1e+300,-1e+300]

# ========= BEGIN Keyboard interface 
keyDict={
"w":False,
"s":False,
"a":False,
"d":False,
"q":False,
"e":False,
"*":False,
"/":False,
" ":False,
}

GO_ON = True

def on_press(key):
    try:
        if key == keyboard.Key.space:
            #print "shoot ON"
            keyDict[" "] = True
        elif key.char == 'c': 
            print " !!!! exit !!!!"
            global GO_ON 
            GO_ON = False
            return False
        elif key.char in keyDict:
            #print key.char
            keyDict[key.char] = True
    except:
        pass

def on_release(key):
    try:
        if key == keyboard.Key.space:
            #print "shoot OFF"
            keyDict[" "] = False
        elif key.char in keyDict:
            #print "release ", key.char
            keyDict[key.char] = False
    except:
        pass

listener = keyboard.Listener( on_press=on_press, on_release=on_release)
listener.start()
# ========= END Keyboard interface 


iframe = 0

tragetMoveSpeed = 10.0

while GO_ON: # repeat until pressed key "c"

    if tragetMoveSpeed > 0.0:
        #targets[:,0] += ( np.random.random(len(targets)) - 0.5 )*tragetDiffuseSpeed;
        #targets[:,2] += ( np.random.random(len(targets)) - 0.5 )*tragetDiffuseSpeed;
        phi = 0.01*iframe + np.array( range(len(targets)) )
        targets[:,0] += np.sin(phi)*tragetMoveSpeed
        targets[:,2] += np.cos(phi)*tragetMoveSpeed

    #print keyDict
    #print [k for k, v in keyDict.items() if v]
    print "# frame ", iframe
    # set controls by keyboard
    controlBuff[:]          = 0.0   # 0.0 means not change - preserve state of controls from previous step
    if   keyDict["w"]: 
        controlBuff[iPitch] = +1.0    # move elevator up,    value +1.0 is maximal rate of change
    elif keyDict["s"]: 
        controlBuff[iPitch] = -1.0    # move elevator down
    elif keyDict["a"]: 
        controlBuff[iRoll]  = +1.0    # move ailerons left
    elif keyDict["d"]: 
        controlBuff[iRoll]  = -1.0    # move ailerons right
    elif keyDict["q"]: 
        controlBuff[iYaw]   = -1.0    # move rudder left
    elif keyDict["e"]: 
        controlBuff[iYaw]   = +1.0    # move rudder right
    elif keyDict["*"]: 
        controlBuff[iThrottle] = +1.0  # increase engine power
    elif keyDict["/"]: 
        controlBuff[iThrottle] = -1.0  # decrease engine power
    #elif keyDict[keyboard.Key.space]: 
    elif keyDict[" "]:
        #print "python shoot"
        controlBuff[iTrigger ] =  1.0  # shoot when >0.5  
    if not ( keyDict["a"] or keyDict["d"] ):               # if roll keys not active, retract ailerons to neutral position
        controlBuff[iRoll ] = stateBuff[12+iRoll]*-10.0    # dRoll/dt =  roll * -10
    #if(niter==0):
    #    controlBuff[iTrigger] = 1
    #flight.fly( poss, vels, rots, nsub=10, dt=0.001 )
    
    # !!!!! HERE WE CALL THE SIMULATION !!!!!
    flight.flyAndShootTargets( controlBuff, stateBuff, targetHits, nsub=10, dt=0.003 )
    #print "control state ", stateBuff[12:]
    
    if( stateBuff[1]<0.0 ):
        print " pos.y < 0.0  =>  YOU CRASHED TO GROUND !!!"
        exit()
    
    # write out which targets where hit by how old projectiles ( age of projectiles is useful for feedback propagation to history ) 
    for i in xrange(len(targetHits)):
        if targetHits[i,0]<1e+200:
            print "hit target ", i," by projectile of age interval" , targetHits[i,:]
    
    # visualize - ( optional, just for debugging )
    fview.draw()
    time.sleep(0.05)  # delay for visualization - may be removed for actual training
    iframe+=1
    #if iframe > 500:
    #    exit()

