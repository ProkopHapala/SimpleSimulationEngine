#!/usr/bin/python

import matplotlib.pyplot as plt 
import numpy as np
from pySimE.KosmoSuiteCpp import main as KS

datadir = "/home/prokop/git/SimpleSimulationEngine/cpp/apps/OrbitalWar/data/"

def makeLaunch( ts, maxThrust =3.5e+07 ):
    n = len(ts)
    thetaCPs  = 0.5*np.pi*np.exp( -0.15*ts**4 );
    thrustCPs = np.ones(n)*maxThrust
    dirCPs    = np.zeros((n,3))
    dirCPs[:,0] = np.cos(thetaCPs)
    dirCPs[:,1] = np.sin(thetaCPs)
    return thetaCPs,thrustCPs, dirCPs

aeroData = np.transpose( np.genfromtxt( datadir+'Cd.dat', skip_header=1 ) ).copy()
atmoData = np.transpose( np.genfromtxt( datadir+'standard_atmosphere.dat', skip_header=1 ) ).copy()

#print aeroData
#print atmoData

planet     = KS.Planet(6371e+3,3.9860044189e+14, KS.Vec3d(0.0,-6371e+3,0.0) )
rocket     = KS.Rocket(4500.0,2.970e+6,130e+3 + 140e+3 + 40e+3 + 13.5e+3,35e+6, 320.0)
atmosphere = KS.Atmosphere  ( atmoData[0,1]-atmoData[0,0], atmoData[3] )
aero       = KS.Aerodynamics( aeroData[0,1]-aeroData[0,0], aeroData[1]  )

ts = (np.arange(20)-1)*0.5
thetaCPs,thrustCPs, dirCPs = makeLaunch( ts )

launch     = KS.Launch ( thetaCPs, thrustCPs, dirCPs, 50.0, 7800, 200e+3 )
logTrig    = KS.LogTrig( 1.0, 1000.0 )

#print atmoData[3]
#KS.printStruc(atmosphere)
#KS.printStruc(launch)

KS.SpaceLaunchODE_init( planet, rocket, aero, atmosphere )
outbuff = KS.SpaceLaunchODE_run( 500, 1000000, launch, logTrig )
#print outbuff[0]
#print outbuff

# --- cut out not filled part of buffer
iend = np.argmax(outbuff[:,0])
outbuff = outbuff[:iend]

t  = outbuff[:,0]
m  = outbuff[:,4]
h  = outbuff[:,5]
v  = outbuff[:,6]
vx = outbuff[:,7]
vy = outbuff[:,8]
T  = outbuff[:,9]
Tx = outbuff[:,10]
Ty = outbuff[:,11]
G  = outbuff[:,12]


plt.figure(figsize=(8,8))
plt.plot(thetaCPs)

plt.figure(figsize=(8,8))

plt.subplot(3,1,1);
plt.plot(t,h*1e-3);
#plt.xlabel("t[s]")
plt.ylabel("h[m]");
plt.grid()

plt.subplot(3,1,2);
plt.plot(t,v *1e-3, label="v");
plt.plot(t,vx*1e-3, label="vx");
plt.plot(t,vy*1e-3, label="vy");
#plt.xlabel("t[s]")
plt.ylabel("vel[km/s]");
plt.legend()
plt.grid()

plt.subplot(3,1,3);
plt.plot(t,T /m, label="T");
plt.plot(t,Tx/m, label="Tx");
plt.plot(t,Ty/m, label="Ty");
plt.plot(t, G  , label="G");
#plt.plot(t,10*Fd/m, label="Fd*10");
plt.xlabel("t[s]")
plt.ylabel("acc[m/s2]");
plt.legend()
plt.grid()

plt.show()

 
