
import numpy as np
import matplotlib.pyplot as plt


data = np.transpose( np.genfromtxt( "O:\\vbs3_beta\\flyStraight.log" ) )
    
# 0  1      2          3        4     5     6       7        8         9         10      11          12           13             14         15            16              17
# i, t, _c_elevator, cpitch, cvup, dpitch, fy, _c_aileron, croll, _c_propel, _c_break, speed,  _velocity[0], _velocity[1], _velocity[2], _position[0], _position[1], _position[2] );

#i     = data[0]
t      = data[1]
elev   = data[2]
y      = data[16]
vy     = data[13]
cpitch = data[3]
vr     = data[11]
cvup   = data[4]
dpitch = data[5]
fy     = data[6]
speed  = data[11]


plt.figure()
#plt.plot(t,y-y[0], label="y")
#plt.plot(t,vy,     label="vy")
plt.plot(t,cpitch,        label="cpitch")
plt.plot(t,dpitch*10,     label="dpitch")
plt.plot(t,cvup,          label="cvup")
plt.plot(t,fy/np.abs(fy), label="fy")
#plt.plot(i,vr,     label="vr")
plt.plot(t,elev,   label="elev")

#print "elevator", elev
#print "cvup", cvup

plt.legend()
plt.axhline(0,ls="--",c="k")

plt.figure()
plt.plot(t,speed, label="speed")

plt.show()