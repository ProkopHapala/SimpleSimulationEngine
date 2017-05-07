#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt


DATA = np.transpose( np.genfromtxt( "SpaceFlightODE.log" ) )

t  = DATA[0]
h  = DATA[1]

v  = DATA[2]
vx = DATA[3]
vy = DATA[4]

T  = DATA[5]
Tx = DATA[6]
Ty = DATA[7]
m  = DATA[8]
G  = DATA[9] 
Fd = DATA[10]

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
plt.plot(t,10*Fd/m, label="Fd*10");
plt.xlabel("t[s]")
plt.ylabel("acc[m/s2]");
plt.legend()
plt.grid()


plt.show()
