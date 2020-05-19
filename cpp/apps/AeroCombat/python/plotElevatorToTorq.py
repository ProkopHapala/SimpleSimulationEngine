# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt


data = np.transpose( np.genfromtxt( "O:\\vbs3_beta\\elevatorToTorq.log" ) )
    
elev   = -data[1]
AoA    = data[2]
L      = data[3]
D      = data[4]
LD = L/D

plt.figure()

plt.plot(elev,AoA, label="AoA")
plt.plot(elev,L,   label="L")
plt.plot(elev,D,   label="D")
plt.plot(elev,LD, label="L/D")

plt.xlabel("elevator")

plt.legend()
plt.axhline(0,ls="--",c="k")

# ======= Optimal Turn


# Ref 
# https://aviation.stackexchange.com/questions/2871/how-to-calculate-angular-velocity-and-radius-of-a-turn

imax  =np.argmax(LD)
#LDmax = LD[imax]
Lmax  = L[imax]
Dmax  = D[imax] 
 
def evalThrustLimitedTurn( D,L, Thrust, mass, rho=1.2 ): 
    v2   = Thrust/(0.5*rho*D)
    R    = mass/(rho*L)
    return np.sqrt(v2),R

def evalThrustLimitedFlatTurn( D,L, Thrust, mass, rho=1.2, G=9.066 ): 
    v2       = Thrust/(0.5*rho*D)
    #print "v", np.sqrt(v2)
    cosTheta = G*mass/(L*0.5*rho*v2)
    mask = cosTheta>1.0
    cosTheta[mask] = np.NaN
    #print "cosTheta", cosTheta
    sinTheta = np.sqrt(1-(cosTheta**2))
    R        = mass/(rho*L*sinTheta)
    return np.sqrt(v2),R,np.arcsin(sinTheta)

def turnTime( vs, R ):
    return np.pi*2*R/vs
   
def evalThrustLimitedTurnV( v, Ds,Ls, Thrust, mass, rho=1.2 ):
    v2   = v**2    
    pref = 0.5*rho*v2
    #print "D ", Ds
    Dreq = Thrust/pref;              
    mask = Dreq<Ds[-2]
    Dreq = Dreq[mask];        #print "Dreq ", Dreq
    v = v[mask]
    #Dreq[Dreq>Ds[-2]]=np.NaN; print "Dreq ", Dreq
    i    = np.searchsorted(Ds,Dreq); print "i ", i
    D0   = Ds[i  ]
    D1   = Ds[i+1]
    L0   = Ls[i  ]
    L1   = Ls[i+1]
    #print "L0 ", L0
    #print "L1 ", L1
    f    = (Dreq-D0)/(D1-D0); mf = 1-f
    D    = D0*mf + D1*f 
    L    = L0*mf + L1*f 
    #Fl   = L*pref
    R    = mass/(rho*L)
    return v,R,L,D

vs   = np.linspace(100.0,500.0,20)
mass = 15.0e+3
G = 9.8066

vs,R0 = evalThrustLimitedTurn    ( D,L, mass*G*0.36, mass )
vs,R,theta = evalThrustLimitedFlatTurn( D,L, mass*G*0.36, mass )
T0 = turnTime(vs,R0)
T = turnTime(vs,R)

plt.figure()
plt.plot(vs,R0, label="TurnRadius")
plt.plot(vs,R,  label="TurnRadiusFlat")
plt.xlabel("Speed [m/s]")
plt.ylabel("Turn radius [m]")
plt.ylim(0,R[-1]*5)
plt.legend()

plt.figure()
plt.plot(vs,T0, label="TurnTime")
plt.plot(vs,T,  label="TurnTimeFlat")
plt.xlabel("Speed [m/s]")
plt.ylabel("Turn time [s]")
plt.ylim(0,T[-1]*2)
plt.legend()

plt.figure()
plt.plot(vs,theta,         label="Theta")
plt.plot(vs,np.sin(theta), label="cR")
plt.plot(vs,np.cos(theta), label="cUp")
plt.xlabel("Speed [m/s]")
plt.legend()


plt.show()