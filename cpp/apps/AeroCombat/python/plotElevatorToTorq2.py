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
plt.plot(elev,AoA*(180.0/np.pi), label="AoA")
plt.plot(elev,L,   label="Lift")
plt.plot(elev,D,   label="Drag")
plt.plot(elev,LD, label="Lift/Drag")
plt.xlabel("elevator")
plt.legend()
plt.grid()
plt.axhline(0,ls="--",c="k")


# ======= Optimal Turn

# Ref 
# https://aviation.stackexchange.com/questions/2871/how-to-calculate-angular-velocity-and-radius-of-a-turn

imax  =np.argmax(LD)
#LDmax = LD[imax]
Lmax  = L[imax]
Dmax  = D[imax] 


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

def quadraticRoots( a, b, c ):
    D    = b**2 - 4.0*a*c
    mask = D<0
    D[mask] = np.NaN 
    sqD = np.sqrt(D)
    #print "sqD = ", sqD
    return (-b-sqD)/(2.0*a), (-b+sqD)/(2.0*a)
    
def evalThrustLimitedTurn( D,L, Thrust, mass, rho=1.2 ): 
    K    = 0.5*rho
    v2   = Thrust/(K*D)
    R    = mass/(K*L)
    return np.sqrt(v2),R

def evalThrustLimitedFlatTurn( D,L, Thrust, mass, rho=1.2, G=9.066 ): 
    K        = 0.5*rho
    v2       = Thrust/(K*D)
    #print "v", np.sqrt(v2)
    cosTheta = G*mass/(K*L*v2)
    mask = cosTheta>1.0
    cosTheta[mask] = np.NaN
    #print "cosTheta", cosTheta
    sinTheta = np.sqrt(1-(cosTheta**2))
    R        = mass/(K*L*sinTheta)
    return np.sqrt(v2),R,np.arcsin(sinTheta)

def turnTime( vs, R ):
    return np.pi*2*R/vs

def climbBranch( FdTh, Fg, KL, v2 ):
    #mask = (ThFd<0)
    #v2[mask] = np.NaN
    Fl = KL*v2;        
    sa = FdTh/Fg
    ca = Fl/Fg
    v = np.sqrt(v2)
    vx = ca*v
    vy = sa*v
    return vx,vy,np.arcsin(sa)

def evalClimbRate( D,L, Thrust, mass, rho=1.2, G=9.066 ): 
    '''
    (Ft-Fd)^2 + L^2 = G^2
    K=S*rho/2
    (Ft-K*Cd*v2)**2 + ( K*Cl*v2)**2 = (m*g)**2
    Ft**2 + 2*K*Cd*Ft*v2  + (K*Cd*v2)**2 + ( K*Cl*v2)**2 = (m*g)**2
    ((K*(Cd+Cl))**2)* v2**2 - (2*K*Cd*Ft)*v2 + (Ft**2 - (m*g)**2) = 0 
    '''
    #print "L =", L
    #print "D =", D
    K = 0.5*rho
    Fg = G*mass; #print "Fg = ", Fg    
    a = (K*D)**2 + (K*L)**2  ; #print "a = ", a
    b = -2*K*D*Thrust ; #print "b = ", b
    c =  Thrust**2 - Fg**2 ; #print "c = ", c
    vv1,vv2 = quadraticRoots( a, b, c )
    #print "vv1 = ",vv1
    #print "vv2 = ",vv2
    #vx1,vy1,alfa1 = climbBranch( Thrust-K*D*vv1, Fg, K*L, vv1 )
    #vx2,vy2,alfa2 = climbBranch( Thrust-K*D*vv2, Fg, K*L, vv2 )
    
    # back subs:
    #err = (Thrust-K*D*vv2)**2 + (K*L*vv2)**2 - Fg**2
    ##err =  Thrust**2  -2*K*D*vv2 +  (K*D*vv2)**2 +  (K*L*vv2)**2 - Fg**2
    ##err =  (K*D*vv2)**2 + (K*L*vv2)**2  - 2*K*D*Thrust*vv2 + Thrust**2 - Fg**2
    ##err =  ((K*D)**2 + (K*L)**2)*(vv2**2)  - 2*K*D*Thrust*vv2 + Thrust**2 - Fg**2
    #err = a*(vv2**2) + b*vv2 + c
    #print "err", err
    
    vx,vy,alfa = climbBranch( Thrust-K*D*vv2, Fg, K*L, vv2 )
    '''
    #v2 = np.maximum(vv1,vv2)
    #v2 = -np.minimum(vv1,vv2)
    Fd1  = K*D*vv1; mask = (Fd<Thrust)
    Fd2  = K*D*vv2; mask = (Fd>Thrust)
    v = np.sqrt(v2)
    print "v = ", v
    Fd = K*D*v2; print "Fd = ",Fd; print " Thrust-Fd ", Thrust-Fd
    #mask = (Fd>Thrust)
    #Fd[mask] = np.NaN;
    #v2[mask] = np.NaN;
    Fl = K*L*v2; print "Fl = ",Fl;
    sa = (Thrust-Fd)/Fg
    ca = Fl/Fg
    vx = ca*v
    vy = sa*v
    return vx,vy,np.arcsin(sa)
    '''
    return vx,vy,alfa
    #return (vx1,vy1,alfa1), (vx2,vy2,alfa2)


#plt.figure(figsize=(15,5))
plt.figure(figsize=(21,7))
plt.subplot( 1,3,1 )
plt.subplot( 1,3,2 )
plt.subplot( 1,3,2 )

def plotAero( mass, Thrust, L, D, clr=None, label="Aero1" ):
    
    vs,R0 = evalThrustLimitedTurn    ( D,L, Thrust, mass )
    vs,R,theta = evalThrustLimitedFlatTurn( D,L, Thrust, mass )
    T0 = turnTime(vs,R0)
    T  = turnTime(vs,R)

    #print "vs ",vs
    #print "R ",R    
    
    #print "T=",T
    vx,vy,alfa = evalClimbRate( D,L, Thrust, mass )
    #(vx1,vy1,alfa1), (vx2,vy2,alfa2) = evalClimbRate( D,L, Thrust, mass )
    #vx,vy,salfa = evalClimbRate( D[::10],L[::10], Thrust, mass, rho=1.2, G=9.066 )
    plt.subplot( 1,3,1 )
    plt.plot(vs,R0, c=clr, ls=":" )
    plt.plot(vs,R,  c=clr, ls="-", label=label )
    #plt.ylim(0,R[-1]*10)
    Rmin=np.nanmin(R); plt.ylim(0,Rmin*10); #print "Rmin ", Rmin
    plt.subplot( 1,3,2 )
    plt.plot(vs,T0, c=clr, ls=":" )
    plt.plot(vs,T,  c=clr, ls="-", label=label )
    #plt.ylim(0,T[-1]*5)
    Tmin=np.nanmin(T); plt.ylim(0,Tmin*5); #print " Tmin ", Tmin
    #plt.subplot( 1,2,3 )
    plt.subplot( 1,3,3 )
    plt.plot(vx,vy,  c=clr, ls="-", label=label )
    #plt.plot( np.sqrt(vx**2+vy**2), vy,  c=clr, ls="-", label=label )
    vymax = np.nanmax(vy); plt.ylim(-vymax*1.1,vymax*1.1); #print "vymax ", vymax
    #plt.ylim(-60,60) # [m/s]
    #plt.plot(vx2,vy2,  c=clr, ls="-", label=label )
    #plt.plot(vx1,vy1,  c=clr, ls=":" )    
    #plt.ylim(0,T[-1]*5)
    
    '''
    plt.figure()
    v = 350 # [m/s]
    v2= v**2
    K = 1.2*0.5
    Fg  = mass * 9.066
    Fd  = K*D*v2
    Fl  = K*L*v2
    Fx  = Thrust - Fd 
    Fy  = Fl - Fg
    plt.plot( L, Fx, label="Fx" )
    plt.plot( L, Fy, label="Fy" )
    plt.grid()
    plt.legend()
    '''
    
    
mass = 15e+3
G    = 9.066

plotAero( mass,   mass*G*0.36, L,   D,   clr='k', label="Aero1"      )
plotAero( mass*2, mass*G*0.36, L,   D,   clr='g', label="Aero1.Mass*=2"   )
plotAero( mass,   mass*G*0.36, L*2, D*2, clr='b', label="Aero1.Area*=2"   )
plotAero( mass,   mass*G*0.36*2,  L,   D,   clr='r', label="Aero1.Thrust*=2" )

plt.subplot( 1,3,1 )
plt.xlabel("Speed [m/s]")
plt.ylabel("Turn radius [m]")
plt.grid()
plt.legend()

plt.subplot( 1,3,2 )
plt.xlabel("Speed [m/s]")
plt.ylabel("Turn time [s]")
plt.grid()
plt.legend()

plt.subplot( 1,3,3 )
plt.xlabel("Speed [m/s]")
plt.ylabel("Rate of Climb [m/s]")
plt.grid()
plt.legend()

plt.show()