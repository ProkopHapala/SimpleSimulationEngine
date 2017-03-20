#!/usr/bin/python

import numpy as np
import matplotlib.pylab as plt

'''
Mass driver is first limited by maximum acceleration projectile/line can sustain, than by maximum power source can provide
Therefore acceleration can be decomposed to two phases:

1) conctant acceleration:
v(t) =       a * t
s(t) = 0.5 * a * t**2 

2) constant power acceleration:
P = F*v = m*a*v

dv/dt   = P/(m*v)
v*dv    = (P/m) *dt 
v**2 /2 = (P/m) * t + C_
v**2    = 2P/m  * t + C 

for t = 0 :  v0**2 = C

v(t) = sqrt( 2*P/m * t + v0**2 )    


ds/dt = v(t)
ds    = sqrt( 2*P/m * t + v0**2 ) dt

T = 2*P/m * t + v0**2;   dT = 2*P/m * dt
 
 
ds    = 1/(2*P/m) sqrt( T ) dT
s     = (2/3)/(2*P/m) * (  2*P/m * t + v0**2   )**(3/2)   + C
s(t)  = (m/(3*P)) * (  2*P/m * t + v0**2   )**(3/2)    + s0

3) condition where power max is reached
a*t = P/(m*a)
t   = P/(m*a**2)


'''


# ========== Functions

def tswitch( PmaxSpec, a0 ):
    return PmaxSpec / a0**2

def acclelForceLimited( t, a ):
    v =       a * t
    s = 0.5 * a * t**2
    return v, s

def accelPowerLimited( t, PmaxSpec, v0, s0 ):
    v = np.sqrt( 2*PmaxSpec * t + v0**2 )
    s = ((v**3) - (v0**3) )/( 3*PmaxSpec)  + s0 
    return v, s
        
def getAccelCombined( ts, PmaxSpec, a0 ):
    tsw = tswitch( PmaxSpec, a0 )
    print tsw, ts[-1]
    # --- acceleration limited 
    t1    = np.append( ts[ ts < tsw ] , [tsw] )
    v1,s1 = acclelForceLimited(  t1, a0 )
    #1 --- power limited
    if( len(t1) < (len(ts)+1) ):
        t2    = ts[ ts > tsw ]
        v2,s2 = accelPowerLimited( t2-tsw, PmaxSpec, v1[-1], s1[-1] )
        t1 = np.concatenate((t1,t2))
        v1 = np.concatenate((v1,v2))
        s1 = np.concatenate((s1,s2))
    print t1[-1]
    return t1,v1,s1
    
def getTerminalVelocity( PmaxSpec, a0, sMax ):
    tsw = tswitch( PmaxSpec, a0 )
    ssw = 0.5 * a0 * tsw**2
    if(ssw>sMax):  # only Force limited
        tsmax = np.sqrt(2*sMax/a) 
        return a0 * tsmax
    else:          # both
        vsw = a0 * tsw
        #t   = (((sMax - ssw )*(3*PmaxSpec) + vsw**3)**(2.0/3.0) - vsw**2)/(2*PmaxSpec)
        v    = ( (sMax - ssw )*(3*PmaxSpec) + vsw**3 )**(1.0/3.0)
        return v

def main():
    Pspec = 1e+5;     # [W/kg]
    a     = 9.81 * 10 # [m/s^2] 
    t,v,s = getAccelCombined( np.linspace(0, 60.0, 100 ), Pspec, a ) 
    dt =  (t[1]-t[0]);
    plt.subplot(2,1,1); plt.title( "a=%g[G] P=%g[kW/kg]" %(a/9.81,Pspec*0.001) )
    plt.subplot(2,1,1); plt.grid(); plt.xlabel('time [s]'); plt.ylabel('Speed  [km/s]'); plt.plot(t,v*0.001,'.-'); #plt.plot(t[:-1]+dt*0.5,(s[1:]-s[:-1])/dt)   
    plt.subplot(2,1,2); plt.grid(); plt.xlabel('time [s]'); plt.ylabel('distance [km]'); plt.plot(t,s*0.001,'.-');  


if __name__ == "__main__":
    main()
    plt.show()

