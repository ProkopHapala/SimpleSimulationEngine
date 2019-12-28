#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

def linearPowerLimitedRocket_dist( t, F, ve=10000, m0=1.0 ):
    A =  m0*ve 
    B = np.log( F*t - A ) 
    return F*ve*( t/F + A*B/F**2 ) - t*ve*B

def getAccFT( F=None, ve=10000, P=1e+9, m0=1.0, fmass=0.1 ):
    if F is None:
        F    = P/ve   ;print "F[N]"   , F
    tend = (1-fmass)*ve*m0/F           ;print "tend[s]", tend
    return tend, F


def linearPowerLimitedRocket( t, ve=10000, P=1e+9, m0=1.0  ):
    #F    = P/ve                         ;print "F[N]"   , F
    #tend = m0*ve*(fmass - 1)/(F*fmass) ;print "tend[s]", tend
    #tend =  (1-fmass)*ve*m0/F           ;print "tend[s]", tend
    tend, F = getAccFT( ve=ve, P=P, m0=m0 )
    a  = F/( m0 -F*t/ve )   #;primt "a[G]", a/9.81 
    v0 =  ve*np.log( m0*ve )
    v  = -ve*np.log( np.abs( m0*ve - F*t ) ) +  v0   
    #s = (  ve*t + t*v -  v*(m0*ve/F) )
    s = ve*t + v * ( t -  m0*ve/F )

    #s = linearPowerLimitedRocket_dist( t, F, ve=ve, m0=m0 )
    return s,v,a

P     = 10e+9
ve    = 10e+3
fmass = 0.1
m0    = 1.0


tend, F = getAccFT( ve=ve, P=P, m0=m0, fmass=fmass )

ts = np.linspace(0,tend,1000)

s,v,a = linearPowerLimitedRocket( ts, ve=ve, P=P, m0=m0 )

plt.figure( figsize=(5,9) )
plt.subplot(3,1,1); plt.plot( ts, a ); plt.ylabel('a'); plt.xlabel('t[s]'); plt.grid()
plt.axvline( tend, ls="--")
plt.subplot(3,1,2); plt.plot( ts, v ); plt.ylabel('v [m/s]'); plt.xlabel('t[s]') ; plt.grid()
plt.axvline( tend, ls="--")
plt.subplot(3,1,3); plt.plot( ts, s ); plt.ylabel('s [m]  '); plt.xlabel('t[s]') ; plt.grid()

plt.show()