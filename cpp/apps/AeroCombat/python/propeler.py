#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

def quadratic_roots(a,b,c):
    D     = b*b - 4*a*c;
    sqrtD = np.sqrt( D );
    ia    = -0.5/a;
    x1    = ( b - sqrtD )*ia;
    x2    = ( b + sqrtD )*ia;
    # swap
    mask     = ia<0
    xtmp     = x1[mask]
    x1[mask] = x2[mask]
    x2[mask] = xtmp 
    return x1,x2

def propelerThrust( v0, power, area=1.0, efficiency=1.0, CD=0.0, dens=1.22 ):
    vstatic=  (power/area*dens)**0.33333
    print "vstatic: ", vstatic
    dm = dens*area*(v0+vstatic);
    a  = 0.5*dm;
    b  = dm*v0;
    c  = -power;
    Dv1, Dv2  = quadratic_roots( a, b, c );
    return dm * Dv2 * efficiency #- dm*v0*CD

power = 1.0e+6;

vs = np.linspace(0.0, 300.0, 100)

T = propelerThrust( vs, power, area=1.0,   )

plt.plot( vs, T       , "-b", label="quadratic" )
plt.plot( vs, power/vs, ":r" , label="naive" )
plt.legend()
plt.ylim(0.0,np.max(T))

plt.show()
