# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 22:50:34 2013

@author: asiJa
"""

import numpy as np
from pylab import *

# copied from java script  CoffeeScript 
# http://blog.mackerron.com/2011/01/01/javascript-cubic-splines/
def Spline4_clamped(x, a, d0, dn):
    np1 = len(x)
    n = np1 - 1
    #print n
    h = [0.0]*np1; y = [0.0]*np1; l = [0.0]*np1; u = [0.0]*np1; z = [0.0]*np1; 
    c = [0.0]*np1; b = [0.0]*np1; d = [0.0]*np1; k = [0.0]*np1; s = [0.0]*np1;
    for i in range(0,n):
        h[i] = x[i + 1] - x[i]
        k[i] = a[i + 1] - a[i]
        s[i] = k[i] / h[i]
    #if clamped
    y[0] = 3 * (a[1] - a[0]) / h[0] - 3 * d0
    y[n] = 3 * dn - 3 * (a[n] - a[n - 1]) / h[n - 1]
    for i in range(1,n):
        y[i] = 3 / h[i] * (a[i + 1] - a[i]) - 3 / h[i - 1] * (a[i] - a[i - 1])
    #if clamped
    l[0] = 2 * h[0]
    u[0] = 0.5
    z[0] = y[0] / l[0]
    for i in range(1,n):
        l[i] = 2 * (x[i + 1] - x[i - 1]) - h[i - 1] * u[i - 1]
        u[i] = h[i] / l[i]
        z[i] = (y[i] - h[i - 1] * z[i - 1]) / l[i]
    #if clamped
    l[n] = h[n - 1] * (2 - u[n - 1])
    z[n] = (y[n] - h[n - 1] * z[n - 1]) / l[n]
    c[n] = z[n]
    for i in xrange(n-1, -1, -1):
        c[i] = z[i] - u[i] * c[i + 1]
        b[i] = (a[i + 1] - a[i]) / h[i] - h[i] * (c[i + 1] + 2 * c[i]) / 3
        d[i] = (c[i + 1] - c[i]) / (3 * h[i])
    #return x,a,b,c,d
    splines = []
    for i in xrange(n+1):
        splines.append((a[i],b[i],c[i],d[i],x[i]))
    #print " len(x) ",len(x)," len(splines) ",len(splines)
    return splines,a[n]


def evalSpline4( xn, splines, perSpline):
    #print len(splines)
    n = len(splines)-1
    m = n*perSpline
    X= zeros(m); Y= zeros(m); dY= zeros(m); ddY= zeros(m);
    uts = arange(0,1.0, 1.0/perSpline)
    for i in range(n):
        S  = splines[i]
        t  = uts*(splines[i+1][4]-S[4])
        t2 = t**2
        t3 = t2*t
        j0=perSpline*i
        j1=perSpline*(i+1)
        #print j0,j1,j0-j1, len(t), len(Y[j0:j1])
        X  [j0:j1]=t+S[4]
        Y  [j0:j1]=  S[0]+S[1]*t +   S[2]*t2 +   S[3]*t3
        dY [j0:j1]=       S[1]   + 2*S[2]*t  + 3*S[3]*t2
        ddY[j0:j1]=                2*S[2]    + 6*S[3]*t
    #print " t : ",t
    #print " X : ",X
    
    return X,Y,dY,ddY

def evalPolarForces( T, Rs, Os,  perSpline ):
    ts,O,dO,ddO = evalSpline4(T,Os,perSpline)
    ts,R,dR,ddR = evalSpline4(T,Rs,perSpline)
    #print ts
    Os = [ O, dO, ddO ]
    Rs = [ R, dR, ddR ]
    # Gravitational force
    R2 = R*R
    G   =  -1.0 / R2
    # Total force
    F_O  = R * ddO  + 2 * dR * dO         # Angular Kinematic Force = Angular engine thrust 
    F_R  = ddR      - R *  dO**2               # Radial Kinematic force
    FTR  = F_R - G                           # Radial engine thrust
    FT   = sqrt( F_O**2 + FTR**2 )           # Total  engine Trust Force
    Fs = [F_O,F_R,G,FTR, FT]
    return ts, Os, Rs, Fs