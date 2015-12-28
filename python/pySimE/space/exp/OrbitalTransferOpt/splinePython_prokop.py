# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 22:50:34 2013

@author: asiJa
"""

import numpy as np
from pylab import *

 
# spline defined as  
# S[i](x) = a[i] + b[i](x-x[i]) + c[i](x-x[i])**2 + d[i](x-x[i])**3 
def Splines(data):
    np1=len(data)
    n=np1-1
    X,Y = zip(*data)
    X = [float(x) for x in X]
    Y = [float(y) for y in Y]
    a = Y[:]
    b = [0.0]*(n)
    d = [0.0]*(n)
    h = [X[i+1]-X[i] for i in xrange(n)]   #  h[i]
    alpha = [0.0]*n
    for i in xrange(1,n):
        alpha[i] = 3/h[i]*(a[i+1]-a[i]) - 3/h[i-1]*(a[i]-a[i-1])   # right had side
    c = [0.0]*np1
    L = [0.0]*np1
    u = [0.0]*np1
    z = [0.0]*np1
    L[0] = 1.0; u[0] = z[0] = 0.0   # S''(x[0]) = 0
    for i in xrange(1,n):
        L[i] = 2*(X[i+1]-X[i-1]) - h[i-1]*u[i-1]
        u[i] = h[i]/L[i]
        z[i] = (alpha[i]-h[i-1]*z[i-1])/L[i]
    L[n] = 1.0; z[n] = c[n] = 0.0   # S''(x[N]) = 0
    for j in xrange(n-1, -1, -1):
        c[j] = z[j] - u[j]*c[j+1]
        b[j] = (a[j+1]-a[j])/h[j] - (h[j]*(c[j+1]+2*c[j]))/3
        d[j] = (c[j+1]-c[j])/(3*h[j])
    splines = []
    for i in xrange(n):
        splines.append((a[i],b[i],c[i],d[i],X[i]))
    return splines,X[n]

# copied from java script  CoffeeScript 
# http://blog.mackerron.com/2011/01/01/javascript-cubic-splines/
def Splines_clamped(x, a, d0, dn):
    np1 = len(x)
    n = np1 - 1
    print n
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
    for i in xrange(n):
        splines.append((a[i],b[i],c[i],d[i],x[i]))
    return splines,a[n]
 
def splinesToPlot(splines,xn,res):
    n=len(splines)
    perSpline = int(res/n)
    if perSpline < 3: perSpline = 3
    X=[]; Y=[]; dY=[]; ddY=[]; 
    for i in xrange(n-1):
        S = splines[i]
        x0 = S[4]
        x1 = splines[i+1][4]
        x = np.linspace(x0,x1,perSpline)
        for xi in x:
            X.append(xi)
            h=(xi-S[4])
            #Y.append(S[0]+S[1]*h + S[2]*h**2 + S[3]*h**3)
            Y  .append(S[0]+S[1]*h +   S[2]*h**2 +   S[3]*h**3)
            dY .append(     S[1]   + 2*S[2]*h    + 3*S[3]*h**2)
            ddY.append(            + 2*S[2]      + 6*S[3]*h   )
    S=splines[n-1]
    x=np.linspace(S[4],xn,perSpline)
    for xi in x:
        X.append(xi)
        h=(xi-S[4])
        Y  .append(S[0]+S[1]*h +   S[2]*h**2 +   S[3]*h**3)
        dY .append(     S[1]   + 2*S[2]*h    + 3*S[3]*h**2)
        ddY.append(            + 2*S[2]      + 6*S[3]*h   )
    return X,Y,dY, ddY
 
x = lambda n: np.linspace(-1.0, 0.0,n)**3 + np.linspace(-1.0, 0.0,n)
#f = lambda x: 2*np.cos(np.sin(2*np.pi*x)) -0.8

f   = lambda x:  np.cos(2*np.pi*x)
df  = lambda x: -np.sin(2*np.pi*x)*2*np.pi
ddf = lambda x: -np.cos(2*np.pi*x)*(2*np.pi)**2

n = 12
E=200
data = zip(x(n),f(x(n)))
splines1,xn1 = Splines(data)
splines2,xn2 = Splines_clamped(x(n),f(x(n)), 0, 0)
X1,Y1,dY1,ddY1 = splinesToPlot(splines1,xn1,E)
X2,Y2,dY2,ddY2 = splinesToPlot(splines2,xn1,E)

figure(num=None, figsize=(16, 6))
subplot(1,3,1)
plot(x(300),f(x(300)),'k-',  label="Exact")
plot(x(n),f(x(n)),    'ok',  label="control points")
plot(X1,Y1,           'r--', label="Natrual spline")
plot(X2,Y2,           'b--', label="Clamped spline")
ylim(-1.5,2.0)
legend(loc=9)

subplot(1,3,2)
plot(x(300),df(x(300)),'k-',  label="Exact")
plot(X1,dY1,           'r--', label="Natrual spline")
plot(X2,dY2,           'b--', label="Clamped spline")

subplot(1,3,3)
plot(x(300),ddf(x(300)),'k-',  label="Exact")
plot(X1,ddY1,           'r--', label="Natrual spline")
plot(X2,ddY2,           'b--', label="Clamped spline")

show()