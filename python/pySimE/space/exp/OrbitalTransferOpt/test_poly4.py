#!/usr/bin/env python

from pylab import *
from Poly4th_numeric import *
figure(num=None, figsize=(14, 6))

xs = [ 1, 2, -1, 4, 5 ]
vs = [ 1, 2, 3, 4, 5 ]
a  = [ 1, 2, 3, 4, 5 ]
N = len(xs)

ts = arange(0,1.0001,0.02)

# ======= Velocity conditions =========

for i in range(len(xs)-1):
    As = Coefs_x0x1v0v1(xs[i], xs[i+1], vs[i], vs[i+1])
    P,dP,ddP = evalPoly4(ts,As)
    subplot(1,3,1); plot( ts+i,   P,'-' )
    subplot(1,3,2); plot( ts+i,  dP,'-' )
    subplot(1,3,3); plot( ts+i, ddP,'-' )


# ======= Acceleration conditions =========

Ass = []
Ass.append(     Coefs_x0x1v0a1(xs[0]  , xs[1]  , vs[0] , a[1])   )
for i in range(1,N-2):
    Ass.append( Coefs_x0x1a0a1(xs[i]  , xs[i+1], a[i]  , a[i+1]) )
Ass.append(     Coefs_x0x1a0v1(xs[N-2], xs[N-1], a[N-2], vs[N-1])   )
for i in range(N-1):
    P,dP,ddP = evalPoly4(ts,Ass[i])
    subplot(1,3,1); plot( ts+i,   P,'--' )
    subplot(1,3,2); plot( ts+i,  dP,'--' )
    subplot(1,3,3); plot( ts+i, ddP,'--' )


show()
