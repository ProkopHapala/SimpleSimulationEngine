#!/usr/bin/env python

from pylab import *
from Poly6th_numeric_simple import *
figure(num=None, figsize=(14, 6))


def gravity(R):
	return 1.0/ R**2

def omega(R):
	return sqrt(1.0/ R**3)


xs = array([ 1, 2, 3, 4, 5 ])
#vs = array([ 1, 2, 3, 4, 5 ])
#a  = array([ 1, 2, 3, 4, 5 ])

vs = omega(xs)
a  = gravity(xs)

N = len(xs)

ts = arange(0,1.0001,0.02)

# ======= Velocity conditions =========

for i in range(len(xs)-1):
    As = poly6Coefs(   xs[i], vs[i],a[i],     xs[i+1],vs[i+1], a[i+1] )
    P,dP,ddP = evalPoly6( ts, As)
    subplot(1,3,1); plot( ts+i,   P,'-' )
    subplot(1,3,2); plot( ts+i,  dP,'-' )
    subplot(1,3,3); plot( ts+i, ddP,'-' )


show()
