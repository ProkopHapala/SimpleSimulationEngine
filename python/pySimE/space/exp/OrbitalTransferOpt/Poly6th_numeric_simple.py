from pylab import *

# =====  6th order polynom of time which fit position, velocty, and acceleration and end points

def poly6Coefs(p0,v0,a0, p1,v1,a1):
    As = zeros(6)
    As[0] =                                          p0
    As[1] =                                          v0
    As[2] =                                      0.5*a0
    As[3] = -1.5*a0 + 0.5*a1 - 10*p0 + 10*p1 - 6*v0 - 4*v1
    As[4] =  1.5*a0 -     a1 + 15*p0 - 15*p1 + 8*v0 + 7*v1
    As[5] = -0.5*a0 + 0.5*a1 -  6*p0 +  6*p1 - 3*v0 - 3*v1
    return As

def evalPoly6( t, As):
    x   = As[0] + As[1]*t +    As[2]*t**2 +    As[3]*t**3 +    As[4]*t**4 +    As[5]*t**5
    dx  =         As[1]   +  2*As[2]*t    +  3*As[3]*t**2 +  4*As[4]*t**3 +  5*As[5]*t**4
    ddx =                    2*As[2]      +  6*As[3]*t    + 12*As[4]*t**2 + 20*As[5]*t**3
    return x,dx,ddx
