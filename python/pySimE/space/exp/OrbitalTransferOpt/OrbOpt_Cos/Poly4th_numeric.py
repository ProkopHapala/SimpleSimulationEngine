from pylab import *

def Coefs_x0x1v0v1(x0, x1, v0, v1):
    As = zeros(4)
    As[0] =    x0
    As[1] =                v0
    As[2] = -3*x0 +3*x1 -2*v0 - v1  
    As[3] = +2*x0 -2*x1 +  v0 + v1
    return As
    
def Coefs_x0x1a0a1(x0, x1, a0, a1):
    As = zeros(4)
    As[0] =  x0
    As[1] = -x0 +x1   -a0/3.0 - a1/6.0 
    As[2] =            a0/2.0  
    As[3] =           -a0/6.0 + a1/6.0
    return As
    
def Coefs_x0x1v0a1(x0, x1, v0, a1):
    As = zeros(4)
    As[0] =    x0
    As[1] =                        v0
    As[2] = -1.5*(x0 -x1 +v0) -a1/4.0     
    As[3] = +0.5*(x0 -x1 +v0) +a1/4.0         
    return As
    
def Coefs_x0x1a0v1(x0, x1, a0, v1):
    As = zeros(4)
    As[0] =       x0
    As[1] = -1.5*(x0 - x1) -0.5*v1/2 -a0/4.0 
    As[2] =                           a0/2.0    
    As[3] =  0.5*(x0 - x1 +     v1)  -a0/4.0          
    return As  

def evalPoly4( t, As):
    x   = As[0] + As[1]*t +    As[2]*t**2 +    As[3]*t**3
    dx  =         As[1]   +  2*As[2]*t    +  3*As[3]*t**2
    ddx =                    2*As[2]      +  6*As[3]*t
    return x,dx,ddx
    
def evalPolarForces(t, T, Rs, Os):
    T2 = T*T
    # Powers of time
    t2 = t*t
    t3 = t2*t
    # Polynom and its time-derivatives
    O   =  Os[0] + Os[1]*t +    Os[2]*t2   +    Os[3]*t3   
    dO  = (        Os[1]   +  2*Os[2]*t    +  3*Os[3]*t2 ) /T
    ddO = (                   2*Os[2]      +  6*Os[3]*t  ) /T2
    R   =  Rs[0] + Rs[1]*t +    Rs[2]*t2   +    Rs[3]*t3   
    dR  = (        Rs[1]   +  2*Rs[2]*t    +  3*Rs[3]*t2 ) /T
    ddR = (                   2*Rs[2]      +  6*Rs[3]*t  ) /T2
    Os = [ O, dO, ddO ]
    Rs = [ R, dR, ddR ]
    # Gravitational force
    R2 = R*R
    G   =  1 / R2
    # Total force
    F_O  = R * ddO     + 2 * dR * dO         # Angular Kinematic Force = Angular engine thrust 
    F_R  = R *  dO**2  +    ddR              # Radial Kinematic force
    FTR  = F_R - G                           # Radial engine thrust
    FT   = sqrt( F_O**2 + FTR**2 )           # Total  engine Trust Force
    Fs = [F_O,F_R,G,FTR, FT]
    return Os, Rs, Fs