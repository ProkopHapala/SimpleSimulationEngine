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

def evalPoly6_array( ts, As,  xs,dxs, ddxs  ):
    xs  [:] = As[0] + As[1]*ts[:] +    As[2]*ts[:]**2 +    As[3]*ts[:]**3 +    As[4]*ts[:]**4 +    As[5]*ts[:]**5
    dxs [:] =         As[1]       +  2*As[2]*ts[:]    +  3*As[3]*ts[:]**2 +  4*As[4]*ts[:]**3 +  5*As[5]*ts[:]**4
    ddxs[:] =                        2*As[2]          +  6*As[3]*ts[:]    + 12*As[4]*ts[:]**2 + 20*As[5]*ts[:]**3

def ForceEstimate_Linear( p0,v0, p1, v1 ):
    a0 = -6.0*p0 + 6.0*p1 - 4.0*v0 - 2.0*v1
    a1 =  6.0*p0 - 6.0*p1 + 2.0*v0 + 4.0*v1
    return a0,a1

def ForceEstimate_Cos( p0,v0, p1, v1 ):
    a0 = -5.0*p0 + 5.0*p1 - 3.5*v0 - 1.5*v1
    a1 =  5.0*p0 - 5.0*p1 + 1.5*v0 + 3.5*v1
    return a0,a1

# ======  Forces and variations in polar coordinates

def evalPolarForces(t, T, Rs, Os):
    T2 = T*T
    # Powers of time
    t2 = t*t
    t3 = t2*t
    t4 = t3*t
    t5 = t4*t
    # Polynom and its time-derivatives
    O   =  Os[0] + Os[1]*t +    Os[2]*t2   +    Os[3]*t3   +    Os[4]*t4 +    Os[5]*t5
    dO  = (        Os[1]   +  2*Os[2]*t    +  3*Os[3]*t2   +  4*Os[4]*t3 +  5*Os[5]*t4 ) /T
    ddO = (                   2*Os[2]      +  6*Os[3]*t    + 12*Os[4]*t2 + 20*Os[5]*t3 ) /T2
    R   =  Rs[0] + Rs[1]*t +    Rs[2]*t2   +    Rs[3]*t3   +    Rs[4]*t4 +    Rs[5]*t5
    dR  = (        Rs[1]   +  2*Rs[2]*t    +  3*Rs[3]*t2   +  4*Rs[4]*t3 +  5*Rs[5]*t4 ) /T
    ddR = (                   2*Rs[2]      +  6*Rs[3]*t    + 12*Rs[4]*t2 + 20*Rs[5]*t3 ) /T2
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


# FIXME - What to do with time scaling? 1/T   1/T*T ????????????????
def evalPolarForceAndVariations(t,Rs, Os):
    # Powers of time
    t2 = t*t
    t3 = t2*t
    t4 = t3*t
    t5 = t4*t
    # Polynom and its time-derivatives
    O   = Os[0] + Os[1]*t +    Os[2]*t2   +    Os[3]*t3   +    Os[4]*t4 +    Os[5]*t5
    dO  =         Os[1]   +  2*Os[2]*t    +  3*Os[3]*t2   +  4*Os[4]*t3 +  5*Os[5]*t4
    ddO =                    2*Os[2]      +  6*Os[3]*t    + 12*Os[4]*t2 + 20*Os[5]*t3
    R   = Rs[0] + Rs[1]*t +    Rs[2]*t2   +    Rs[3]*t3   +    Rs[4]*t4 +    Rs[5]*t5
    dR  =         Rs[1]   +  2*Rs[2]*t    +  3*Rs[3]*t2   +  4*Rs[4]*t3 +  5*Rs[5]*t4
    ddR =                    2*Rs[2]      +  6*Rs[3]*t    + 12*Rs[4]*t2 + 20*Rs[5]*t3
    Os = [ O, dO, ddO ]
    Rs = [ R, dR, ddR ]
    # Gravitational force
    R2 = R*R
    G   = 1 / R2
    # Total force
    F_O  = R * ddO     + 2 * dR * dO         # Angular Kinematic Force = Angular engine thrust 
    F_R  = R *  dO**2  +    ddR              # Radial Kinematic force
    FTR  = F_R - G                           # Radial engine thrust
    FT   = sqrt( F_O**2 + FTR**2 )           # Total  engine Trust Force
    Fs = [F_O,F_R,G,FTR, FT]
    # Variational derivatives of polynoms by parameter a0 and a1 
    P_a0   = -0.5*t5 + 1.5*t4 - 1.5*t3   + 0.5*t2
    P_a1   =  0.5*t5 -     t4 + 0.5*t3        
    dP_a0  = -2.5*t4 +   6*t3 - 4.5*t2/2 +     t
    dP_a1  =  2.5*t4 -   4*t3 + 1.5*t2
    ddP_a0 =  -10*t3 +  18*t2 -   9*t    +     1
    ddP_a1 =   10*t3 -  12*t2 +   3*t
    # Variation of Kinematic force in polar coordinates  
    dRtwo = 2*dR    
    dF_O_ao0   = R*ddP_a0 + dRtwo*dP_a0
    dF_O_ao1   = R*ddP_a1 + dRtwo*dP_a1
    dOtwo = 2*dO
    dF_O_ar0   = dOtwo*dP_a0 + ddO*P_a0
    dF_O_ar1   = dOtwo*dP_a1 + ddO*P_a1
    RdO = 2*R*dO
    dF_R_ao0   = RdO*dP_a0
    dF_R_ao1   = RdO*dP_a1
    dO2 = dO*dO
    dF_R_ar0   = dO2*P_a0 + ddP_a0
    dF_R_ar1   = dO2*P_a1 + ddP_a1
    dFRPs = [ [[dF_O_ao0, dF_O_ao1],[dF_O_ar0,dF_O_ar1]],[[dF_R_ao0,dF_R_ao1][dF_R_ar0,dF_R_ar1]] ] 
    # Graviatioanl force Variations of total force
    R3 = R2*R
    dG_a0 = -2*P_a0 / R3
    dG_a1 = -2*P_a1 / R3
    dGs = [dG_a0,dG_a1]
    # Total Force Variations
    dFT_ao0 =  (  F_O*dF_O_ao0 + FTR*  dF_R_ao0           ) / FT
    dFT_ao1 =  (  F_O*dF_O_ao1 + FTR*  dF_R_ao1           ) / FT
    dFT_ar0 =  (  F_O*dF_O_ar0 + FTR*( dF_R_ar0 - dG_a0 ) ) / FT
    dFT_ar1 =  (  F_O*dF_O_ar1 + FTR*( dF_R_ar1 - dG_a1 ) ) / FT
    dFTs = [[dFT_ao0,dFT_ao1],[dFT_ar0,dFT_ar1 ]]
    return Os,Rs,Fs,dFRPs,dGs,dFTs