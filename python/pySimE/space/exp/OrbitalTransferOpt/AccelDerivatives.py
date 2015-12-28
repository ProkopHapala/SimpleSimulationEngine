from pylab import *
from Poly6th_numeric import *



def diff_P_a0(t):
    return -t**5/2 + 3*t**4/2 - 3*t**3/2 + t**2/2

def diff_P_a1(t):
    return t**5/2 - t**4 + t**3/2

def diff_dP_a0(t):
    return -5*t**4/2 + 6*t**3 - 9*t**2/2 + t

def diff_dP_a1(t):
    return 5*t**4/2 - 4*t**3 + 3*t**2/2

def diff_ddP_a0(t):
    return -10*t**3 + 18*t**2 - 9*t + 1

def diff_ddP_a1(t):
    return 10*t**3 - 12*t**2 + 3*t


def dF_Phi_ap(R, dR,   dP, ddP ):
    # R(ar0, ar1)*Derivative(ddPhi(ap0, ap1), ap0) + 2*dR(ar0, ar1)*Derivative(dPhi(ap0, ap1), ap0)
    return R*ddP + 2*dR*dP

def dF_Phi_ar(Phi, dPhi,  P,  dP):
    # 2*dPhi(ap0, ap1)*Derivative(dR(ar0, ar1), ar0) + ddPhi(ap0, ap1)*Derivative(R(ar0, ar1), ar0)
    return 2*dPhi*dP + ddPhi*P

def dF_R_ap(R, dR,   dP, ddP ):
    # 2*R(ar0, ar1)*dPhi(ap0, ap1)*Derivative(dPhi(ap0, ap1), ap0)
    2*R*dPhi*dP

def dF_R_ar(R, dR,   dP, ddP ):
    # dPhi(ap0, ap1)**2*Derivative(R(ar0, ar1), ar0) + Derivative(ddR(ar0, ar1), ar0)
    dPhi**2*P + ddP



def evalAll(t,Rs, Os):
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
    # Variational derivatives by parameter a0 and a1 
    P_a0   = -0.5*t5 + 1.5*t4 - 1.5*t3   + 0.5*t2
    P_a1   =  0.5*t5 -     t4 + 0.5*t3        
    dP_a0  = -2.5*t4 +   6*t3 - 4.5*t2/2 +     t
    dP_a1  =  2.5*t4 -   4*t3 + 1.5*t2
    ddP_a0 =  -10*t3 +  18*t2 -   9*t    +     1
    ddP_a1 =   10*t3 -  12*t2 +   3*t
    # Variation of Kinematic force in polar coordinates  
    dRtwo = 2*dR    
    dF_O_ao0 = R*ddP_a0 + dRtwo*dP_a0
    dF_O_ao1 = R*ddP_a1 + dRtwo*dP_a1
    dOtwo = 2*dO
    dF_O_ar0 = dOtwo*dP_a0 + ddO*P_a0
    dF_O_ar1 = dOtwo*dP_a1 + ddO*P_a1
    RdO = 2*R*dO
    dF_R_ao0   = RdO*dP_a0
    dF_R_ao1   = RdO*dP_a1
    dO2 = dO*dO
    dF_R_ar0   = dO2*P_a0 + ddP_a0
    dF_R_ar1   = dO2*P_a1 + ddP_a1
    # Gravitational force
    R2 = R*R
    G   = -1 / R2
    R3 = R2*R
    dG_a0 = -2*P_a0 / R3
    dG_a1 = -2*P_a1 / R3
    # Total force
    F_O  = R * ddO     + 2 * dR * dO         # Angular Kinematic Force = Angular engine thrust 
    F_R  = R *  dO**2  +    ddR              # Radial Kinematic force
    FTR  = F_R - G                           # Radial engine thrust
    FT   = sqrt( F_O**2 + FTR**2 )           # Total  engine Trust Force
    # Variation of total force
    dF_ao0 =  (  F_O*dF_O_ao0 + FTR*  dF_R_ao0           ) / FT
    dF_ao1 =  (  F_O*dF_O_ao1 + FTR*  dF_R_ao1           ) / FT
    dF_ar0 =  (  F_O*dF_O_ar0 + FTR*( dF_R_ar0 - dG_a0 ) ) / FT
    dF_ar1 =  (  F_O*dF_O_ar1 + FTR*( dF_R_ar1 - dG_a1 ) ) / FT




'''
def dFPhi_da0(t,    r0, r1, dr0, dr1, ddr0,ddr1 ):
    T1 = -t*(5*t**3 - 12*t**2 + 9*t - 2)
    L2 =  (-5*ddr0 + 5*ddr1 - 30*dr0 - 30*dr1 - 60*r0 + 60*r1)
    L3 =  (12*ddr0 - 8*ddr1 + 64*dr0 + 56*dr1 + 120*r0 - 120*r1)
    L4 =  (-9*ddr0 + 3*ddr1 - 36*dr0 - 24*dr1 - 60*r0 + 60*r1) 
    T5 =  10*t**3 - 18*t**2 + 9*t - 1
    L6 =  ddr0*t**2 + 2*dr0*t + 2*r0
    L7 = -ddr0 + ddr1 - 6*dr0 - 6*dr1 - 12*r0 + 12*r1
    L8 = 3*ddr0 - 2*ddr1 + 16*dr0 + 14*dr1 + 30*r0 - 30*r1
    L9 = (-3*ddr0 + ddr1 - 12*dr0 - 8*dr1 - 20*r0 + 20*r1)
    return T1*(2*ddr0*t + 2*dr0 + t**4*L2 + t**3*L3 + t**2*L4)/2 - T5*(L6 + t**5*L7 + t**4*L8 + t**3*L9)/2

def dFPhi_da1(t,    r0, r1, dr0, dr1, ddr0,ddr1 ):
    L1 = t*(t*(5*t**2 - 8*t + 3)
    L2 = -5*ddr0 + 5*ddr1 - 30*dr0 - 30*dr1 - 60*r0 + 60*r1
    L3 = 12*ddr0 - 8*ddr1 + 64*dr0 + 56*dr1 + 120*r0 - 120*r1
    L4 = -9*ddr0 + 3*ddr1 - 36*dr0 - 24*dr1 - 60*r0 + 60*r1
    T1 = 10*t**2 - 12*t + 3
    L1*(2*ddr0*t + 2*dr0 + t**4*L2 + t**3*L3 + t**2*L4) + T1*(ddr0*t**2 + 2*dr0*t + 2*r0 + t**5*(-ddr0 + ddr1 - 6*dr0 - 6*dr1 - 12*r0 + 12*r1) + t**4*(3*ddr0 - 2*ddr1 + 16*dr0 + 14*dr1 + 30*r0 - 30*r1) + t**3*(-3*ddr0 + ddr1 - 12*dr0 - 8*dr1 - 20*r0 + 20*r1)))/2

def dFR_da0 =   -10*t**3 + t**2*(-t**3/8 + 3*t**2/8 - 3*t/8 + 1/8)*(2*ddh0*t + 2*dh0 - 5*t**4*(ddh0 - ddh1 + 6*dh0 + 6*dh1 + 12*h0 - 12*h1) + 4*t**3*(3*ddh0 - 2*ddh1 + 16*dh0 + 14*dh1 + 30*h0 - 30*h1) - 3*t**2*(3*ddh0 - ddh1 + 12*dh0 + 8*dh1 + 20*h0 - 20*h1))**2 + 18*t**2 - 9*t + 1

def dFR_da1 =   t*(t**2*(t**2 - 2*t + 1)*(2*ddh0*t + 2*dh0 - 5*t**4*(ddh0 - ddh1 + 6*dh0 + 6*dh1 + 12*h0 - 12*h1) + 4*t**3*(3*ddh0 - 2*ddh1 + 16*dh0 + 14*dh1 + 30*h0 - 30*h1) - 3*t**2*(3*ddh0 - ddh1 + 12*dh0 + 8*dh1 + 20*h0 - 20*h1))**2 + 80*t**2 - 96*t + 24)/8
'''