#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

import scipy.special as spc

'''

Note: It seems that H2 Molecule cannot be sable without varying Kinetic Energy

see:  

[1] https://link.aps.org/doi/10.1103/PhysRevLett.99.185003
Excited Electron Dynamics Modeling of Warm Dense Matter
Julius T. Su, William A. Goddard,

[2] http://aip.scitation.org/doi/10.1063/1.3272671
The dynamics of highly excited electronic systems: Applications of the electron force field
Julius T. Su, William A. Goddard

[3] http://dx.doi.org/10.1016/j.mechmat.2015.02.008
Non-adiabatic dynamics modeling framework for materials in extreme conditions
Hai Xiao, Andres Jaramillo-Botero, Patrick L. Theofanis, William A. Goddard


To check and obtain constants:

https://en.wikipedia.org/wiki/Hydrogen_atom#Bohr%E2%80%93Sommerfeld_Model
https://en.wikipedia.org/wiki/Fine-structure_constant

'''


# ==== constants in  SI Units


# see https://en.wikipedia.org/wiki/Fine-structure_constant

const_hbar  = 1.054571817e-34 # [J.s]  #6.582119569e-16 # [eV/s]
const_Me    = 9.10938356e-31  # [kg]
const_e     = 1.602176620898e-19  # [Coulomb]
const_eps0  = 8.854187812813e-12 # [F.m = Coulomb/(Volt*m)]
const_eV    = 1.602176620898e-19 # [J]
const_Angstroem = 1.0e-10 

const_K     =  const_hbar**2/const_Me
const_El    =  const_e**2/(4.*np.pi*const_eps0)

const_Ry     = 0.5 * const_El**2/const_K
const_Ry_eV  = 13.6056925944
const_El_eVA = const_El/( const_e*const_Angstroem )

const_K_eVA  = (const_El_eVA**2)/(2*const_Ry_eV)

print "const_El, const_El_eVA ", const_El, const_El_eVA
print "const_Ry const_Ry_eV ", const_Ry, const_Ry/const_eV
print "const_K, const_K_eVA ", const_K, const_K_eVA

#exit()

#const_K  =   const_hbar**2/const_Me   #   [ eV * A^2 ]
#const_K  = 0.1* 30.0824137226  # [eV*A^2] hbar[J.s]^2/(Me [kg])   /  (  eV[J]*A^2[m])    # (6.62607015e-34^2/9.10938356e-31)/1.602176620898e-19/10e-20
#const_Ke =  1.5*const_K

const_Ke_eVA = const_K_eVA*1.5
print "const_Ke_eVA ", const_Ke_eVA

#const_El =  14. # 14 (1./((4*np.pi*const_eps0))

sqrt2 = np.sqrt(2.)

# ================ Kinetic

def Kinetic( s ):
    '''
    Ek = (hbar^2/me) (3./2.) 1/s^2
    '''
    return const_Ke_eVA/(s**2)

def Kinetic_ds( s ):
    '''
    Ek = (hbar^2/me) (3./2.) 1/s^2
    '''
    return -2.*const_Ke_eVA/(s**3)

# ================ Electrostatic

def El( r, qq, si, sj ):
    s = np.sqrt( si**2 + sj**2 )
    #print "El ", s.shape
    #print "El ",  const_El_eVA * (qq/r) * spc.erf( sqrt2 * r/s )
    return const_El_eVA * (qq/r) * spc.erf( sqrt2 * r/s )

def El_dr( r, qq, si, sj ):
    s2 = si**2 + sj**2
    s  = np.sqrt(s2)
    e1 = (qq/r)
    e2 = spc.erf( sqrt2 * r/s )
    f1 = -qq/(r*r)
    f2 = np.exp( -2.*r*r/s2 ) *(2.*np.sqrt(2./np.pi))/s
    #print "El_dr ", -const_El_eVA * ( f1*e2 + e1*f2 )
    return const_El_eVA * ( f1*e2 + e1*f2 )

def El_ds( r, qq, si, sj ):
    s2 = si**2 + sj**2
    s  = np.sqrt(s2) 
    f  = const_El_eVA * (qq/r) * np.exp( -2.*r*r/s2 ) * (-2.*np.sqrt(2./np.pi)) * r/s2
    #print "El_ds ", s2.shape
    return f*(si/s),  f*(sj/s)

'''
def El_aa( r, qq ):
    return const_El_eVA * (qq/r)

def El_ae( r, qq, s ):
    return const_El_eVA * (qq/r) * spc.erf( sqrt2 * r/s )

def El_ee( r, qq, si, sj ):
    s = np.sqrt( si**2 + sj**2 )
    return const_El_eVA * (qq/r) * spc.erf( sqrt2 * r/s )
'''

# ================ delta-T

def getT( r, si, sj ):
    #print "getT si, sj ", si, sj
    #   r = r * 1.125
    #   s = s*0.9
    si2 = si**2
    sj2 = sj**2
    r2  = r**2

    s2 = si2 + sj2
    s4 = s2*s2

    #return const_K * ( 1.5*( 1./si2 + 1./sj2 ) )   - 2.*( 3.*(si2+sj2) - 2.*r2 )/( si2 + sj2 )**2 )
    return const_K_eVA * ( 1.5*s2/(si2*sj2) - 2.*( 3.*s2 - 2.*r2 )/s4 )
    #return const_K_eVA * 1.5*( s2/(si2*sj2) )
    #return const_K_eVA * -2. *( 3.*s2 - 2.*r2 )/s4

def getT_ds( r, si, sj ):
    #print "getT si, sj ", si, sj
    #   r = r * 1.125
    #   s = s*0.9
    si2 = si**2
    sj2 = sj**2
    r2  = r**2

    s2  = si2 + sj2
    b  =  4.*( 3.*s2 - 4.*r2 )/(s2*s2*s2)
    return const_K_eVA*( -3./(si2*si) + b*si ), const_K_eVA*( -3./(sj2*sj) + b*sj )

def getT_dr( r, si, sj ):
    #print "getT si, sj ", si, sj
    #   r = r * 1.125
    #   s = s*0.9
    si2 = si**2
    sj2 = sj**2
    s2 = si2 + sj2
    s4 = s2*s2
    #r2  = r**2
    return const_K_eVA * (8./s4) * r

# ================ delta-S

def getS( r, si, sj ):
    #print "getS si, sj ", si, sj
    #   r = r * 1.125
    #   s = s*0.9
    si2 = si**2
    sj2 = sj**2
    s2  = si2 + sj2
    r2  = r**2
    return ( 2.*(si*sj)/s2 )**1.5 * np.exp( -r2/s2 )

def getS_ds( r, si, sj ):
    #print "getS si, sj ", si, sj
    #   r = r * 1.125
    #   s = s*0.9
    si2 = si**2
    sj2 = sj**2
    s2  = si2 + sj2
    s4  = s2*s2
    r2  = r**2
    # E = ( 2.*(si*sj)/s2 )**1.5 * np.exp( -r2/s2 )

    e1 = ( 2.*(si*sj)/s2 )**1.5
    e2 = np.exp( -r2/s2 )
    f1 = 3.*( 2.*(si*sj)/s2 )**0.5 * ( si2 - sj2 )/s4 # * sj
    f2 = np.exp( -r2/s2 ) * (2.*r2/s4)               # * si
    #f  = e1*f2 + e2*f1
    return  e1*f2*si - e2*f1*sj, e1*f2*sj + e2*f1*si

def getS_dr( r, si, sj ):
    #print "getS si, sj ", si, sj
    #   r = r * 1.125
    #   s = s*0.9
    si2 = si**2
    sj2 = sj**2
    s2  = si2 + sj2
    r2  = r**2
    # E = ( 2.*(si*sj)/s2 )**1.5 * np.exp( -r2/s2 )
    e1 = ( 2.*(si*sj)/s2 )**1.5
    f2 = np.exp( -r2/s2 ) * (-2./s2) * r         # * si
    return  e1*f2

# ================ Pauli

def EPaul_S2( S2, anti=False, rho=0.2 ):
    if anti:
        return rho*S2/(1.+S2)
    else:
        return S2/(1.-S2) + ( 1.-rho )*(S2/(1.+S2))

def EPaul_dS( S, anti=False, rho=0.2 ):
    S2 = S*S
    if anti:
        #e = rho*S2/(1.+S2)
        return 2.*rho*S/(1.+S2)**2
    else:
        #e =  (S2/(1.-S2)) + (( 1.-rho )*(S2/(1.+S2))
        return 2.*S/(1-S2)**2  + ( 1.-rho )*2.*S/(1.+S2)**2 

def EPauli( r, si, sj, anti=False, rho=0.2, kr=1.125, ks=0.9 ):
    r  = r*kr
    si = si*ks
    sj = sj*ks
    T = getT( r, si, sj )
    S = getS( r, si, sj )
    S2=S*S
    return T * EPaul_S2( S2, anti=anti, rho=0.2 )
    #if anti:
    #    return T * ( rho*S2/(1.+S2) )
    #else:
    #    return T * ( (S2/(1.-S2))   + ( 1.-rho )*(S2/(1.+S2)) )

def EPauli_ds( r, si, sj, anti=False, rho=0.2, kr=1.125, ks=0.9 ):
    r  = r *kr
    si = si*ks
    sj = sj*ks
    T  = getT( r, si, sj )
    S  = getS( r, si, sj )
    dTi,dTj = getT_ds( r, si, sj )
    dSi,dSj = getS_ds( r, si, sj )
    ES  = EPaul_S2( S*S, anti=anti, rho=rho )*ks
    dES = EPaul_dS(   S, anti=anti, rho=rho )*ks
    return dTi*ES + T*dES*dSi,   dTj*ES + T*dES*dSj

def EPauli_dr( r, si, sj, anti=False, rho=0.2, kr=1.125, ks=0.9 ):
    r  = r *kr
    si = si*ks
    sj = sj*ks
    T  = getT( r, si, sj )
    S  = getS( r, si, sj )
    dT = getT_dr( r, si, sj )
    dS = getS_dr( r, si, sj )
    ES  = EPaul_S2( S*S, anti=anti, rho=rho )*kr
    dES = EPaul_dS(   S, anti=anti, rho=rho )*kr
    return dT*ES + T*dES*dS

def numDeriv( x, y ):
    dy =  (y[2:]-y[:-2]) / (x[2:]-x[:-2])
    print dy.shape, x[1:-1].shape 
    return dy,x[1:-1]

def plotVsNum( x,y,ynum, name ):
    xnum=x[1:-1]
    plt.plot(x, y ,'-', label=name )
    plt.plot(xnum, ynum,':',               label=name+"_num" )
    #plt.plot(xnum, abs(y[1:-1]-ynum), '--', label=name+"_err" )

def checkNumDeriv( x, func, dfunc, name ):
    dy         = dfunc( x )
    y          = func(x)
    dynum,xnum = numDeriv( x, y )
    #print y
    #print "y.shape, ynum.shape ", y.shape, ynum.shape
    plotVsNum( x,dy,dynum, name )
    plt.plot(x, y,'-.', label=name+"_F" )

if __name__ == "__main__":

    #extent=( 0.5,8.0,  0.5,4.5 )
    #xs = np.arange( extent[0], extent[1], 0.05 )
    #ys = np.arange( extent[2], extent[3], 0.1  )

    # ==== With respect to blob-size  's'

    s  = np.arange( 0.1, 5.0, 0.05 )
    rs = np.arange( 0.1, 5.0, 0.05 )
    S  = np.arange( 1.25, 5.0, 0.05 )
    r  = 1.5 + 0.*s
    #sj = 0.7 + 0.*s
    sj = 0.7
    si = 0.58
    qq = 1.0

    #checkNumDeriv( s, Kinetic, Kinetic_ds, "dKinetic" )
    #checkNumDeriv( s, lambda x: El(r,qq,x,sj), lambda x : El_ds(r,qq,x,sj)[0], "dEl_ds" )
    #checkNumDeriv( s, lambda x: getT(r,x,sj), lambda x : getT_ds(r,x,sj)[0], "dT_ds" )
    #checkNumDeriv( s, lambda x: getS(r,x,sj), lambda x : getS_ds(r,x,sj)[0], "dS_ds" )
    #checkNumDeriv( S, lambda x: EPaul_S2(x*x,anti=True), lambda x : EPaul_dS(x,anti=True), "dES_ds" )
    #checkNumDeriv( S, lambda x: EPaul_S2(x*x,anti=False), lambda x : EPaul_dS(x,anti=False), "dES_ds" )
    #checkNumDeriv( s, lambda x: EPauli( r, x, sj, anti=False), lambda x : EPauli_ds( r, x, sj, anti=False)[0], "dPauli66_ds" )
    #checkNumDeriv( s, lambda x: EPauli( r, x, sj, anti=True), lambda x : EPauli_ds( r, x, sj, anti=True)[0], "dPauli69_ds" )

    #checkNumDeriv( rs, lambda x: El(x,qq,si,sj), lambda x : El_dr(x,qq,si,sj), "dEl_dr" )
    #checkNumDeriv( rs, lambda x: getT(x,si,sj), lambda x : getT_dr(x,si,sj), "dT_dr" )
    #checkNumDeriv( rs, lambda x: getS(x,si,sj), lambda x : getS_dr(x,si,sj), "dS_dr" )

    checkNumDeriv( rs, lambda x: EPauli( x, si, sj, anti=False), lambda x : EPauli_dr( x, si, sj, anti=False), "dPauli66++_dr" )
    checkNumDeriv( rs, lambda x: EPauli( x, si, sj, anti=True), lambda x : EPauli_dr( x, si, sj, anti=True), "dPauli69_dr" )

    #ylim=100.0; plt.ylim(-ylim,ylim)

    plt.grid()
    plt.legend()
    plt.show()

