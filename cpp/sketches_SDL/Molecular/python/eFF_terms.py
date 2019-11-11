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

'''

const_hbar  = 6.62607015e-34 # [J.s]  #6.582119569e-16 # [eV/s]
const_Me    = 9.10938356e-31 # [kg]
const_eps0  = 8.854187812813e-12 # [F.m]

const_eV = 1.602176620898e-19 # [J]

#const_K  =   const_hbar**2/const_Me   #   [ eV * A^2 ]
const_K  = 0.1* 30.0824137226  # [eV*A^2] hbar[J.s]^2/(Me [kg])   /  (  eV[J]*A^2[m])    # (6.62607015e-34^2/9.10938356e-31)/1.602176620898e-19/10e-20
const_Ke =  1.5*const_K
const_El =  14. # 14 (1./((4*np.pi*const_eps0))

sqrt2 = np.sqrt(2.)

def Kinetic( s ):
    '''
    Ek = (hbar^2/me) (3./2.) 1/s^2
    '''
    return const_Ke/(s**2)


def El( r, qq, si=0, sj=0 ):
    if si>0:
        if sj>0:
            si = np.sqrt( si**2 + sj**2 )
        return const_El * (qq/r) * spc.erf( sqrt2 * r/s )
    else:
        return const_El * (qq/r)

def El_aa( r, qq ):
    return const_El * (qq/r)

def El_ae( r, qq, s ):
    return const_El * (qq/r) * spc.erf( sqrt2 * r/s )

def El_ee( r, qq, si, sj ):
    s = np.sqrt( si**2 + sj**2 )
    return const_El * (qq/r) * spc.erf( sqrt2 * r/s )

def getT( r, si, sj ):
    #   r = r * 1.125
    #   s = s*0.9
    si2 = si**2
    sj2 = sj**2
    r2  = r**2
    #return const_K * ( 1.5*( (si2+sj2)/(si2*sj2) )   - 2.*( 3.*(si2+sj2) - 2.*r2 )/( si2 + sj2 )**2 )
    return const_K * ( 1.5*( 1./si2 + 1./sj2 )   - 2.*( 3.*(si2+sj2) - 2.*r2 )/( si2 + sj2 )**2 )

def getS( r, si, sj ):
    #   r = r * 1.125
    #   s = s*0.9
    si2 = si**2
    sj2 = sj**2
    r2  = r**2
    return ( 2.*(si*sj)/(si2-sj2) )**1.5 * np.exp(  -r2/( si2 + sj2 ) )

def EPauli( r, qq, si, sj, rho=0.2 ):
    T = getT( r, si, sj )
    S = getS( r, si, sj )
    S2 = S**2
    # ( S2*(1+S2) + (1-rho)* S2*(1-S2) )  / (1-S2*S2 )
    # ( S2+S2*S2 + (1-rho)*(S2-S2*S2) )  / (1-S2*S2 )
    # ( ( (2-rho)*S2 +rho*S2*S2 )  / (1-S2*S2 )
    return T * ( (S2/(1.-S2))   + ( 1-rho )*(S2/(1+S2))     )

def EPauli_pair( r, qq, si, sj, rho=0.2 ):
    T  = getT( r, si, sj )
    S  = getS( r, si, sj )
    S2 = S**2
    return T * ( rho*S2/(1+S2) )


if __name__ == "__main__":
    xs = np.arange( 0.5, 6.0, 0.05 )
    ys = np.arange( 0.5, 2.5, 0.1  )

    Xs,Ys = np.meshgrid( xs,ys )

    #s = 1.5

    y = 0.1

    Ek = Kinetic( Ys )
    Eaa = El_aa( Xs*2., 1. )
    Eae = El_ae( np.sqrt(Xs**2+y**2), -1., Ys )
    Eee = El_ee( 2.*y, 1., Ys, Ys )

    #Eaa = El_aa( Xs*2., 1. )
    #Eae = El_ae( np.sqrt(Xs**2+y**2), -1., Ys )
    #Eee = El_ee( 2.*y, 1., Ys, Ys )

    Etot = Ek + Eaa + Eae

    print Eae.min(), Etot.min()
    vmin=Etot.min()
    #vmin=Eae.min()

    plt.figure(figsize=(20,5))
    plt.subplot(1,4,1); plt.imshow( Etot, vmin=vmin,vmax=-vmin ) ;plt.colorbar()  ;plt.title('Etot')
    #plt.subplot(1,4,2); plt.imshow( Ek  , vmin=vmin,vmax=-vmin ) ;plt.colorbar()  ;plt.title('Ek'  )
    #plt.subplot(1,4,3); plt.imshow( Eaa , vmin=vmin,vmax=-vmin ) ;plt.colorbar()  ;plt.title('Eaa' )
    #plt.subplot(1,4,4); plt.imshow( Eae , vmin=vmin,vmax=-vmin ) ;plt.colorbar()  ;plt.title('Eel' )
    plt.subplot(1,4,2); plt.imshow( Ek  ) ;plt.colorbar()  ;plt.title('Ek'  )
    plt.subplot(1,4,3); plt.imshow( Eaa ) ;plt.colorbar()  ;plt.title('Eaa' )
    plt.subplot(1,4,4); plt.imshow( Eae ) ;plt.colorbar()  ;plt.title('Eel' )

    #plt.legend()
    plt.grid()
    plt.show()





