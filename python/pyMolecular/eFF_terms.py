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

const_bohr_eVA =  0.52917724

const_K     =  const_hbar**2/(2*const_Me)    # see schroedinger equation : https://en.wikipedia.org/wiki/Schr%C3%B6dinger_equation#Preliminaries
const_El    =  const_e**2/(4.*np.pi*const_eps0)

const_Hartree = const_El**2/const_K
const_Ry      = 0.5 * const_Hartree
const_Ry_eV   = 13.6056925944
const_El_eVA  = const_El/( const_e*const_Angstroem )

# prefactor ( h^2/(2me)) from Schroedinger equation can be derived from Hartree energy (https://en.wikipedia.org/wiki/Hartree)
# Eh= 2Ry = me * (e^2/(4pi eps0 hbar ))^2 = (me/(hbar^2)) * (e^2/(4pi eps0))^2 = (2*me/(hbar^2))/2 * const_El^2 = (const_El^2) / (2*const_K)
# 2Eh=4Ry =  const_El^2/ const_K
# const_K = const_El^2 / ( 4*Ry ) 

#const_K_eVA  = (const_El_eVA**2)/(2*const_Ry_eV)  # is wrong
const_K_eVA  = (const_El_eVA**2)/(4*const_Ry_eV)     # should be, see 

print "const_Ry const_Ry_eV ",   const_Ry, const_Ry/const_eV
print "const_El, const_El_eVA ", const_El, const_El_eVA,  (const_El/const_El_eVA )/const_Angstroem
print "const_K, const_K_eVA   ", const_K,  const_K_eVA,   (const_K/const_K_eVA   )/(const_Angstroem**2)


const_K_eVA_  = const_K  / ( const_eV * const_Angstroem**2 )
const_El_eVA_ = const_El / ( const_eV * const_Angstroem    )

print "const_El_eVA const_El_eVA_ ",  const_El_eVA, const_El_eVA_ 
print "const_K_eVA  const_K_eVA_   ", const_K_eVA, const_K_eVA_ 



#const_K  =   const_hbar**2/const_Me   #   [ eV * A^2 ]
#const_K  = 0.1* 30.0824137226  # [eV*A^2] hbar[J.s]^2/(Me [kg])   /  (  eV[J]*A^2[m])    # (6.62607015e-34^2/9.10938356e-31)/1.602176620898e-19/10e-20
#const_Ke =  1.5*const_K

#const_El =  14. # 14 (1./((4*np.pi*const_eps0))

sqrt2 = np.sqrt(2.)

def Kinetic( s ):
    '''
    Ek = (hbar^2/me) (3./2.) 1/s^2      ... is wrong
    Ek = (hbar^2/(2*me)) (3./2.) 1/s^2  ... should be, see https://en.wikipedia.org/wiki/Schr%C3%B6dinger_equation#Definition
    '''
    return const_K_eVA*(3/2.)/(s**2)

def El( r, qq, si=0, sj=0 ):
    if si>0:
        if sj>0:
            si = np.sqrt( si**2 + sj**2 )
        return const_El_eVA * (qq/r) * spc.erf( r/s )
    else:
        return const_El_eVA * (qq/r)

def El_aa( r, qq ):
    return const_El_eVA * (qq/r)

def El_ae( r, qq, s ):
    return const_El_eVA * (qq/r) * spc.erf( r/s )

def El_ee( r, qq, si, sj ):
    s = np.sqrt( si**2 + sj**2 )
    return const_El_eVA * (qq/r) * spc.erf( r/s )

def getT( r, si, sj ):
    #print "getT si, sj ", si, sj
    #   r = r * 1.125
    #   s = s*0.9
    si2 = si**2
    sj2 = sj**2
    r2  = r**2
    #return const_K * ( 1.5*( (si2+sj2)/(si2*sj2) )   - 2.*( 3.*(si2+sj2) - 2.*r2 )/( si2 + sj2 )**2 )
    return const_K_eVA * ( 1.5*( 1./si2 + 1./sj2 )   - 2.*( 3.*(si2+sj2) - 2.*r2 )/( si2 + sj2 )**2 )

def getAmp( si, sj ):
    si2 = si**2
    sj2 = sj**2
    #return const_K_eVA * ( 1.5*( 1./si2 + 1./sj2 )   - 2.*( 3.*(si2+sj2) - 2.*0 )/( si2 + sj2 )**2 )
    #return const_K_eVA * ( 1.5*( 1./si2 + 1./sj2 )   - 2.*( 1.*(si2+sj2) )/( si2 + sj2 )**2 )
    #return const_K_eVA * 2.2*( 1.5*( 1/si2 + 1/sj2 ) - 4.9/( si2 + sj2 ) )
    #return const_K_eVA * 2.2*( 1.5*( (si2 + sj2)/(si2*sj2) ) - 4.9/( si2 + sj2 ) )
    #return const_K_eVA * 2.2*( 1.5*(si2*si2 + sj2*sj2) - 1.9*(si2*sj2)  )/((si2*sj2)*(si2+sj2))
    #return const_K_eVA * 2.2*( 1.5*(si2*si2 + sj2*sj2) - 1.9*(si2*sj2)  )/((si2*sj2)*(si2+sj2))
    return const_K_eVA * 3.3*( si2*si2 + sj2*sj2 - 1.25*(si2*sj2) )/((si2*sj2)*(si2+sj2))
    #return const_K_eVA * 3.14*( si2*si2 + sj2*sj2 - 1.25*(si2*sj2) )/((si2*sj2)*(si2+sj2))
    #return const_K_eVA * ( 1.5*( 1./si2 + 1./sj2 ) ) 
    #return const_K_eVA * ( 1.5*( 1./si2 + 1./sj2 )  - 2.*3./( si2 + sj2 ) )

def getS( r, si, sj ):
    #print "getS si, sj ", si, sj
    #   r = r * 1.125
    #   s = s*0.9
    si2 = si**2
    sj2 = sj**2
    r2  = r**2
    return ( 2.*(si*sj)/(si2+sj2) )**1.5 * np.exp(  -r2/( si2 + sj2 ) )

'''
def EPauli( r, si, sj, rho=0.2 ):
    T = getT( r, si, sj )
    S = getS( r, si, sj )
    S2 = S**2
    # ( S2*(1+S2) + (1-rho)* S2*(1-S2) )  / (1-S2*S2 )
    # ( S2+S2*S2 + (1-rho)*(S2-S2*S2) )  / (1-S2*S2 )
    # ( ( (2-rho)*S2 +rho*S2*S2 )  / (1-S2*S2 )
    return T * ( (S2/(1.-S2))   + ( 1.-rho )*(S2/(1.+S2))     )

def EPauli_pair( r, si, sj, rho=0.2 ):
    T  = getT( r, si, sj )
    S  = getS( r, si, sj )
    S2 = S**2
    return T * ( rho*S2/(1.+S2) )
'''

def EPauli( r, si, sj, anti=False, rho=0.2, kr=1.125, ks=0.9 ):
    r  = r*kr
    si = si*ks
    sj = sj*ks
    T = getT( r, si, sj )
    S = getS( r, si, sj )
    S2 = S**2
    if anti:
        return T * ( rho*S2/(1.+S2) )
    else:
        return T * ( (S2/(1.-S2))   + ( 1.-rho )*(S2/(1.+S2)) )

def DensOverlap( r, si, sj, amp=10 ):
    s2 = si**2+sj**2
    #amp *= 1.4/s2
    #amp *= 0.7/(si*sj)
    #amp *= (1/si**2 + 1/sj**2)
    #amp  *= (si**2+sj**2)/(si*sj)**2
    #amp  *= (si+sj)**2/(si*sj)**2
    #amp  *= (1+(si-sj)**2)/min(si,sj)**2
    #amp  *= 0.5*(1+4*(si-sj)**2) *( 1/si**2 + 1/sj**2 )
    a  = 2*si*sj/s2
    e1 = amp*a**3
    e2 = np.exp( -2*(r**2/s2) )
    return e1*e2

def Hatom( s ):
    Ek  = Kinetic( s )
    Eae = El_ae( 0.01, -1., s )
    #Etot = Ek+Eae
    return Ek,Eae

def Hatom_au_( s ):
    # eq.10 from http://aip.scitation.org/doi/10.1063/1.3272671
    Ek  =  1.5 * (1/s**2)
    #Eae = 1.35*   -np.sqrt(8./np.pi) * (1/s)
    Eae = -np.sqrt(8./np.pi) * (1/s)
    return Ek,Eae

def H2cation( rHH, s, cr=0.5 ):
    Ek   = Kinetic( s )                  # kinetic energy of electron
    Eaa  = El_aa( rHH,  1. )             # Coulomb repulsion  nuclei_1 + nuclei_2
    Eae  = El_ae( rHH*(cr   ), -1., s )  # Coulomb attraction electron + nuclei_1
    Eae += El_ae( rHH*(1.-cr), -1., s )  # Coulomb attraction electron + nuclei_2
    return Ek, Eae, Eaa

def H2molecule( r, s, cr=0.5 ):
    Ek    = 2*Kinetic( s )                      # kinetic energy of electron_1 and electron_2
    Eaa   =   El_aa( r,  +1. )                  # Coulomb repulsion   nuclei_1 * nuclei_2
    Eae   = 2*El_ae( r*(cr   ), -1., s )        # Coulomb attraction (electron_1 * nuclei_1)   +   (electron_2 * nuclei_2)
    Eae  += 2*El_ae( r*(1.-cr), -1., s )        # Coulomb attraction (electron_1 * nuclei_2)   +   (electron_2 * nuclei_1)
    Eee   =   El_ee( r*(1.-2.*cr), +1., s, s )  # Coulomb repulsion   electron_1 * electron_2
    EPaul =   EPauli( r*(1.-2.*cr), s, s, anti=True )   # Pauli repulsion electron_1 * electron_2
    return Ek, Eae, Eaa, Eee, EPaul

# =========== Macro Functions

def run_ee_onsite():
    # ============= e-e onsite 
    r0 = 0.01
    ss = np.arange( 0.25, 5.0, 0.1 )
    rho=0.2; kr=1.125; ks=0.9
    r_ = r0*kr
    s_ = ss*ks
    T = getT( r_, s_, s_ )
    S = getS( r_, s_, s_ )
    S2 = S**2
    EPminus  =  T * ( rho*S2/(1.+S2) )
    EPplus   =  T * ( (S2/(1.-S2))  + ( 1.-rho )*(S2/(1.+S2)) )

    plt.figure()
    plt.title( 'Onsite (R= %g [A])' %r0 )
    plt.xlabel('sigma[A]')
    plt.plot( ss, S,   ':', label="S" )
    plt.plot( ss, T,   ':', label="T" )
    plt.plot( ss, EPplus, 'b', label="EP+" )
    plt.plot( ss, EPminus,'c', label="EP-" )
    plt.legend()
    plt.grid()
    #plt.show(); exit()

def run_ee():
    rs = np.arange( 0.1, 6.0, 0.05 )
    #ss = [0.5, 1.0, 1.5 ]
    ss = [0.25, 1.0, 2.5 ]

    rho=0.2; kr=1.125; ks=0.9

    plt.figure(figsize=(13,10))
    for i,si in enumerate(ss):
        for j,sj in enumerate(ss):
            Eee = El_ee( rs, +1., si, sj )
            r_ = rs*kr
            #s_ = s*ks
            T = getT( r_, si*ks, sj*ks )
            S = getS( r_, si*ks, sj*ks )
            S2 = S**2
            EPminus  =  T * ( rho*S2/(1.+S2) )
            EPplus   =  T * ( (S2/(1.-S2))  + ( 1.-rho )*(S2/(1.+S2)) )

            #amp  = 10*(1+(si-sj)**2)/min(si,sj)**2
            #amp  = 10/min(si,sj)**2
            #amp  = 10*(1+0.6*abs(si-sj))/min(si,sj)**2
            #amp  = 10*(si/sj+sj/si)
            #amp  = 10
            #amp  = T*1.8
            amp   = getAmp( si, sj )

            EPdens   =   DensOverlap( rs, si, sj, amp=amp )
            plt.subplot(3,3,3*j+i+1)

            #plt.plot( rs, S,   ':', label="S" )
            #plt.plot( xs, T,   ':', label="T" )
            #plt.plot( rs, Eee ,   'r', label="Eee" )
            plt.plot( rs, EPplus, 'b', label="EP+" )
            #plt.plot( rs, EPminus,'c', label="EP-" )
            plt.plot( rs, EPdens,  'm', label="EPdens"  )
            plt.title( 'sigma (%g,%g)' %(si,sj) )
            plt.legend()
            plt.grid()
            #plt.plot( ys, Etot, 'k', label="Etot" )
    #plt.show(); exit()

def run_Hatom():
    #Ek  = Kinetic( ys )
    #Eae = El_ae( 0.01, -1., ys )

    Ek,Eae = Hatom( ys )
    Etot = Ek+Eae
    plt.figure()
    plt.plot( ys, Ek ,  'r', label="Ek" )
    plt.plot( ys, Eae,  'b', label="Eae" )
    plt.plot( ys, Etot, 'k', label="Etot" )

    imin = np.argmin( Etot )
    print "H-atom Rmin Emin(Ek,Eel) ", ys[imin], Etot[imin], Ek[imin], Eae[imin] 
    print "should be Rmin=%g[A](1.88[bohr]) Emin=%g[eV](-0.424[hartree])" %(1.88*0.52917724, -0.424*const_Hartree/const_eV )

    EHatom = Etot[imin]

    plt.title("H-atom")
    plt.xlabel("s [A]")
    plt.legend()
    plt.grid()
    #plt.show(); exit()

def run_Hcation():
    Xs,Ys = np.meshgrid( xs,ys )

    Ek, Eae, Eaa = H2cation( Xs, Ys, cr=0.5 )

    Etot = Ek + Eaa + Eae

    #Emin = Etot.min()
    imin = np.unravel_index( np.argmin(Etot), Etot.shape )
    Emin = Etot[imin]
    Rmin = xs[imin[0]]
    Smin = ys[imin[1]]
    print "H2cation Rmin, Smin Emin Ebond ", Rmin, Smin, Emin, Emin-EHatom
    vmin=-20.0 # [eV]
    vmax=-vmin

    plt.figure(figsize=(20,5))
    plt.subplot(1,4,1); plt.imshow( Etot, origin='image', extent=extent, vmin=vmin,vmax=vmax ) ;plt.title('Etot')
    plt.subplot(1,4,2); plt.imshow( Ek  , origin='image', extent=extent, vmin=vmin,vmax=vmax ) ;plt.title('Ek'  )
    plt.subplot(1,4,3); plt.imshow( Eaa , origin='image', extent=extent, vmin=vmin,vmax=vmax ) ;plt.title('Eaa' )
    plt.subplot(1,4,4); plt.imshow( Eae , origin='image', extent=extent, vmin=vmin,vmax=vmax ) ;plt.title('Eel' )
    #plt.subplot(1,4,2); plt.imshow( Ek  , origin='image', extent=extent ) ;plt.colorbar()  ;plt.title('Ek'  )
    #plt.subplot(1,4,3); plt.imshow( Eaa , origin='image', extent=extent ) ;plt.colorbar()  ;plt.title('Eaa' )
    #plt.subplot(1,4,4); plt.imshow( Eae , origin='image', extent=extent ) ;plt.colorbar()  ;plt.title('Eel' )

def run_Hmolecule():
    Ek, Eae, Eaa, Eee, EPaul = H2molecule( Xs, Ys, cr=0.49 )
    Etot = Ek + Eae + Eaa + Eee + EPaul
    #Emin = Etot.min()
    imin = np.unravel_index( np.argmin(Etot), Etot.shape )
    Emin = Etot[imin]
    Rmin = xs[imin[0]]
    Smin = ys[imin[1]]
    print "H2molecule Rmin, Smin Emin Ebond ", Rmin, Smin, Emin, Emin - 2*EHatom
    vmin=-50.0 # [eV]
    vmax= 0.0 # [eV]

    plt.figure( figsize=(18,3) )
    plt.subplot(1,6,1); plt.imshow( Etot, origin='image', extent=extent, vmin=vmin,vmax=vmax ) ;plt.colorbar()  ;plt.title('Etot')
    #plt.subplot(1,6,2); plt.imshow( Ek  , origin='image', extent=extent, vmin=vmin,vmax=vmax ) ;plt.colorbar()  ;plt.title('Ek'  )
    #plt.subplot(1,6,3); plt.imshow( Eaa , origin='image', extent=extent, vmin=vmin,vmax=vmax ) ;plt.colorbar()  ;plt.title('Eaa' )
    #plt.subplot(1,6,4); plt.imshow( Eae , origin='image', extent=extent, vmin=vmin,vmax=vmax ) ;plt.colorbar()  ;plt.title('Eea' )
    #plt.subplot(1,6,5); plt.imshow( Eee , origin='image', extent=extent, vmin=vmin,vmax=vmax ) ;plt.colorbar()  ;plt.title('Eee' )
    #plt.subplot(1,6,6); plt.imshow( EPaul, origin='image', extent=extent, vmin=vmin,vmax=vmax ) ;plt.colorbar()  ;plt.title('EPaul')
    plt.subplot(1,6,2); plt.imshow( Ek  , origin='image', extent=extent  ) ;plt.colorbar()  ;plt.title('Ek'  )
    plt.subplot(1,6,3); plt.imshow( Eaa , origin='image', extent=extent ) ;plt.colorbar()  ;plt.title('Eaa' )
    plt.subplot(1,6,4); plt.imshow( Eae , origin='image', extent=extent ) ;plt.colorbar()  ;plt.title('Eea' )
    plt.subplot(1,6,5); plt.imshow( Eee , origin='image', extent=extent  ) ;plt.colorbar()  ;plt.title('Eee' )
    plt.subplot(1,6,6); plt.imshow( EPaul, origin='image', extent=extent ) ;plt.colorbar()  ;plt.title('EPaul')

if __name__ == "__main__":

    extent=( 0.25,6.0,  0.25,4.0 )
    xs = np.arange( extent[0], extent[1], 0.05 )
    ys = np.arange( extent[2], extent[3], 0.1  )

    #run_ee_onsite()
    #run_ee()
    run_Hatom()
    #run_Hcation()
    #run_Hmolecule()
    plt.show()





