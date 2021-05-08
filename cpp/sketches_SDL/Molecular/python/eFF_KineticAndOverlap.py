#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

import scipy.special as spc

'''

To evaluate Change of Kinetic energy due to orthogonalization (Deltat E_kin ) which is Valence-Bond model for Pauli repulsion 
see:
[1] eq.4 in http://aip.scitation.org/doi/10.1063/1.3272671
            Jaramillo-Botero, A., Su, J., Qi, A. & Goddard, W. A. Large-scale, long-term nonadiabatic electron molecular dynamics for describing material properties and phenomena in extreme environments. J. Comput. Chem. 32, 497-512 (2011).
[2] eq.2 in http://doi.wiley.com/10.1002/jcc.21637
            Su, J. T. & Goddard, W. A. The dynamics of highly excited electronic systems: Applications of the electron force field. J. Chem. Phys. 131, 244501 (2009).

Derivation in  Gauss_KineticAndOverlap.wxmx
/SimpleSimulationEngine/cpp/sketches_SDL/Molecular/notes/Gauss_KineticAndOverlap.wxmx

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

# ======= Functions

# ToDo : Kinetic and Overlap share much of calculations => make sense to calculate them together in one function

def overlap_(r,si,sj):
    # NOTE : this gaussian is not normalized !!!
    const = (2*np.pi)**(1.5)
    s2    = si**2 + sj**2
    r2    = r**2
    g     = np.exp( -r**2/(2*s2) )
    E     = const * (si*sj)**3/( s2**1.5 ) * g
    return E

def kinetic_(r,si,sj):
    # NOTE : this gaussian is not normalized !!!
    const = (2*np.pi)**(1.5)
    s2    = si**2 + sj**2
    r2    = r**2
    g     = np.exp( -r**2/(2*s2) )
    tau   = -(r2 - 3 * s2)/(s2**2)
    E     = const * (si*sj)**3/( s2**1.5 ) * g
    return E * tau


def overlap(r,si,sj):
    # NOTE : this gaussian IS normalized !!!
    #const = (2*np.pi)**(1.5) / (np.pi)**(1.5)
    const = 2**1.5
    s2    = si**2 + sj**2
    r2    = r**2
    norm  = 1/(si*sj)**(1.5)
    g     = np.exp( -r**2/(2*s2) )
    E     = const * norm * (si*sj)**3/( s2**1.5 ) * g 
    return E

def kinetic(r,si,sj):
    const = 2**1.5
    s2    = si**2 + sj**2
    r2    = r**2
    norm  = 1/(si*sj)**(1.5)
    g     = np.exp( -r**2/(2*s2) )
    tau   = -(r2 - 3 * s2)/(s2**2)
    E     = const * norm * (si*sj)**3/( s2**1.5 ) * g
    return E * tau

def tau(r,si,sj): # Kinetic/Overlap Tij/Sij
    s2    = si**2 + sj**2
    r2    = r**2
    tau   = -(r2 - 3 * s2)/(s2**2)
    return tau




if __name__ == "__main__":
    rs = np.arange( 0.0, 6.0, 0.05 )
    #ss = [0.25, 1.0, 2.5 ]
    #si=1.0; sj=1.0;
    #si=0.5; sj=0.5;
    #si=2.0; sj=2.0;
    si=0.5**0.5; sj=0.5**0.5;
    #si=0.5; sj=2.0;
    #si=2.0; sj=0.5;
    
    Ss   = overlap( rs, si, sj )   ;print( Ss[0] )
    Ts   = kinetic( rs, si, sj )   ;print( Ts[0] )
    taus = tau    ( rs, si, sj )   ;print( taus[0] )
    
    plt.plot( rs, Ss,   label="Overlap S12" )
    plt.plot( rs, Ts,   label="Kinetic T12" )
    plt.plot( rs, taus, label="tau T12/S12" )
    
    plt.axhline(0,c='k',ls='--',lw=2 )
    plt.ylim(-2,6)
    plt.grid()
    plt.legend()
    plt.show()


