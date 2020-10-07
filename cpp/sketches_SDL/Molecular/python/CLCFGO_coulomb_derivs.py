#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

import scipy.special as spc


'''

We have oribtals 
Phi_1 = c_a * chi_a + c_b * chi_b
Phi_2 = c_c * chi_c + c_d * chi_d 

where chi_i are gaussian type basis functions
and   c_i   are expansion coefficients

electron density of molecular orbitals Rho_1 = <phi_1|phi_1>
can be expressed using auxuliary gaussian basisfunctionc 
rho_ab = chi_a * chi_b
Rho_1  = sim_ab c_a*c_b*S_ab * rho_ab   
       = sim_ab q_ab * rho_ab

where q_ab = c_a*c_b*S_ab is charge of the auxuliary electron blob
with  S_ab being overlap integral between the basis functions chi_a*chi_b
we can use collective index i=ab and j=cd 

qi = Sab*ca*cb
qj = Scd*cc*cd

The repulsion between blobs qi,qj can be expressed as
qi*qj /  = 

'''


const_hbar_SI      = 1.054571817e-34;    #< [J.s]  #6.582119569e-16 # [eV/s]
const_Me_SI        = 9.10938356e-31;     #< [kg]
const_e_SI         = 1.602176620898e-19; #< [Coulomb]
const_eps0_SI      = 8.854187812813e-12; #< [F.m = Coulomb/(Volt*m)]
const_eV_SI        = 1.602176620898e-19; #< [J]
const_Angstroem_SI = 1.0e-10;

const_K_SI     =  const_hbar_SI*const_hbar_SI/const_Me_SI;
const_El_SI    =  const_e_SI*const_e_SI/(4.*np.pi*const_eps0_SI);
const_Ry_SI    = 0.5 * const_El_SI*const_El_SI/const_K_SI;

const_Ry_eV  = 13.6056925944;
const_El_eVA = const_El_SI/( const_e_SI*const_Angstroem_SI );
const_K_eVA  = (const_El_eVA*const_El_eVA)/(2*const_Ry_eV);
const_Ke_eVA = const_K_eVA*1.5;


def Coulomb( r, s ):
    # ToDo: maybe we can do without s=sqrt(s2) and r=sqrt(r2)
    #constexpr const double const_F2 = -2.*sqrt(2./np.pi);
    #const_F2 = M_2_SQRTPI * M_SQRT2;
    M_SQRT2 = 1.41421356237 
    const_F2 = 2*np.sqrt(2/np.pi)
    ir   = 1./r                         #(r+1.e-8);
    is_  = 1./s                         #(s+1.e-8);
    r_s  = r*is_
    r_2s = M_SQRT2 * r_s
    e1   = ir * const_El_eVA
    e2   = spc.erf(  r_2s      )
    g    = np.exp( -r_2s*r_2s ) * const_F2
    f1   = -e1*ir
    f2   = g*is_
    e1f2 = e1*f2
    fr = (f1*e2 + e1f2)*ir
    fs =          e1f2 *r_s * is_
    E  = e1 * e2
    return E,fr,fs

def product3D_s_deriv( si,pi, sj,pj ):
    ''' returns
    S,   p,
    dSsi, dSsj,
    dXsi, dXsj,
    dXxi, dXxj,
    dCsi, dCsj, dCr
    '''
    si2   = si*si
    sj2   = sj*sj
    s2    = si2 + sj2
    is2   = 1/s2
    is4   = is2*is2
    sqrtis2 = np.sqrt(is2)

    s      =  si*sj*sqrtis2               # size
    p      =  pj*(si2*is2) + pi*(sj2*is2) # position
    #X      =  ( si2*xj + sj2*xi )*inv;

    inv3_2 = sqrtis2*is2
    dSsi   = sj*sj2*inv3_2
    dSsj   = si*si2*inv3_2

    dp     = pi-pj
    dXsi   = dp*(-2*si*sj2*is4)
    dXsj   = dp*( 2*sj*si2*is4)

    dXxi   = sj2*is2
    dXxj   = si2*is2

    #r2 = dp.norm2()
    r2  = dp*dp

    a2   = 2.*(si*sj)*is2
    a    = np.sqrt(a2)
    e1   = a2*a
    e2   = np.exp( -r2*is2 )

    f1   = 3.*a  * (si2-sj2)*is4
    f2   = 2.*e2 * r2*is4

    dCsi = e1*f2*si - e2*f1*sj
    dCsj = e1*f2*sj + e2*f1*si
    dCr  = e1*e2*(-2.*is2)          # derivative is correct, tested !

    #double logC =  wxi*xi + wxj*xj - wx*X;
    #double C   = np.exp(-logC) * Ci * Cj

    S = e1*e2 # Overlap
    return S,s,p, dCr*dp, (dSsi,dXsi,dXxi,dCsi), (dSsj,dXsj,dXxj,dCsj)











def checkNumDeriv( x, func, dfunc, name ):
    dy         = dfunc( x )
    y          = func(x)
    dynum,xnum = numDeriv( x, y )
    #print y
    #print "y.shape, ynum.shape ", y.shape, ynum.shape
    plotVsNum( x,dy,dynum, name )
    plt.plot(x, y,'-.', label=name+"_F" )




if __name__ == "__main__":

    #s  = np.arange( 0.1, 5.0, 0.05 )
    #rs = np.arange( 0.1, 5.0, 0.05 )
    #S  = np.arange( 1.25, 5.0, 0.05 )
    #r  = 1.5 + 0.*s

    ca = 1.0
    cb = 1.0
    cc = 1.0
    cd = 1.0

    sa = 1.0
    sb = 1.0
    sc = 1.0
    sd = 1.0

    dx =  0.05
    xa =  np.arange( 0.01, 5.0, dx )
    xb =  0.0
    xc =  -0.5
    xd =   0.0


    # overlaps
    Sab, si, xab, dQab, dA, dB = product3D_s_deriv( sa,xa, sb,xb )
    Scd, sj, xcd, dQcd, dC, dD = product3D_s_deriv( sc,xc, sd,xd )
    # coulomb
    s2   = si*si + sj*sj;
    s    = np.sqrt(s2);
    E, fx, fs = Coulomb( xab-xcd, s )

    #plt.plot( xa, Sab , label='Sab' )
    #plt.plot( xa, dQab, label='dSab_ana' )
    #plt.plot( (xa[1:]+xa[:-1])*0.5, (Sab[1:]-Sab[:-1])/dx,':', label='dSab_num' )

    # Q: Why we dont need derivatives of charge ????
    #Fx = -fx*0.5*dA[1] # This works for zero initial distance between blobs
    Fx = -fx*0.5*dA[1]

    plt.plot( xa, E , label='E' )
    plt.plot( xa, Fx, label='dEdx_ana' )
    plt.plot( (xa[1:]+xa[:-1])*0.5, (E[1:]-E[:-1])/dx,':', label='dEdx_num' )


    plt.grid()
    plt.legend()
    plt.show()

