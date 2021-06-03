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

M_SQRT2 = 1.41421356237
M_SQRT1_2 = 1/M_SQRT2 

def Coulomb( r, s ):
    '''
    double ir   = 1./r; //(r+1.e-8);
    double is   = 1./s; //(s+1.e-8);
    double r_s  = r*is;
    double r_2s = M_SQRT1_2 * r_s; // This is for charge-density blobs (assuming si,sj comes from charge denisty)
    //double r_2s = r_s;
    //double r_2s = M_SQRT2   * r_s; // This is for wavefunction blobs (assuming si,sj comes from wavefunction)
    double e1   = ir * const_El_eVA;
    double e2   = erf(  r_2s      );            // ToDo : this should be possible to compute together !!!
    double g    = exp( -r_2s*r_2s ) * const_F2;
    double f1   = -e1*ir;
    double f2   = g*is*0.5;
    double e1f2 = e1*f2;
    fr          = (f1*e2 + e1f2)*ir;
    fs          =          e1f2 *r_s * is;
    return e1 * e2;
    '''
    # ToDo: maybe we can do without s=sqrt(s2) and r=sqrt(r2)
    #constexpr const double const_F2 = -2.*sqrt(2./np.pi);
    #const_F2 = M_2_SQRTPI * M_SQRT2;
    const_F2 = 2*np.sqrt(2/np.pi)
    ir   = 1./r                         #(r+1.e-8);
    is_  = 1./s                         #(s+1.e-8);
    r_s  = r*is_
    r_2s = M_SQRT1_2 * r_s; # This is for charge-density blobs (assuming si,sj comes from charge denisty)
    #r_2s = r_s;
    #r_2s = M_SQRT2   * r_s; # This is for wavefunction blobs (assuming si,sj comes from wavefunction)
    e1   = ir * const_El_eVA
    e2   = spc.erf(  r_2s      )
    g    = np.exp( -r_2s*r_2s ) * const_F2
    f1   = -e1*ir
    #f2   = g*is_        # This is for wavefunction blobs (assuming si,sj comes from wavefunction)
    f2   = g*is_*0.5     # This is for charge-density blobs (assuming si,sj comes from charge denisty)
    e1f2 = e1*f2
    fr = (f1*e2 + e1f2)*ir
    fs =          e1f2 *r_s * is_
    E  = e1 * e2

    #for i in range(len(r)):
    #    print "Gauss::Coulomb r %g s %g E %g fr %g " %(r[i],s, E[i], fr[i]  )


    return E,fr,fs

def product3D_s_deriv( si,pi, sj,pj ):
    ''' returns
    S,   p,
    dSsi, dSsj,
    dXsi, dXsj,
    dXxi, dXxj,
    dCsi, dCsj, dCr
    '''

    '''
    #define _Gauss_sij_aux( si, sj ) \
    double si2  = si*si; \
    double sj2  = sj*sj; \
    double s2   = si2 + sj2; \
    double is2  = 1/s2; \
    double is4  = is2*is2; \
    '''
    si2   = si*si
    sj2   = sj*sj
    s2    = si2 + sj2
    is2   = 1/s2
    is4   = is2*is2

    ''''
    #define _Gauss_product(pi,pj,si,sj)    \
    double sqrt_is2 = sqrt(is2);       \
    double size_ij  = si*sj*sqrt_is2;            \
    Vec3d  pos_ij   = pj*(si2*is2) + pi*(sj2*is2); \
    '''
    sqrt_is2 = np.sqrt(is2)
    s        = si*sj*sqrt_is2               # size
    p        = pj*(si2*is2) + pi*(sj2*is2) # position

    dp     = pi-pj
    r2  = dp*dp  # r2 = dp.norm2()    in 1D

    '''
    #define _Gauss_product_derivs(dSsi,dSsj,dXsi,dXsj,dXxi,dXxj){ \
    double inv3_2 = is2 * sqrt_is2; \
    dSsi   = sj*sj2*inv3_2; \
    dSsj   = si*si2*inv3_2; \
    dXsi   = dp*(-2*si*sj2*is4); \
    dXsj   = dp*( 2*sj*si2*is4); \
    dXxi   = sj2*is2; \
    dXxj   = si2*is2; \
    '''
    inv3_2 = is2*sqrt_is2
    dSsi   = sj*sj2*inv3_2
    dSsj   = si*si2*inv3_2
    dXsi   = dp*(-2*si*sj2*is4)
    dXsj   = dp*( 2*sj*si2*is4)
    dXxi   = sj2*is2
    dXxj   = si2*is2

    '''
    #define _Gauss_overlap( r2, si, sj ) \
    double sisj    = si*sj; \
    double inv_sisj= 1/sisj; \
    double g       = exp( -r2/(2*s2) );    \
    double S       =  (2*M_SQRT2) * g *  sisj*sisj*is2 * sqrt( inv_sisj*is2 ) ;  \
    double dS_dr   = -S * is2; \
    double inv_si  = sj*inv_sisj; \
    double inv_sj  = si*inv_sisj; \
    double S_s4    = S * is4; \
    double dS_dsi  = S_s4 * ( si2*r2 + 3*sj2*s2 ) * inv_si; \
    double dS_dsj  = S_s4 * ( sj2*r2 + 3*si2*s2 ) * inv_sj; \
    dS_dsi        -= 1.5*S * inv_si;    \
    dS_dsj        -= 1.5*S * inv_sj;    \
    '''

    sisj    = si*sj
    inv_sisj= 1/sisj
    g       = np.exp( -r2/(2*s2) )
    S       = (2*M_SQRT2) * g *  sisj*sisj*is2 * np.sqrt( inv_sisj*is2 )
    dS_dr   = -S * is2
    inv_si  = sj*inv_sisj
    inv_sj  = si*inv_sisj
    S_s4    = S * is4
    dS_dsi  = S_s4 * ( si2*r2 + 3*sj2*s2 ) * inv_si
    dS_dsj  = S_s4 * ( sj2*r2 + 3*si2*s2 ) * inv_sj
    dS_dsi -= 1.5*S * inv_si
    dS_dsj -= 1.5*S * inv_sj

    '''
    # ==== Overlaps OLD
    a2   = 2.*(si*sj)*is2
    a    = np.sqrt(a2)
    e1   = a2*a
    e2   = np.exp( -r2*is2 )
    f1   = 3.*a  * (si2-sj2)*is4
    f2   = 2.*e2 * r2*is4
    dCsi = e1*f2*si - e2*f1*sj
    dCsj = e1*f2*sj + e2*f1*si
    C    = e1*e2        # Overlap
    dCr  = C*(-2.*is2)  # derivative is correct, tested !
    '''

    
    return S,s,p, dS_dr*dp, (dSsi,dXsi,dXxi,dS_dsi), (dSsj,dXsj,dXxj,dS_dsj)











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

    #ks = np.sqrt(0.5)
    sa = 0.2
    sb = 0.2
    sc = 0.2
    sd = 0.2

    #n  = 40; x0 = 0.00001
    n  = 10; x0 = 0.75
    #n  = 1;  x0 = 1.0
    
    dx =  0.05
    xa =  np.arange( x0, x0+dx*n, dx )
    xb =  1.0
    xc =  1.0
    xd =  2.0

    xs_ = (xa[1:]+xa[:-1])*0.5

    # overlaps
    Sab, si, xab, dSab, dA, dB = product3D_s_deriv( sa,xa, sb,xb )
    Scd, sj, xcd, dScd, dC, dD = product3D_s_deriv( sc,xc, sd,xd )
    # coulomb
    s2        = si*si + sj*sj
    s         = np.sqrt(s2)
    r         = xab-xcd
    e, fx, fs = Coulomb( r, s )
    dXxi = dA[2] + xa*0 

    #plt.plot( xa, Sab , label='Sab' )
    #plt.plot( xa, r   , label='r'   )
    #plt.plot( xa, dQab, label='dSab_ana' )
    #plt.plot( xs_, (Sab[1:]-Sab[:-1])/dx,':', label='dSab_num' )

    qi = 2*Sab
    qj = 2*Scd

    qij  = qi*qj
    dQij = qj* 2*dSab

    #Vec3d dSdp = Rij*(dSr*cij);
    #Vec3d Fq   = dSdp*dEdQ;
    dSdp = dSab*2;
    dEdQ = e*qj;

    # Q: Why we dont need derivatives of charge ????
    #Fx = -fx*0.5*dA[1] # This works for zero initial distance between blobs
    #Fx  = fx*r*dXxi
    Fpi = fx*r*qij # see 
     
    #print "Scd, 4*Scd ", Scd, 4*Scd
    #print "For some reason each charge is scaled by 2.0"

    E = e*qij
    #F = Fpi*dXxi + e*dQij   # total derivtive F = dE/dx = d(e*qi)/dx
    F = Fpi*dXxi + dEdQ*dSdp


    # Note:   e,fx=de/dx   are NOT multiplied by charge Qij
    #         Total force Fx = dE/dx = d(e*q)/dx = q*(de/dx) + e*(dq/dx)

    for i in range(len(r)):
        print "fromRho x %g  F %g =(Fpi %g * dXxi %g)+(dSdp %g * dEdQ %g)"    %(   xa[i], F[i],  Fpi[i], dXxi[i], dSdp[i], dEdQ[i] ) ;
        #print "fromRho x %g  dSdp %g =(dSab %g * cab %g)" %(  xa[i], dSdp[i], dSab[i], 1 ) ;
        #print "fromRho x %g  dEdQ %g =(e %g * qj %g)" %(  xa[i], dEdQ[i], e[i], qj ) ;
        pass

    # ==== Derivative of Coulomb term without considering changes of Charges 
    #plt.plot( xa, e ,  label='e' )
    #plt.plot( xa, Fx,  label='dedx_ana' )
    #plt.plot( xs_, (e[1:]-e[:-1])/dx,':', label='dedx_num' )
    
    # ==== Derivative of Coulomb term with considering the Charges
    plt.plot( xa, E,  label='E' )
    plt.plot( xa, F,  label='dEdx_ana' )
    plt.plot( xs_, (E[1:]-E[:-1])/dx,':', label='dEdx_num', lw=3 )
    #plt.plot( xa, fxi, label='fxi' )

    
    #plt.plot( xa, fx,   label='fx' )
    #plt.plot( xa, dXxi, label='dXxi' )


    plt.grid()
    plt.legend()
    plt.show()

