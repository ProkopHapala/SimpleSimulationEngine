#!/usr/bin/python

import numpy as np
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
    #    #print  " r %g fr %g = (f1 %g * e2 %g )+(e1 %g *f2 %g) r_2s %g r %g s %g " %(r[i], fr[i], f1[i],e2[i], e1[i],f2[i], r_2s[i], r[i], s ) 
    #    #print "Gauss::Coulomb r %g s %g E %g fr %g fs %g " %((r+0*s)[i],(s+0*r)[i], E[i], fr[i], fs[i]  )
    #    print "Gauss::Coulomb r,s %g,%g  f1,2 %g,%g e1,2 %g,%g ir %g " %( (r+0*s)[i],(s+0*r)[i], (f1+s*0)[i], (f2+s*0)[i], (e1+s*0)[i], (e2+s*0)[i],   (ir+s*0)[i]  )


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


'''
inline double erfx_e6( double x ){
    if( x>4.5 ){ return 1./x; }
    double xx = x*x;
    double even =  0.9850156202961753  +xx*(-0.02756061032579559  +xx*(-0.00188409579491924  +xx*(-0.003098629936170076 +xx*(-0.001348858853909826  +xx*(-3.98946569988845e-05 ) ) ) ) );
    double odd  = -0.13893350387140332 +xx*(-0.007664292475021448 +xx*( 0.003046826535877866 +xx*( 0.002879338499080343 +xx*( 0.0003260490382458129 +xx*( 1.97093650414204e-06 ) ) ) ) );
    double  t = even + x*odd;
    t*=t; t*=t; t*=t; // ^8
    return 1./( t + x );
}

inline double erfx_e9( double x ){
    if( x>4.5 ) return 1./x;
    double xx = x*x;
    if(x<1.){
        return 1.1283791662308296 +xx*(-0.3761262972953429 +xx*(0.1128363404233098 +xx*(-0.02685603827999912 +xx*(0.005192885862299865 +xx*(-0.0008053004722300972 +xx*(8.004020068129447e-05 ) ) ) ) ) );
    }
    double even =   0.9903386741213333  +xx*( 0.08180278811069948 +xx*( 0.219787883285348  +xx*( 0.0893543139653664  +xx*( 0.0071698531450102   +xx*( 8.644883946761633e-05 ) ) ) ) );
    double odd  =  -0.17511814497584813 +xx*(-0.2010794452848663  +xx*(-0.1692686167813105 +xx*(-0.03129254573733003 +xx*(-0.001037968593234627 +xx*(-3.164137211658646e-06 ) ) ) ) );
    double t = even + x*odd;
    t*=t; t*=t; t*=t; // ^8
    return 1./( t + x );
}
'''

def erfx_approx( x, s ):
    x=np.abs(x)
    ys   = 1./x 
    invs = 1/s
    x *= invs
    xx = x*x;
    even =  0.9850156202961753  +xx*(-0.02756061032579559  +xx*(-0.00188409579491924  +xx*(-0.003098629936170076 +xx*(-0.001348858853909826  +xx*(-3.98946569988845e-05 ) ) ) ) );
    odd  = -0.13893350387140332 +xx*(-0.007664292475021448 +xx*( 0.003046826535877866 +xx*( 0.002879338499080343 +xx*( 0.0003260490382458129 +xx*( 1.97093650414204e-06 ) ) ) ) );
    t    = even + x*odd;
    t*=t; t*=t; t*=t; # ^8
    mask = np.abs(x)<4.5
    ys[mask] = (invs/( t + x ))[mask]
    #ys[mask] = ( t )[mask]
    return ys

def derfx_approx( x, s ):
    '''
    df(t(x))/dx = (df/dt)*(dt/dx)
    df/dt = (  ) 
    '''
    x=np.abs(x)
    ys  =  1./x 
    dys  = -1./(x*x) 
    invs = 1/s
    x *= invs
    xx = x*x;
    even =  0.9850156202961753  +xx*(-0.02756061032579559  +xx*(-0.00188409579491924  +xx*(-0.003098629936170076 +xx*(-0.001348858853909826  +xx*(-3.98946569988845e-05 ) ) ) ) )
    odd  = -0.13893350387140332 +xx*(-0.007664292475021448 +xx*( 0.003046826535877866 +xx*( 0.002879338499080343 +xx*( 0.0003260490382458129 +xx*( 1.97093650414204e-06 ) ) ) ) )

    #deven =                           -0.02756061032579559*2  +xx*(-0.00188409579491924*4  +xx*(-0.003098629936170076*6 +xx*(-0.001348858853909826*8  +xx*(-3.98946569988845e-05*10 ) ) ) ) 
    #dodd  = -0.13893350387140332 +xx*(-0.007664292475021448*3 +xx*( 0.003046826535877866*5 +xx*( 0.002879338499080343*7 +xx*( 0.0003260490382458129*9 +xx*( 1.97093650414204e-06*11 ) ) ) ) )

    deven =                          -0.05512122065159118 +xx*(-0.00753638317967696 +xx*(-0.01859177961702045 +xx*(-0.01079087083127861  +xx*(-0.000398946569988845 ) ) ) )      
    dodd  = -0.1389335038714033 +xx*(-0.02299287742506434 +xx*( 0.01523413267938933 +xx*( 0.0201553694935624  +xx*( 0.002934441344212316 +xx*(2.168030154556244e-05 ) ) ) ) )

    print "%.16g %.16g %.16g %.16g %.16g      " %( -0.02756061032579559*2, -0.00188409579491924*4,  -0.003098629936170076*6, -0.001348858853909826*8, -3.98946569988845e-05*10 )
    print "%.16g %.16g %.16g %.16g %.16g %.16g" %(-0.13893350387140332   , -0.007664292475021448*3,  0.003046826535877866*5,  0.002879338499080343*7, 0.0003260490382458129*9, 1.97093650414204e-06*11 )

    t    =  even   + x*odd;
    dt   = deven*x +  dodd;
    t2 = t*t
    t4 = t2*t2
    t8 = t4*t4
    dt8_dx = 8*dt*t*t2*t4
    
    D = 1/(t8 + x)
    mask = np.abs(x)<4.5
    ys [mask] = (      invs           *D   )[mask]
    dys[mask] = (-invs*invs*(dt8_dx+1)*D*D )[mask]
    #ys[mask] = ( t )[mask]
    return ys, dys

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    x=np.linspace(0.001,6,1000)
    s = 1.0
    #y     = erfx_approx( x, s )
    y,dy     = derfx_approx( x, s )
    y_ref = spc.erf( x/s )/x
    #dy_ref = (  (2/np.sqrt(np.pi))*np.exp(-x*x/(s*s)) - spc.erf( x/s )  )/(x**2)
    dy_ref = (  (2/np.sqrt(np.pi))*np.exp(-x**2)*x - spc.erf( x )  )/(x**2)

    plt.plot(x,y    , label="erfx(x,s)")
    plt.plot(x,y_ref, label="erf(x/s)/x")

    x_= x[1:-1]
    dx=x[1]-x[0]
    #dy_num = (y_ref[2:]-y_ref[:-2])/(2.*dx)
    
    plt.plot(x      ,dy     , label="derfx(x,s)")
    plt.plot(x      ,dy_ref , label="derfx(x/s)/x")
    #plt.plot(x_,dy_num, label="derf(x/s)/x")

    plt.plot(x,(y -y_ref)*1e+6, label="error * 1e+6")
    plt.plot(x,(dy -dy_ref)*1e+6, label="derror * 1e+6")
    #plt.plot(x_,(dy[1:-1]-dy_num)*1e+6, label="derror * 1e+6")
    plt.legend(); plt.grid(); plt.show()
    exit(0)

    #s  = np.arange( 0.1, 5.0, 0.05 )
    #rs = np.arange( 0.1, 5.0, 0.05 )
    #S  = np.arange( 1.25, 5.0, 0.05 )
    #r  = 1.5 + 0.*s

    #ca = 1.0
    #cb = 1.0
    #cc = 1.0
    #cd = 1.0

    ca = 0.4
    cb = 0.3
    cc = 0.6
    cd = 0.8

    #ks = np.sqrt(0.5)
    sa = 1.2
    sb = 1.2
    sc = 1.2
    sd = 1.2

    n  = 40; x0 = 0.00001
    #n  = 10; x0 = 0.75
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
        #print "fromRho x %g  Fpi %g fx %g qij %g r %g "    %(   xa[i], Fpi[i],  fx[i], qij[i], r[i] ) ;
        #print "fromRho x %g  dSdp %g =(dSab %g * cab %g)" %(  xa[i], dSdp[i], dSab[i], 1 ) ;
        #print "fromRho x %g  dEdQ %g =(e %g * qj %g)" %(  xa[i], dEdQ[i], e[i], qj ) ;
        pass

    # ==== Derivative of Coulomb term without considering changes of Charges 
    #plt.plot( xa, e ,  label='e' )
    #plt.plot( xa, Fx,  label='dedx_ana' )
    #plt.plot( xs_, (e[1:]-e[:-1])/dx,':', label='dedx_num' )
    
    # ==== Derivative of Coulomb term with considering the Charges
    plt.plot( xa, E,  label='E' )
    plt.plot( xa, -F,  label='dEdx_ana' )
    #plt.plot( xs_, (E[1:]-E[:-1])/-dx,':', label='dEdx_num', lw=3 )
    plt.plot( xa[1:-1], (E[2:]-E[:-2])/(-2*dx),':', label='dEdx_num', lw=3 )
    #plt.plot( xa, fxi, label='fxi' )

    
    #plt.plot( xa, fx,   label='fx' )
    #plt.plot( xa, dXxi, label='dXxi' )


    plt.grid()
    plt.legend()
    plt.show()

