#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

import scipy.special as spc
import copy


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
    M_SQRT2   = 1.41421356237 
    M_SQRT1_2 = 0.70710678118
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
    fs =       -  e1f2 *r_s * is_ * M_SQRT1_2
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
    C    = e1*e2        # Overlap
    dCr  = C*(-2.*is2)  # derivative is correct, tested !
    #   TODO : How is it possible that derivative is the same as value (just rescaled) ???? 

    #double logC =  wxi*xi + wxj*xj - wx*X;
    #double C   = np.exp(-logC) * Ci * Cj

    #try:
    #    for i in range(len(r2)):
    #        print "product3D_s_deriv r %g s %g S %g dS %g " %(np.sqrt(r2[i]),s, S[i], dCr[i]  )
    #except: 
    #    pass
    return C,s,p, dCr*dp, (dSsi,dXsi,dXxi,dCsi), (dSsj,dXsj,dXxj,dCsj)

def getCoulombEF( r, si, sj, qi, qj, dSi=None, dA=None, ci=None, out=None, i=0, j=0 ):

    o = si*0 + r*0

    s =  combSize(si,sj);

    e, fx, fs = Coulomb( r, s ) 
    qij = qi * qj
    E   = e * qij

    # rhofS[i] -= fs*si;   rhofS[j] -= fs*sj;
    # double fsi = Fs*dssi - Fp.dot( dxsi );       # fromRho()

    Fp = fx * r  * qij    # pure derivative of coulombic forcefield
    Fs = fs * si * qij

    #print "qi,qj ", qi[0],qj[0]
    #if( dA is  None  ):
    #    print "e Fs %g qij %g Fsq %g "  %( , (fs*si+o)[0], (qij+o)[0], (Fs+o)[0] );

    #print "q(%g,%g)  E %g fs %g fr %g s %g r %g \n"  %( (qi+0*si)[0],(qj+0*si)[0], (e+0*si)[0], (fs+0*si)[0], (fx+0*si)[0], (s+0*si)[0], (r+0*si)[0] )

    if( dA is not None  ):
        (dSsi,dXsi,dXxi,dCsi) = dA
        eqj = e*qj
        #print "eqj %g E %g Fs %g dSsi %g dCsi %g cij %g "  %( (eqj+o)[0], (E+o)[0], (Fs+o)[0], (dSsi+o)[0], (dCsi+o)[0], (ci+o)[0] )
        Fs  = Fp*dXsi + Fs*dSsi  + eqj  *dCsi*ci 
        Fp  = Fp*dXxi            + eqj*2*dSi *ci  # total derivative due to charge change

        #Fs = ( fs*si*dSsa +  fx*r*dXsa ) * qi*qj   +  e*( dCsa*ci*qj )

        #print "e %g qij %g Fs %g " %( (e+o)[0], (qi*qj+o)[0], (Fs+o)[0]  )
        #print "e %g E %g s %g q %g r %g fx %g F %g dS %g dSr %g cc %g dEdQ %g " %( e[0], E[0], s[0], (qi*qj)[0], r[0], fx[0], F[0], (2*dSi*ci)[0], dSi[0], ci[0], eqj[0]  )
        #if out is not None:
        #    out[0] += eqj
        #    out[1] += Fq
        #    out[2] += Fp

    #print "e %g E %g s %g(%g,%g) q %g(%g,%g) r %g fx %g F %g dSi %g " %( e[0], E[0], s[0],si[0],sj[0], (qi*qj)[0],qi[0],qj[0], r[0], fx[0], F[0], 2*dSi*ci )
    
    print "[%i,%i] E %g qij %g Fs %g " %( i,j, (E+o)[0], (qi*qj+o)[0], (Fs+o)[0]  )

    return E,Fp, Fs
    #outs[0] += E
    #outs[1] += F 


def combSize(si,sj):
    return np.sqrt(si*si + sj*sj)

def evalEFtot( ecoef, esize, eXpos, xa=0., sa=0. ):
    o = xa*0 + sa*0
    if isinstance(xa,float):
        xa = eXpos[0][0] + o
    if isinstance(sa,float):
        sa = esize[0][0] + o

    eXpos = copy.deepcopy(eXpos);  eXpos[0][0] = xa
    esize = copy.deepcopy(esize);  esize[0][0] = sa

    qs = [[0,0,0],[0,0,0]]
    xs = [[0,0,0],[0,0,0]]
    ss = [[0,0,0],[0,0,0]]
    auxs = [None,None,None,None] 
    isqrt2 = np.sqrt(0.5)
    for io in range(2):
        for ib in range(2):
            xs[io][ib] = eXpos[io][ib]
            qs[io][ib] = ecoef[io][ib]**2
            ss[io][ib] = esize[io][ib]*isqrt2
        S, s, x, dS, dA, dB = product3D_s_deriv( esize[io][0]+o, eXpos[io][0]+o, esize[io][1]+o, eXpos[io][1]+o )
        cc = ecoef[io][0]*ecoef[io][1]*2
        xs[io][2] = x
        qs[io][2] = S*cc
        ss[io][2] = s
        auxs[io]  = ( dS, dA, cc )

    #out = [0.,0.]
    Etot  = 0
    Fptot = 0
    Fstot = 0
    # -- from Diagonal charges of orb #1 
    for i in range(2):
        for j in range(3):
            E,Fp,Fs = getCoulombEF( xs[0][i]-xs[1][j]+o,  ss[0][i],ss[1][j]+o, qs[0][i]+o, qs[1][j]+o, i=i,j=j )
            Etot += E
            #if (i==0)and(j==1):
            if (i==0):
                #Etot  += E
                Fptot += Fp
                Fstot += Fs
    # -- from overlap (off-diagonal) charges of orb #1
    i = 2
    out = [0.,0.,0.]
    for j in range(3):
    #for j in range(2,3):
        r = xs[0][i]-xs[1][j]
        E,Fp,Fs = getCoulombEF( r+o,  ss[0][i]+o,ss[1][j]+o, qs[0][i]+o, qs[1][j]+o, dSi=auxs[0][0], dA=auxs[0][1], ci=auxs[0][2], out=out, i=i,j=j )
        #print "[%i,%i] E %g r %g " %(i,j,(o+E)[0],(o+r)[0])
        Etot += E
        Fptot += Fp
        Fstot += Fs
    print " sum: dEdQ %g Fq %g F %g " %((o+out[0])[0],(o+out[1])[0],(o+out[2])[0])

    return Etot,Fptot,Fstot


def evalEF_off( xa, ecoef, esize, eXpos ):
    eXpos[0][0] = xa
   
    Sab, si, xab, dSab, dA, dB = product3D_s_deriv( esize[0][0], xa         , esize[0][1],eXpos[0][1] )
    Scd, sj, xcd, dScd, dC, dD = product3D_s_deriv( esize[1][0], eXpos[1][0], esize[1][1],eXpos[1][1] )
    
    #out = [0.,0.]
    
    ci = ecoef[0][0]*ecoef[0][1]*2
    cj = ecoef[1][0]*ecoef[1][1]*2
    qi = Sab*ci
    qj = Scd*cj

    E,Fr,Fs = getCoulombEF( xab-xcd, si, sj, qi, qj,dSi=dSab, dA=dA, ci=ci )

    return E,Fr

def evalEF_S_off( sa, ecoef, esize, eXpos ):
    #eXpos[0][0] = xa
    #esize[0][0] = sa

    Sab, si, xab, dSab, dA, dB = product3D_s_deriv( sa,          eXpos[0][0], esize[0][1],eXpos[0][1] )
    Scd, sj, xcd, dScd, dC, dD = product3D_s_deriv( esize[1][0], eXpos[1][0], esize[1][1],eXpos[1][1] )
    
    #out = [0.,0.]
    
    ci = ecoef[0][0]*ecoef[0][1]*2
    cj = ecoef[1][0]*ecoef[1][1]*2
    qi = Sab*ci
    qj = Scd*cj

    E,Fr,Fs = getCoulombEF( xab-xcd, si, sj, qi, qj, dSi=dSab, dA=dA, ci=ci )

    return E,Fs



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

    dx =  0.1
    xa =  np.arange( 0.01, 3.0, dx )
    xb =  0.0
    xc = -1.5
    xd =  0.0

    xs_ = (xa[1:]+xa[:-1])*0.5

    # overlaps
    Sab, si, xab, dQab, dA, dB = product3D_s_deriv( sa,xa, sb,xb )
    Scd, sj, xcd, dQcd, dC, dD = product3D_s_deriv( sc,xc, sd,xd )
    # coulomb
    s2        = si*si + sj*sj
    s         = np.sqrt(s2)
    r         = xab-xcd
    e, fx, fs = Coulomb( r, s )
    dXxi = dA[2] + xa*0 

    Qab = Sab*ca*cb*2; # factor 2 because offdiagonal is there twice 

    plt.plot( xa, r   , label='r'   )
    #plt.plot( xa, xab , label='Xab' )
    #plt.plot( xa, Sab , label='Sab' )
    plt.plot( xa, Qab , label='Qab' )
    #plt.plot( xa, dQab, label='dSab_ana' )
    #plt.plot( xs_, (Sab[1:]-Sab[:-1])/dx,':', label='dSab_num' )

    qij  = 4*Scd*Sab
    #qij  = Sab
    dQij = 4*Scd*dQab

    # Q: Why we dont need derivatives of charge ????
    #Fx = -fx*0.5*dA[1] # This works for zero initial distance between blobs
    Fx  = fx*r*dXxi
    Fpi = fx*r*qij # see 
    fxi = Fpi*dXxi
    
    print "Scd, 4*Scd ", Scd, 4*Scd
    print "For some reason each charge is scaled by 2.0"

    E = e*qij
    F = fxi + e*dQij   # total derivtive F = dE/dx = d(e*qi)/dx

    # Note:   e,fx=de/dx   are NOT multiplied by charge Qij
    #         Total force Fx = dE/dx = d(e*q)/dx = q*(de/dx) + e*(dq/dx)

    # ==== Derivative of Coulomb term without considering changes of Charges 
    #plt.plot( xa, e ,  label='e' )
    #plt.plot( xa, Fx,  label='dedx_ana' )
    #plt.plot( xs_, (e[1:]-e[:-1])/dx,':', label='dedx_num' )
    
    # ==== Derivative of Coulomb term with considering the Charges
    plt.plot( xa, E,  label='E' )
    #plt.plot( xa, F,  label='dEdx_ana' )
    #plt.plot( xs_, (E[1:]-E[:-1])/dx,':', label='dEdx_num', lw=3 )
    #plt.plot( xa, fxi, label='fxi' )

    
    #plt.plot( xa, fx,   label='fx' )
    #plt.plot( xa, dXxi, label='dXxi' )


    plt.grid()
    plt.legend()

    '''
    for i in range(len(r)):
        #print "Gauss::Coulomb r %g s %g E %g Fx %g fx %g " %(r[i], s, E[i], Fx[i], fx[i] )
        #print "fromRho r %g s %g E %g Fx %g fx %g " %((xa-xb)[i], s, E[i], Fx[i], fx[i] )
        #print "CoublombElement r %g s %g E %g fr %g qij %g  frq %g fij %g" %((xa-xb)[i], s, e[i], fx[i], qij[i], (fx*qij)[i], (fx*qij*r)[i] )
        #print "fromRho r %g s %g | E %g e %g qij %g(%g) | Fx %g Fpi %g dQij %g " %((xa-xb)[i], si, E[i],e[i]*2*Scd,qij[i],Sab[i], Fx[i], Fpi[i],dQij[i] )
        print "fromRho r %g  Eqi %g Cij %g | Fpi %g dXxi %g fxi %g Fxi %g "    %((xa-xb)[i],          e[i]*2*Scd,       Sab[i], Fpi[i], dXxi[i], fxi[i], F[i] );
        pass
    '''

    plt.show()

