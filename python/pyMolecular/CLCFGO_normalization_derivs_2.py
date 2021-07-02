#!/usr/bin/python

import numpy as np
import scipy.special as spc
import CLCFGO_coulomb_derivs as clc
import matplotlib.pyplot as plt

'''

Forces for constrained minimization
===================================
We variationally optimize oparameters (ca,cb,xa,xb) 
by moving these parameters along energy gradient dE/d(ca,cb,xa,xb)
with constrain Q=1, to make sure these forces dE/d(ca,cb,xa,xb) does not break the constrain 
we have to subtract component in direction of dQ/d(ca,cb,xa,xb)


Lets have simple orbital composed of two gaussians  
Phi_1(x) = ca * Ga(x) + cb * Gb(x)
where Ga(x)=G(x-xa) and Gb(x)=G(x-xb)

Electron denity from this orbital is 
Rho_1(x) = Phi_1(x)^2 
         = ca^2*Ga^2 + cb^2*Gb^2 +  2*cb*cb*Ga*Gb

Total charge is:
Q(ca,cb,xa,xb) = ca^2 + cb^2 + Sab(xa-xb)*ca*cb
               = qaa  + qbb  + qab  
where overlap 
Sab(xa-xb) = <G(x-xa)|G(x-xb)>

For simplicity define total energy operating only on charges
E(ca,cb,xa,xb) = qaa *V(xa) + qbb**V(xb) +  qab             *V(xab(xa,xb))
               = ca^2*V(xa) + cb^2*V(xb) +  Sab(xa-xb)*ca*cb*V(xab(xa,xb))

Charge constrain derivatives dQ/d(ca,cb,xa,xb):
dQ/dca = 2*ca + Sab*cb
dQ/dcb = 2*cb + Sab*ca
dQ/dxa = ca*cb*dSab/dxa
dQ/dxb = ca*cb*dSab/dxb

dE/dca = 2*ca*V(xa) + Sab*cb*V(xab(xa,xb))
dE/dcb = 2*cb*V(xb) + Sab*ca*V(xab(xa,xb))
dE/dxa = ca^2*dV(xa)/dxa + ca*cb*( V(xab)*(dSab/dxa) + Sab*(dV(xab)/dxa) )
dE/dxb = cb^2*dV(xb)/dxb + ca*cb*( V(xab)*(dSab/dxb) + Sab*(dV(xab)/dxb) )

Ei     = Hii/Sii =  <psi_i|H|psi_i>/<psi_i|psi_i>
dEi/dx = d(Hii/Sii) = (1/Sii)*(dH/dx) + (Hii/Sii^2)*(dS/dx)
Because we re-normalize wavefunctions in each step, we can assume S=1, which we substitute to obtain
dEi/dx = (dH/dx) + Hii*(dS/dx) = (dH/dx) + Ei*(dS/dx)
'''

def potV( x, s=0, K=-1.0, x0=0, s0=0.2, y0=0.2 ):
    #print " s ", s
    x_   = x-x0
    r    = np.sqrt( x_**2 + y0**2  )
    dr_x = x_/r
    #s_   = np.sqrt( s**2 + s0**2 )
    #ds_s = s/s_
    #E,fr,fs = clc.Coulomb( r, s )
    E,fr,fs = clc.Coulomb_new(r, s )
    return E*K, fr*r*dr_x*K, fs*-s*K

def plotNumDeriv( xs, E, F, F_=None, title="", bNewFig=True ):
    if bNewFig:
        plt.figure()
    #dx   = xs[1]-xs[0]
    dx2   = xs[2:]-xs[:-2]
    xs_  = xs[1:-1]
    #Fnum = (E[2:]-E[:-2])/(2*dx)
    Fnum = (E[2:]-E[:-2])/(dx2)
    plt.plot( xs , E   ,'--k', lw=3, label='E'    )
    plt.plot( xs , F   ,'-r', label='Fana' )
    plt.plot( xs_, Fnum,':y', label='Fnum' )
    if F_ is not None:
        plt.plot( xs , F_,'-m', label='Fana_' )
    plt.grid()
    plt.legend()
    plt.title(title)
    #plt.xlael()

def evalCharge(  A, B ):
    (xa,sa,ca) = A
    (xb,sb,cb) = B
    Sab, si, xab, dSab, (dSsa,dXsa,dXxa,dS_dsa), (dSsb,dXsb,dXxb,dS_dsb) = clc.product3D_s_deriv( sa,xa, sb,xb )
    cab = ca*cb
    qaa = ca*ca
    qbb = cb*cb
    qab = 2*Sab*cab
    Q     =  qaa + qbb + qab
    dQdxa =  cab*dSab*2
    dQdxb = -cab*dSab*2
    dQdca = 2*ca + 2*Sab*cb
    dQdcb = 2*cb + 2*Sab*ca
    dQdsa = 2*cab*dS_dsa
    dQdsb = 2*cab*dS_dsb
    return Q, (dQdxa,dQdsa,dQdca),(dQdxb,dQdcb,dQdsb)

def evalEnergy( A, B ):
    M_SQRT1_2 = 1/np.sqrt(2)
    (xa,sa,ca) = A
    (xb,sb,cb) = B
    Sab, sab , xab, dSab, (dSsa,dXsa,dXxa,dS_dsa), (dSsb,dXsb,dXxb,dS_dsb) = clc.product3D_s_deriv( sa,xa, sb,xb )
    Va ,dVxa ,dVsa  = potV( xa,  sa*M_SQRT1_2  )
    Vb ,dVxb ,dVsb  = potV( xb,  sb*M_SQRT1_2  )
    Vab,dVxab,dVsab = potV( xab, sab )
    qaa = ca**2
    qbb = cb**2
    cab = ca*cb
    qab = 2*Sab*cab
    
    E     = qaa*Va   + qbb*Vb          + qab*Vab
    dEdxa = qaa*dVxa + cab*Vab*dSab*2  + qab*dVxab*dXxa
    dEdxb = qbb*dVxb - cab*Vab*dSab*2  + qab*dVxab*dXxb
    dEdca = 2*ca*Va  + 2*cb*Sab*Vab
    dEdcb = 2*cb*Vb  + 2*ca*Sab*Vab
    dEdsa = cab*Vab*dS_dsa*2 + qab*dVxab*dXsa  + qab*dVsab*dSsa + qaa*dVsa
    dEdsb = cab*Vab*dS_dsb*2 + qab*dVxab*dXsb  + qab*dVsab*dSsb + qaa*dVsb
    
    '''
    E     = qab*Vab
    dEdxa = +cab*Vab*dSab*2  + qab*dVxab*dXxa
    dEdxb = -cab*Vab*dSab*2  + qab*dVxab*dXxb
    dEdca = + 2*cb*Sab*Vab
    dEdcb = + 2*ca*Sab*Vab
    dEdsa = cab*Vab*dS_dsa*2 + qab*dVxab*dXsa  + qab*dVsab*dSsa
    dEdsb = cab*Vab*dS_dsb*2 + qab*dVxab*dXsb  + qab*dVsab*dSsb
    '''
    '''
    E     = qaa*Va      +Sab*0
    #for i in range(len(xa)):
    #    #print "py xa %g E %g  qaa %g Va %g " %(xa[i], E[i], qaa, Va[i] )
    #    print "py xa %g E %g  e %g qij %g qi %g qj %g " %(xa[i], E[i], Va[i], qaa, ca, ca )
    dEdxa = qaa*dVxa    +Sab*0
    dEdxb = 0           +Sab*0
    dEdca = 2*ca*Va     +Sab*0
    dEdcb = 0           +Sab*0
    dEdsa =  + qaa*dVsa +Sab*0
    dEdsb = 0           +Sab*0
    '''
    return E,(dEdxa,dEdsa,dEdca),(dEdxb,dEdsb,dEdcb)

def outprojectNormalForce( dEdx, dQdx, E, Q=1 ):
    if isinstance(Q,int):
        return dEdx   - E*dQdx
    else:
        return dEdx/Q - E*dQdx/Q**2

def evalTest( what="xa", bNormalize=True,     xa=-0.4,sa=0.35,ca=1.6,     xb=+0.5,sb=0.55,cb=-0.4  ):
    if what=="xa":
        #xa  =  np.arange( -2.0, 3.0, 0.01 );   xs  = xa
        xa  =  np.arange( -2.0, 3.0, 0.1 );   xs  = xa
    elif what=="sa":
        sa =  np.arange(  0.25, 2.0, 0.01 );   xs  = sa
    elif what=="ca":
        ca =  np.arange( -2.0, 2.0, 0.01 );    xs  = ca.copy()

    if bNormalize:
        Q,_,_ = evalCharge( [xa,sa,ca], [xb,sb,cb] )
        rescale = 1./np.sqrt(Q)
        ca *= rescale
        cb *= rescale
    A=[xa,sa,ca]; B=[xb,sb,cb]

    Q,(dQdxa,dQdsa,dQdca),(dQdxb,dQdcb,dQdsb) = evalCharge( A,B )
    E,(dEdxa,dEdsa,dEdca),(dEdxb,dEdsb,dEdcb) = evalEnergy( A,B )

    if bNormalize:
        Q = 1; E_ = E
        dEdxa = outprojectNormalForce( dEdxa, dQdxa, E, Q )
        dEdsa = outprojectNormalForce( dEdsa, dQdsa, E, Q )
        dEdca = outprojectNormalForce( dEdca, dQdca, E, Q )*rescale
        #dEdxb = outprojectNormalForce( dEdxb, dQdxb, E, Q )
        #dEdsb = outprojectNormalForce( dEdsb, dQdsb, E, Q )
        #dEdcb = outprojectNormalForce( dEdcb, dQdcb, E, Q )*rescale
    
    return Q,E, (dEdxa,dEdsa,dEdca),(dQdxa,dQdsa,dQdca),xs

def run_test( what="xa", bNormalize=True, xa=-0.4,sa=0.35,ca=1.6,     xb=+0.5,sb=0.55,cb=-0.4 ):
    Q,E, (dEdxa,dEdsa,dEdca),(dQdxa,dQdsa,dQdca),xs= evalTest( what=what, bNormalize=bNormalize,     xa=xa,sa=sa,ca=ca,     xb=xb,sb=sb,cb=cb  )
    #Sab, sab, xab, dSab, (dSsa,dXsa,dXxa,dS_dsa), (dSsb,dXsb,dXxb,dS_dsb) = clc.product3D_s_deriv( sa,xa, sb,xb )
    #Va ,dVxa ,dVsa  = potV( xa, s=sa )
    #Vab,dVxab,dVsab = potV( xab, sab )
    if   what=="xa":
        #plotNumDeriv( xs, Va, dVxa, title="Va(xa)" )
        #plotNumDeriv( xs, Vab, dVxab*dXxa, title="Vab(xa)" )
        plotNumDeriv( xs, E, dEdxa, title="E(xa)" )
    elif what=="sa":
        #plotNumDeriv( xs, Sab, dS_dsa, title="Sab(sa)" )
        #plotNumDeriv( xs, xab, dXsa  , title="Xab(sa)" )
        #plotNumDeriv( xs, sab, dSsa  , title="sab(sa)" )
        #plotNumDeriv( xs, Va,  dVsa  , title="Va(sa)" )
        #plotNumDeriv( xs, Vab, dVsab*dSsa + dVxab*dXsa, title="Vab(sa)" )
        plotNumDeriv( xs, E, dEdsa, title="E(sa)" )
    elif what=="ca":
        plotNumDeriv( xs, E, dEdca, title="E(ca)" )

if __name__ == "__main__":
    bNormalize=False
    #bNormalize=True
    run_test( what="xa", bNormalize=bNormalize )
    run_test( what="sa", bNormalize=bNormalize )
    run_test( what="ca", bNormalize=bNormalize )
    plt.legend()
    plt.show()

