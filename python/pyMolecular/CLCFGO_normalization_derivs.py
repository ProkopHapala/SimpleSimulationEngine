#!/usr/bin/python

import numpy as np
import scipy.special as spc
import CLCFGO_coulomb_derivs as clc


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

def potV( x, K=1.0 ):
    #E    = K*x*x 
    #dEdx = K*2*x
    E    = K*np.sin(x) 
    dEdx = K*np.cos(x)
    return E, dEdx

'''
def potV( x, K=1.0, s=0 ):
    x0=0
    E,fr,fs = clc.Coulomb( r, s )
    return E, fr
'''

def plotNumDeriv( xs, E, F, F_=None, title="" ):
    plt.figure()
    #dx   = xs[1]-xs[0]
    dx2   = xs[2:]-xs[:-2]
    xs_  = xs[1:-1]
    #Fnum = (E[2:]-E[:-2])/(2*dx)
    Fnum = (E[2:]-E[:-2])/(dx2)
    plt.plot( xs , E   ,'-k', label='E'    )
    plt.plot( xs , F   ,'-r', label='Fana' )
    plt.plot( xs_, Fnum,':y', label='Fnum' )
    if F_ is not None:
        plt.plot( xs , F_,'-m', label='Fana_' )
    plt.grid()
    plt.legend()
    plt.title(title)
    #plt.xlael()

def evalCharge( xa,xb,ca,cb ):
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
    return Q, (dQdxa,dQdxb,dQdca,dQdcb,dQdsa,dQdsb)

def evalEnergy( xa,xb,ca,cb ):
    Sab, si, xab, dSab, (dSsa,dXsa,dXxa,dS_dsa), (dSsb,dXsb,dXxb,dS_dsb) = clc.product3D_s_deriv( sa,xa, sb,xb )
    Va ,dVxa  = potV( xa  )
    Vb ,dVxb  = potV( xb  )
    Vab,dVxab = potV( xab )
    qaa = ca**2
    qbb = cb**2
    cab = ca*cb
    qab = 2*Sab*cab
    E     = qaa*Va   + qbb*Vb          + qab*Vab
    dEdxa = qaa*dVxa + cab*Vab*dSab*2  + qab*dVxab*dXxa
    dEdxb = qbb*dVxb - cab*Vab*dSab*2  + qab*dVxab*dXxb
    dEdca = 2*ca*Va  + 2*cb*Sab*Vab
    dEdcb = 2*cb*Vb  + 2*ca*Sab*Vab
    dEdsa = cab*Vab*dS_dsa*2 + qab*dVxab*dXsa
    dEdsb = cab*Vab*dS_dsb*2 + qab*dVxab*dXsb
    return E,(dEdxa,dEdxb,dEdca,dEdcb,dEdsa,dEdsb)

def outprojectNormalForce( dEdx, dQdx, E, Q=1 ):
    if isinstance(Q,int):
        return dEdx   - E*dQdx
    else:
        return dEdx/Q - E*dQdx/Q**2

def test( what="xa", ca=1.6, cb=-0.4, sa=0.35, sb=0.65, xa=-0.5, xb=+0.5, bNormalize=True ):
    if   what=="xa":
        xa =  np.arange( -2.0, 3.0, 0.01 );   xs  = xa
    elif what=="xb":
        xb =  np.arange( -2.0, 3.0, 0.01 );   xs  = xb
    elif what=="sa":
        sa =  np.arange( 0.5, 1.5, 0.01 );    xs  = sa
    elif what=="sb":
        sb =  np.arange( 0.5, 1.5, 0.01 );    xs  = sb
    elif what=="ca":
        ca =  np.arange( -2.0, 2.0, 0.01 );   xs  = ca.copy()
    elif what=="cb":
        cb =  np.arange( -2.0, 2.0, 0.01 );   xs  = cb.copy()
    if bNormalize:
        Q,_ = evalCharge( xa,xb,ca,cb )
        rescale = 1./np.sqrt(Q)
        ca *= rescale
        cb *= rescale
    Q,(dQdxa,dQdxb,dQdca,dQdcb,dQdsa,dQdsb) = evalCharge( xa,xb,ca,cb )
    E,(dEdxa,dEdxb,dEdca,dEdcb,dEdsa,dEdsb) = evalEnergy( xa,xb,ca,cb )
    Sab, si, xab, dSab, (dSsa,dXsa,dXxa,dS_dsa), (dSsb,dXsb,dXxb,dS_dsb) = clc.product3D_s_deriv( sa,xa, sb,xb )
    if   what=="xa":
        plotNumDeriv( xs, E, dEdxa, title="E(xa)" )
    elif what=="xb":
        plotNumDeriv( xs, E, dEdxb, title="E(xb)" )
    elif what=="sa":
        plotNumDeriv( xs, E, dEdsa, title="E(sa)" )
    elif what=="sb":
        plotNumDeriv( xs, E, dEdsb, title="E(sb)" )
    elif what=="ca":
        plotNumDeriv( xs, E, dEdca, title="E(ca)" )
    elif what=="cb":
        plotNumDeriv( xs, E, dEdcb, title="E(cb)" )


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    #x0 = 0.00001; x1 = 0.00001; dx =  0.05

    ca = 1.6
    cb = -0.4
    sa = 0.35
    sb = 0.65
    xa =  -0.5
    xb =  +0.5

    #xa =  np.arange( -2.0, 3.0, 0.01 );     xs  = xa
    #xb =  np.arange( -2.0, 3.0, 0.01 );     xs  = xb
    #ca =  np.arange( -2.0, 2.0, 0.01 );     xs  = ca
    ca =  np.arange( -2.0, 2.0, 0.01 );     xs  = ca.copy()
    #cb =  np.arange( -2.0, 2.0, 0.01 );     xs  = cb.copy()
    #sa =  np.arange( 0.5, 1.5, 0.01 );     xs  = sa
    #sb =  np.arange( 0.5, 1.5, 0.01 );     xs  = sb


    #print "xs ", xs

    # --- Charge constrain 1=Q=<psi|psi> 

    #bNormalize = False
    bNormalize = True

    if bNormalize:
        Q,_ = evalCharge( xa,xb,ca,cb )
        rescale = 1./np.sqrt(Q)
        ca *= rescale
        cb *= rescale

    Q,(dQdxa,dQdxb,dQdca,dQdcb,dQdsa,dQdsb) = evalCharge( xa,xb,ca,cb )
    E,(dEdxa,dEdxb,dEdca,dEdcb,dEdsa,dEdsb) = evalEnergy( xa,xb,ca,cb )

    Sab, si, xab, dSab, (dSsa,dXsa,dXxa,dS_dsa), (dSsb,dXsb,dXxb,dS_dsb) = clc.product3D_s_deriv( sa,xa, sb,xb )
    #plotNumDeriv( xs, xab, dXxa+xs*0, title="xab(xa)" )

    #plotNumDeriv( xs, Q, dQdxa, title="Q(xa)" )
    #plotNumDeriv( xs, Q, dQdxb, title="Q(xb)" )
    #plotNumDeriv( xs, Q, dQdca, title="Q(ca)" )
    #plotNumDeriv( xs, Q, dQdcb, title="Q(cb)" )

    #plotNumDeriv( xs, Sab, dS_dsb, title="Sab(sb)" )
    #plotNumDeriv( xs, Q, dQdsb, title="Q(sb)" )
    #plotNumDeriv( xs, E, dEdsb, title="E(sb)" )

    #plotNumDeriv( xs, Sab, dS_dsa, title="Sab(sa)" )
    #plotNumDeriv( xs, Q, dQdsa, title="Q(sa)" )
    #plotNumDeriv( xs, E, dEdsa, title="E(sa)" )

        
    if bNormalize:
        Q = 1; E_ = E
        #E_ = E/Q; 
        dEdxa_ = outprojectNormalForce( dEdxa, dQdxa, E, Q )
        dEdxb_ = outprojectNormalForce( dEdxb, dQdxb, E, Q )
        dEdca_ = outprojectNormalForce( dEdca, dQdca, E, Q )
        dEdcb_ = outprojectNormalForce( dEdcb, dQdcb, E, Q )
        dEdsa_ = outprojectNormalForce( dEdsa, dQdsa, E, Q )
        dEdsb_ = outprojectNormalForce( dEdsb, dQdsb, E, Q )
        #print "xs ", xs
        #print "E ", E
        #plotNumDeriv( xs, E_, dEdsa_, title="E_(sa)" )
        #plotNumDeriv( xs, E_, dEdsb_, title="E_(sb)" )
        #plotNumDeriv( xs, E_, dEdxa_, title="E_(xa)" )
        #plotNumDeriv( xs, E_, dEdxb_, title="E_(xb)" )
        plotNumDeriv( xs, E_, dEdca_*rescale, title="E_(ca)" )
        #plotNumDeriv( xs, E_, dEdcb_*rescale, title="E_(cb)" )
        #dcadx = (ca[2:]-ca[:-2])/(xs[2:]-xs[:-2])
        #plotNumDeriv( xs, cb, ca, title="cb_(cb)" ); 
        #plt.plot(xs, 2*Sab*ca*cb, label="qab" ); plt.legend()
        #plt.plot(xs, ca**2 + cb**2 + 2*Sab*ca*cb, label="Qtot" ); plt.legend()
        #plotNumDeriv( xs[1:-1], E_[1:-1], dEdcb_[1:-1]*dcadx, title="E_(cb)" )
       
    else:    
        #plotNumDeriv( xs, E, dEdsa, title="E(sa)" )
        #plotNumDeriv( xs, E, dEdsb, title="E(sb)" )
        #plotNumDeriv( xs, E, dEdxa, title="E(xa)" )
        #plotNumDeriv( xs, E, dEdxb, title="E(xb)" )
        #plotNumDeriv( xs, E, dEdca, title="E(ca)" )
        #plotNumDeriv( xs, E, dEdcb, title="E(cb)" )
        pass


    plt.legend()
    plt.show()

