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

'''

def potV( x, K=1.0 ):
    #E    = K*x*x 
    #dEdx = K*2*x

    E    = K*np.sin(x) 
    dEdx = K*np.cos(x)
    return E, dEdx


def plotNumDeriv( xs, E, F, F_=None, title="" ):
    plt.figure()
    dx   = xs[1]-xs[0]
    xs_  = xs[1:-1]
    Fnum = (E[2:]-E[:-2])/(2*dx)
    plt.plot( xs , E   ,'-k', label='E'    )
    plt.plot( xs , F   ,'-r', label='Fana' )
    plt.plot( xs_, Fnum,':y', label='Fnum' )
    if F_ is not None:
        plt.plot( xs , F_,'-m', label='Fana_' )
    plt.grid()
    plt.legend()
    plt.title(title)
    #plt.xlael()


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    #x0 = 0.00001; x1 = 0.00001; dx =  0.05

    ca = 0.5
    cb = 0.5
    sa = 0.5
    sb = 0.5

    xa =  -0.5
    xb =  +0.5

    #xa =  np.arange( x0, x0+dx*n, dx )
    xa =  np.arange( -2.0, 3.0, 0.01 )

    ts = (xa + xb + ca + cb)*0

    xs  = xa
    xs_ = xs[1:-1]

    Sab, si, xab, dSab, (dSsa,dXsa,dXxa,dS_dsa), (dSsb,dXsb,dXxb,dS_dsb) = clc.product3D_s_deriv( sa,xa, sb,xb )

    cab = ca*cb
    qaa = ca*ca
    qbb = ca*cb
    qab = 2*Sab*ca*cb

    # --- Charge constrain 1=Q=<psi|psi> 
    Qtot  =  qaa + qbb + qab
    dQdxa =  cab*dSab*2
    dQdxb = -cab*dSab*2
    dQdca = 2*ca + 2*Sab*cb
    dQdcb = 2*cb + 2*Sab*ca

    renorm = 1/np.sqrt(Qtot)
    ca*=renorm
    cb*=renorm
    cab = ca*cb

    # --- Total energy functional
    Va ,dVxa  = potV( xa  )
    Vb ,dVxb  = potV( xb  )
    Vab,dVxab = potV( xab )
    Etot  = qaa*Va   + qbb*Vb        + qab*Vab
    dEdxa = qaa*dVxa + cab*Vab*dSab*2  + qab*dVxab*dXxa
    dEdxb = qbb*dVxa - cab*Vab*dSab*2  + qab*dVxab*dXxb
    dEdca = 2*ca*Va  + 2*cb*Sab*Vab
    dEdcb = 2*ca*Vb  + 2*ca*Sab*Vab

    #print dEdcb

    #AdQ = (dQdxa**2    + ts) + (dQdxb**2     + ts) + (dQdca**2     + ts) + (dQdcb**2     + ts)
    #dQE = (dQdxa*dEdxa + ts) + (dQdxb**dEdxb + ts) + (dQdca**dEdca + ts) + (dQdcb**dEdcb + ts)

    AdQ = dQdxa**2    +  dQdxb**2   + dQdca**2    + dQdcb**2
    dQE = dQdxa*dEdxa + dQdxb*dEdxb + dQdca*dEdca + dQdcb*dEdcb 
    C   = -dQE/np.sqrt(AdQ) 

    C *=0.3

    dEdxa_ = dEdxa + C*dQdxa
    dEdxb_ = dEdxb + C*dQdxb
    dEdca_ = dEdca + C*dQdca
    dEdcb_ = dEdcb + C*dQdcb
    
    #Etot  = qab*Vab
    #dEdxa = cab*Vab*dSab*2  + qab*dVxab*dXxa
    #dEdxb =-cab*Vab*dSab*2  + qab*dVxab*dXxb
    #dEdca = cb*Sab*Vab
    #dEdcb = ca*Sab*Vab

    #plotNumDeriv( xa, ca, Qtot, "Qtot(ca)" )
    #plotNumDeriv( xa, Sab, dSab, "Qtot(xa)" )
    #plotNumDeriv( xa, xab, dXxa+xa*0, "xab" )
    #plotNumDeriv( xa, Qtot, dQdxa, "Qtot(xa)" )
    plotNumDeriv( xa, Etot, dEdxa, dEdxa_, "Etot(xa)" )
    #plotNumDeriv( xa, Etot, dEdxa, dQdxa, "Etot(xa)" )

    #plt.figure(); plt.plot( xa, C, label="C(xa)" )

    plt.legend()
    plt.show()

