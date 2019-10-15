
import numpy as np

def genPoints( n, sc, dphi0=0.0 ):
    ps = np.empty((n+1,2))
    phiMax = 2*np.pi
    dPhi   = phiMax/n
    #phis   = np.arange(0,phiMax,phiMax/n) + dPhi
    phis   = np.linspace(0,phiMax,n+1) + dphi0*dPhi
    ps[:,0] = np.cos(phis)*sc
    ps[:,1] = np.sin(phis)*sc
    return ps


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    np.set_printoptions(precision=16, linewidth=200 )

    ps  = []
    scs = []
    sc  = 1.0

    ps.append( np.array([[0.0,0.0],]) ); scs.append(sc*2)
    ps.append( genPoints(6 , sc      ) ); 
    scs.append(sc*2)
    sc*=2.0

    ps.append( genPoints(6 , sc*0.9, 0.5       ) );    scs.append(sc)
    ps.append( genPoints(6 , sc*0.95,     ) );  scs.append(sc*0.9)
    sc*=1.53

    for i in xrange(4):
        ps.append( genPoints(12, sc, 0.5*(i+1) ) ); 
        scs.append(sc)
        sc*=1.53

    plt.figure(figsize=(5,5))
    for i,psi in enumerate(ps):
        print i
        plt.plot(psi[:,0],psi[:,1], "o", markersize=(4.5*scs[i]), fillstyle='none' )

    plt.axis('equal')
    plt.show()