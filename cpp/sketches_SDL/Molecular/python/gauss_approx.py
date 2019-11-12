#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

def GaussN( x2, w=1.0, npow=2 ):
    #x2 = (x2+x2*x2)/(1+19./16.*x2*x2)
    W2 = 1./(w*w*2**npow)
    g = 1.-x2*W2
    g[g<0]=0
    for i in range(npow):
        g*=g
    return g

if __name__ == "__main__":
    xs = np.arange( 0.0, 6.0, 0.05 )
    w = 1.0
    y_ref = np.exp( -(xs/w)**2 )
    
    for i in range(1,8):
        ys    = GaussN( xs**2, w=1, npow=i )
        #y_err = ys - y_ref
        plt.plot(xs, ys, ':', label="npow %i" %i  )
        plt.plot(xs, y_ref-ys, label="npow %i" %i  )
    '''
    ys    = GaussN( xs**2, w=1., npow=1 )
    plt.plot(xs, ys, label="approx 2" )
    plt.plot(xs,y_ref-ys, label="approx 2" )
    plt.plot(xs,y_ref,'k', lw=3, label="ref")
    '''
    #x2 = xs**2
    #x2_ = (x2+x2*x2)/(1+19./16.*x2*x2)
    #x2_ = (x2-x2*x2)/(1-(1/16.)*x2*x2)
    #plt.plot(xs, 1-x2_, label="x2_" )

    plt.yscale('log')
    plt.ylim(1e-16,1.0)
    plt.legend()
    plt.grid()
    plt.show()






