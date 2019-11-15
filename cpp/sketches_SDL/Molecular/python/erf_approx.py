#!/usr/bin/python

import numpy as np
import scipy.special as spc
import matplotlib.pyplot as plt

'''

Interesting Functions

exp(-x^2)
(exp(exp(-x^2))-1)/(e-1)
(exp(-exp(x^2)))*e
exp(-1/(1-x^2))*e   # Bump Function https://en.wikipedia.org/wiki/Bump_function

exp(-abs(x))
exp(-exp(abs(x)))*e

Approximate Gaussian like this
1-x^2(1.0-0.43*x^2+0.072*x^4)

'''

def erf_4(x):
    p = 1. + x*( 0.278393 + x*( 0.230389 + x*(0.000972 + x*0.078108 )))
    p=p*p; p=p*p;
    return 1. - 1./p;

def erf_6(x):
    p = 1. + x*( 0.0705230784 + x*( 0.0422820123 + x*( 0.0092705272 + x*( 0.0001520143 + x*( 0.0002765672 + x*0.0000430638 )))))
    p=p*p; p=p*p; p=p*p; p=p*p;
    return 1. - 1./p;

if __name__ == "__main__":

    #xs = np.arange( -1.0, 6.0, 0.05 )
    xs = np.arange( -1.0, 8.0, 0.05 )
    x = xs.copy()
    #n = 10
    #colors = plt.cm.jet(np.linspace(0.,1.,n+1))

    y_ref = spc.erf( x )/x
    y4    = erf_4  ( x )/x
    y6    = erf_6  ( x )/x

    plt.plot( x, y_ref,'k',  lw=3, label="ref")
    plt.plot( x, y4   ,':', label="y4")
    plt.plot( x, y6   ,':', label="y6")
    plt.plot( x, abs(y4-y_ref), '--',label="err4" )
    plt.plot( x, abs(y6-y_ref), '--',label="err6" )

    #plt.yscale('log')
    plt.ylim(1e-16,1.5)
    plt.legend()
    plt.grid()
    plt.show()






