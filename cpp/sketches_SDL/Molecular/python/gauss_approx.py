#!/usr/bin/python

import numpy as np
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


def GaussN( x2, w=1.0, npow=2 ):
    #x2 = (x2+x2*x2)/(1+19./16.*x2*x2)
    W2 = 1./(w*w*2**npow)
    g = 1.-x2*W2
    g[g<0]=0
    for i in range(npow):
        g*=g
    return g

exp_tay_coefs = [ 1., 1., 1./2, 1./6, 1./24, 1./120, 1./720, 1./5040, 1./40320, 1./362880, 1./3628800 ]

def exp_taylor( x, nTay, nPow ):
    Cs=exp_tay_coefs
    x*=1./2**nPow
    y  = 0.*x + Cs[0]
    xn = x.copy()
    for c in Cs[1:nTay+1]:
        y += c*xn
        xn*=x
    for i in range(nPow):
        y*=y
    #print x.shape, y.shape
    return y

def exp_taylor2( x, nTay, nPow ):
    Cs=exp_tay_coefs
    x*=0.5/2**nPow
    y  = 0.*x + Cs[0]
    xn = x.copy()
    for c in Cs[1:nTay]:
        y += c*xn
        xn*=x
    y*=( y + xn*Cs[nTay] )
    #y+=xn*Cs[nTay]
    #for i in range(nPow):
    #    y*=y
    #print x.shape, y.shape
    return y

def exp_taylor3( x, nTay, nPow ):
    Cs=exp_tay_coefs
    x*=0.33333333333/2**nPow
    y  = 0.*x + Cs[0]
    xn = x.copy()
    for c in Cs[1:nTay]:
        y += c*xn
        xn*=x
    y2 = y + xn*Cs[nTay]
    y *= (y2*y2)
    #y+=xn*Cs[nTay]
    #for i in range(nPow):
    #    y*=y
    #print x.shape, y.shape
    return y


if __name__ == "__main__":

    #xs = np.arange( -1.0, 6.0, 0.05 )
    xs = np.arange( -1.0, 8.0, 0.05 )
    x = xs.copy()
    w = 1.0

    #y_ref = np.exp( -(xs/w)**2 )
    
    n = 10
    colors = plt.cm.jet(np.linspace(0.,1.,n+1))

    xs = xs**2

    nTay = 3
    nPow = 3

    y_ref = np.exp( -xs )
    plt.plot( x ,y_ref,'k', lw=3, label="ref")
    for i in range(1,n,2):
        #ys    = GaussN( xs**2, w=1, npow=i )
        npow = n-i-1
        y    =    exp_taylor ( -xs, i, npow )
        #y2    =    exp_taylor2( -xs, i, 0 )
        #y3    =    exp_taylor3( -xs, i, 0 )
        #y13   =    exp_taylor ( -xs/(3.*w), i, 0 )**3
        #y12   =    exp_taylor ( -xs/(2.*w), i, 0 )**2
        yi    = 1./exp_taylor(   xs/w, i, npow )
        #y_err = ys - y_ref
        #plt.plot(x, y, '-',  c=colors[i], label="exp(-x)[%i,%i]"%(i,n-i) )
        #plt.plot(x, y, '-',  c=colors[i], label="exp (-x) tay %i" %i )
        #plt.plot(x, y2, ':',  c=colors[i], label="exp2(-x) tay %i" %i )
        #plt.plot(x, y3, '--',  c=colors[i], label="exp3(-x) tay %i" %i )
        #plt.plot(x, y2, '-',  c=colors[i], label="1/exp(-x) tay %i" %i )
        #plt.plot(x, yi, '-.',  c=colors[i], label="1/exp(x)[%i,%i]"%(i,n-i) )

        plt.plot(x, abs(y-y_ref), '-',  c=colors[i], label="exp(-x)[%i,%i]"%(i,npow) )
        #plt.plot(x, abs(y2-y_ref), ':',  c=colors[i], label="exp2(-x) tay %i" %i )
        #plt.plot(x, abs(y3-y_ref), '--', c=colors[i], label="exp3(-x) tay %i" %i )
        #plt.plot(x, abs(y12-y_ref), '-', c=colors[i], label="exp12(-x) tay %i" %i )
        #plt.plot(x, abs(y13-y_ref), '-', c=colors[i], label="exp13(-x) tay %i" %i )
        plt.plot(x, abs(yi-y_ref), '-.', c=colors[i], label="1/exp(x)[%i,%i]"%(i,npow) )
        #plt.plot(x, y_ref-ys, c=colors[i], label="npow %i" %i  )
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
    plt.ylim(1e-16,1.5)
    plt.legend()
    plt.grid()
    plt.show()






