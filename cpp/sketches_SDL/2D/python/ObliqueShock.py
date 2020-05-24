#!/usr/bin/python

import numpy as np
import sympy as sy


'''

https://en.wikipedia.org/wiki/Oblique_shock#The_%CE%B8-%CE%B2-M_equation

tan(Theta) = 2*cotg(beta)*( M**2 sin(beta)**2 - 1 )/( M**2*(kappa + cos(2*beta)) + 2 )

cos(2*beta) = 1- 2*sin(beta)**2   #  https://mathworld.wolfram.com/TrigonometricAdditionFormulas.html

T = 2*(c/s)*( MM*s*s - 1 )/( MM*(k+1-2*s*s) +2 )
T*( MM*(k+1-2*s*s) +2 )*s = 2*sqrt(1-s*s)*( MM*s*s - 1 )

(2+MM*(k+1)+2)*T*s -  2*MM*T*s*s*s - 2*sqrt(1-s*s)*( MM*s*s - 1 )  = 0



-0.5 *   ( 2 - 2*MM*s*s )/( (MM*(k+1)+2) - 2*MM*s*s )

C =  (MM*(k+1)+2)

-0.5 *   ( 2 - C + C - 2*MM*s*s )/( ( C - 2*MM*s*s )

-0.5 * ( 1 + ( 2 - C )/( C - 2*MM*s*s ) )


'''

'''
s, MM, T, K = sy.symbols( 's MM T K' )

Eq_s = T*( MM*( K + 1 - 2*s*s ) + 2 )*s - 2*sy.sqrt(1-s*s)*( MM*s*s - 1 )


print sy.solve( Eq_s, s )
'''

def thetaOblique_cs( c,s,  M=1.5, k=1.4 ):
    MM=M*M
    return 2*(c/s)*( MM*s**2 - 1 )/( MM*(k+1-2*s**2) +2 )


def tgThetaOblique_cs_new( c,s,  M=1.5, k=1.4 ):
    MM  = M*M
    ctg = c/s
    C   = MM*(k+1) + 2
    #return (c/s)*( 2*MM*s**2 - 2 )/( MM*(k+1-2*s**2) +2 )
    #return (c/s)*( 2*MM*s**2 - 2 )/( C - 2*MM*s*s )
    #return -(c/s)*(  2 - 2*MM*s**2 )/( C - 2*MM*s*s )
    #return -(c/s)*(  2-C+C - 2*MM*s**2 )/( C - 2*MM*s*s )
    return -(c/s)*(1 + (2-C)/( C - 2*MM*s*s ) )

def thetaOblique_cs_new( c,s,  M=1.5, k=1.4 ):
    return np.arctan( tgThetaOblique_cs_new( c,s,  M=M, k=k ) )

if __name__ == "__main__":
    rad2deg = 180./np.pi
    import matplotlib.pyplot as plt
    #xs = np.linspace(0, np.pi,1000)
    xs = np.linspace(0, np.pi/2,1000)
    c  = np.cos(xs)
    s  = np.sin(xs)
    M = 1.5
    MM = M**2
    #nDOFs = 5
    #k = 1 + 2/(nDOFs) 
    #k = 5./3. # diatomic
    #k  = 7./5. # monoatomic
    #tgT = 2*(c/s)*( MM*s**2 - 1 )/( MM*(k+1-2*s**2) +2 )
    #plt.plot( xs, tgT )
    f = c/s
    #plt.plot( xs, f )
    #plt.plot( xs, thetaOblique_cs( c,s,  M=1.1, k=1.4 )*f, label='M=1.1' )
    #plt.plot( xs, thetaOblique_cs    ( c,s,  M=1.5, k=1.4 )/f, label='M=1.5' )
    #plt.plot( xs, thetaOblique_cs_new( c,s,  M=1.5, k=1.4 )/f, label='M=1.5new' )
    
    # compare to wiki 
    # https://en.wikipedia.org/wiki/Oblique_shock#/media/File:ObliqueShockAngleRelation.png
    #Ms = [1.1,1.3,2.0,3.0,6.0]
    #Ms = [3.0]
    Ms = np.exp( np.linspace( 0.,1.,10 )*np.log(10) )
    for M in Ms:
        #plt.plot( xs*rad2deg, rad2deg*thetaOblique_cs_new( c,s,  M=M, k=1.4 ), label=('M=%1.1f' %M) )
        ys = thetaOblique_cs_new( c,s,  M=M, k=1.4 )
        i  = np.argmax(ys)
        i0 = np.searchsorted( ys, -0.5 )
        ymax = ys[i]
        xmax = xs[i]
        #ys[ys<0] = 1000
        #i0 = np.argmin(ys)
        #x0   = xs[i0]; print x0
        i-=50
        st_ = np.sin(ys)
        st  = st_[i0:i]
        s_  = s[i0:i] 
        cfs   = np.polyfit( st , s_, 6 )
        ypoly = np.polyval( cfs, st )
        plt.plot( st_, s, label=('M=%1.1f' %M) )
        plt.plot( st , ypoly, '--', label=('pM=%1.1f' %M) )
        
        plt.plot( st , (s_-ypoly)*1000, ':', label=('pM=%1.1f' %M) )
        
        #plt.plot( np.sin(ys), c, label=('M=%1.1f' %M) )
        #plt.plot( np.sin(ys), xs*rad2deg, label=('M=%1.1f' %M) )
        #plt.plot( np.cos(ys), xs*rad2deg, label=('M=%1.1f' %M) )
        #plt.plot( ys*rad2deg, xs*rad2deg, label=('M=%1.1f' %M) )
        #plt.plot( [ymax*rad2deg],[xmax*rad2deg],'o' )
        #x_ =  (xs-x0)/(xmax-x0)
        #x_ = 1-x_
        #y_ = 1-ys/ymax
        #plt.plot( y_, x_/np.sqrt(y_), label=('M=%1.1f' %M) )
    #plt.plot( xs, c, 'b:', label='cos' )
    #plt.plot( xs, s, 'r:', label='sin' )
    plt.legend()
    #plt.ylim(0.,50); plt.xlim(0.,90)
    #plt.ylim(0.,90); plt.xlim(-1.,1)
    plt.ylim(-1.,2.); plt.xlim(-1.,2)
    #plt.ylim(0.,1.5); plt.xlim(0.,1.)
    #plt.ylabel( "tan(theta)" )
    #plt.ylabel( "tan(theta) / cotan(beta)" )
    #plt.ylabel( "tan(theta)" )
    plt.ylabel( "beta [rad]" ); plt.xlabel( "theta [deg]" )
    plt.grid()
    plt.show()
