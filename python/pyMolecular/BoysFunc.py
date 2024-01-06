
import numpy as np

from scipy.special import hyp1f1
from scipy.special import erf as erf_scipy
#import sympy.mpmath as mpm
#import mpmath as mpm

'''


NOTE:
BoysF is used for Coulomb integrals like :

Vpq = Integral{   exp(-p*(r1)^2) * exp(-q*(r2)^2)  /|r1-r2| }

Vpq = 2*pi^(5/2) / ( p*q*sqrt(p+q) ) * BoysF( a*r^2 )

F0() is related to error-function erf()

F0(x) = sqrt(pi/(4x)) * erf( sqrt(x) ) 

Vpq = (pi/sqrt(p*q))^3 * erf( sqrt(a) * r ) / r

=> It makes sense to tabulate rather F0(x^2)

'''






import splines

# Boys Function (critical for calculation of Gaussian integrals)
#   http://shivupa.github.io/blog/efficient-evaluation-of-the-boys-function/

np.set_printoptions(precision=16, linewidth=200)

# ================== Functions

def  KummerF(a,b,z,order=10):
    '''
    #   https://en.wikipedia.org/wiki/Confluent_hypergeometric_function#Kummer.27s_equation
    F(a,b,z) = Sum_i{ rfac(a,n) * z^n / ( rfac(b,n) * n! )  }
    where 'rfac()' is 'raising factorial'
    rfac(a,n) = a*(a+1)*(a+2) ... (a+n-1)
    '''
    F    = 0*(a+b+z)
    an   = 1 + 0*a
    bn   = 1 + 0*b
    nfac = 1
    zn   = 1 + 0*z
    for i in xrange(1,order):
        zn   *= z
        an   *= (i+a)
        bn   *= (i+b)
        nfac *= i
        F    +=  an*zn/( bn*nfac )

def BoysF( x, n ):
    '''
    BoysF(x,n)  =KummerF(n+0.5,n+1.5,-x) / (2*n+1)
    '''
    return KummerF( n+0.5, n+1.5, -x ) / ( 2*n+1 )

def doubleFactorial(n):
    i = n
    out=1
    while i>0:
        out*=i
        i-=2
    return out

cs_hi = [
8.8622692545274984e-01, 
1.3293403881783550e+00, 
3.3233509704074469e+00, 
1.1631728394811573e+01, 
5.2342777714378713e+01, 
2.8788527514572553e+02, 
1.8712542074918215e+03, 
1.4034403753886485e+04
]

def BoysF_high( x, n ):
    '''
    approximation of Boys function for large argumetn (x>>1)
    http://shivupa.github.io/blog/efficient-evaluation-of-the-boys-function/
    BoysF_hi = ((2n-1)!! / 2^(n+1))  * sqrt(  pi / ( x^2n+1 ) )
    '''
    #return  doubleFactorial(2*n-1) * np.sqrt(np.pi)/ ( ( (2*x)**n ) * 2*np.sqrt(x) )
    return np.sqrt( 1/x**(2*n+1) ) * cs_hi[n]

def BoysF_approx( x, n, b=-0.6 ):
    #x = np.sqrt( x*x + x0**2 )
    #x = np.sqrt(  + x**2 )
    c  = cs_hi[n]
    #x0 = c**(2.0/(2*n+1))
    #x  = np.sqrt( x0**2 + x**2 )
    x0  = c**(2.0/(2*n+1))
    x   = x + x0*np.exp(b*x/(0.3*n+1))
    return np.sqrt( 1./( x**(2*n+1) ) ) * c

def BoysF_approx2( x, n, b=-0.6 ):
    #x = np.sqrt( x*x + x0**2 )
    #x = np.sqrt(  + x**2 )
    c  = cs_hi[n]
    #x0 = c**(2.0/(2*n+1))
    #x  = np.sqrt( x0**2 + x**2 )
    w  = c**2.0
    #w *= np.exp(b*x/(10*n+1))
    #w  *= 1./(1 + x/(10000*n+1)) 
    return np.sqrt( 1./( w + x**(2*n+1) ) ) * c

def BoysF_approx3( x, n, b=-0.6 ):
    c   = cs_hi[n]
    x0  = c**(2.0/(2*n+1))
    x0 *= 1./( 1.+x/(0.75*n+1) )
    x   = x + x0
    return np.sqrt( 1./( x**(2*n+1) ) ) * c

def BoysF_scipy(x, n, bPref=True):
    '''
    approximation of Boys function for large argumetn (x>>1)
    http://shivupa.github.io/blog/efficient-evaluation-of-the-boys-function/
    BoysF_hi = ((2n-1)!! / 2^(n+1))  * sqrt(  pi / ( x^2n+1 ) )
    '''
    out = hyp1f1( n+.5, n+1.5, -x )
    if bPref:
        out/=(2*n+1.)
    return out

'''
def invPowSeries(x, params ):
    x_    = x - params[0]
    invx  = 1./x_
    invxn = 0*x + 1. 
    f     = 0*x + params[1]
    for c in params[2:]:
        invxn *= invx
        f     += c*invxn
    return f

def makePade( f , x0, nP, nQ ):
    ctaylor = mpm.taylor( f,  x0,  nP+nQ )
    cP, cQ  = mpm.pade  ( ctaylor, nP, nQ )
    return cP,cQ 

def evalPade( x, cP, cQ ):
    return mpm.polyval( cP[::-1], x)/mpm.polyval( cQ[::-1], x)

def makePowSeriesBasis(x, nmax=4):
    ys = np.empty( (nmax,len(x)) )
    y  = 1 + 0*x
    for i in xrange(nmax):
        ys[i] = y
        y*=x
    return ys

def fitPowerSeries( x, y_ref, basis=None, nmax=4 ):
    if basis is None:
        basis = makePowSeriesBasis( x, nmax=nmax )
    #coefs = np.linalg.solve( basis, y_ref )
    print "basis.shape, y_ref.shape ", basis.shape, y_ref.shape
    coefs, res, rank, s = np.linalg.lstsq( basis.T, y_ref )
    print "residuals ", res
    return coefs
'''

def erf_p6(x):
    '''
    https://en.wikipedia.org/wiki/Error_function
    Abramowitz and Stegun 
    erf = 1 - 1/(1+a1*x +a2*x^2 +a3*x^3 +a4*x^4 +a5*x^5 +a6*x^6 )^16
    '''
    a1 = 0.0705230784
    a2 = 0.0422820123
    a3 = 0.0092705272
    a4 = 0.0001520143
    a5 = 0.0002765672
    a6 = 0.0000430638
    f  = 1/( 1 + x*( a1 + x*( a2 + x*( a3 + x*( a4 + x*(a5 + x*a6 )  ) ) ) ) )
    f *=f  # f^2
    f *=f  # f^4
    f *=f  # f^8
    f *=f  # f^16
    return 1 - f

# ================== Main

if __name__ == "__main__":
    xs = np.linspace(0,15,1000)
    import matplotlib.pyplot as plt
    
    '''
    nmax = 5
    colors = plt.cm.jet( np.linspace(0,1,nmax) ) 
    for n in xrange(nmax):
        ys_ref = BoysF_scipy(xs, n , bPref=False ); plt.plot(xs,ys_ref , c=colors[n],ls='-' )
        #ys_hi = BoysF_high  (xs, n               ); plt.plot(xs,ys_hi  , c=colors[n],ls=':' )
        #ys_hi = BoysF_approx( xs, n,              ); plt.plot(xs,ys_hi  , c=colors[n],ls=':' )
        ys_hi = BoysF_approx3( xs, n,              ); plt.plot(xs,ys_hi  , c=colors[n],ls=':' )
    '''
    
    bErrFunc = False
    if bErrFunc:
        ys_ref = erf_scipy( xs );   plt.plot(xs,ys_ref, label='erf_ref' )
        ys_p6  = erf_p6   ( xs );   plt.plot(xs,ys_p6,  label='erf_p6'  )
        
        plt.legend(); plt.grid()
        
        plt.figure()
        plt.plot(xs,abs(ys_p6-ys_ref), label='erf_ref' )
        plt.yscale('log')
        plt.legend(); plt.grid()
        plt.show()
        
        exit(0)
    
    

    
    bBoysX = False
    if bBoysX:
        ys_ref = BoysF_scipy( xs, 0, bPref=False ); plt.plot( xs, ys_ref, label='ref' )
        
        
        func = lambda x: BoysF_scipy(x,0)
        
        xstep=1.0
        sp_ys,sp_xs = splines.getSamples( func, 15, xmin=0.0, dx=xstep )
        plt.plot( sp_xs, sp_ys, 'x'  )
        ys_sp = splines.evalValTable( xs, sp_ys, xstep ) ;        
        plt.plot( xs, ys_sp , label="sp" )
        
        sp_ys,sp_xs = splines.getSamples( func, 15*2, xmin=0.0, dx=xstep*0.5 )
        plt.plot( sp_xs, sp_ys, 'x'  )
        ys_sp_2 = splines.evalValTable( xs, sp_ys, xstep*0.5 ) ;
        plt.plot( xs, ys_sp_2 , label="sp" )
        
        sp_ys,sp_dys,sp_xs = splines.getSamplesDeriv( func, 15, xmin=0.0, dx=xstep )
        ys_sp_ = splines.evalValTableDeriv( xs, sp_ys, sp_dys, xstep ) ;
        plt.plot( xs, ys_sp_ , label="spDeriv" )
        
        sp_ys,sp_dys,sp_xs = splines.getSamplesDeriv( func, 15, xmin=0.0, dx=xstep, fddx=0.025 )
        ys_sp__ = splines.evalValTableDeriv( xs, sp_ys, sp_dys, xstep ) ;
        plt.plot( xs, ys_sp__ , label="hi" )
        
        ys_hi = BoysF_high( xs, 0 )
        
        
        plt.legend()
        plt.grid()
        
        plt.figure()
        plt.plot(xs, abs(ys_sp -ys_ref ), label="sp" )
        plt.plot(xs, abs(ys_sp_2 -ys_ref ), label="sp2" )
        plt.plot(xs, abs(ys_sp_-ys_ref ), label="Deriv")
        plt.plot(xs, abs(ys_sp__-ys_ref), label="Deriv")
        plt.plot(xs, abs(ys_hi-ys_ref), label="hi")
        plt.yscale('log')
        plt.legend()
    
    bBoysX2 = True
    if bBoysX2:
        xs  = np.linspace(0.0,3.8,500)
        x2s = xs**2
        ys_ref = BoysF_scipy( x2s, 0, bPref=False ); plt.plot( xs, ys_ref, label='ref' )

        xstep=0.25
        
        func = lambda x: BoysF_scipy(x**2,0)
        
        sp_ys,sp_xs = splines.getSamples( func, 15, xmin=0.0, dx=xstep )
        plt.plot( sp_xs, sp_ys, 'x'  )
        ys_sp = splines.evalValTable( xs, sp_ys, xstep ) ;        
        plt.plot( xs, ys_sp , label="sp" )
        
        
        sp_ys,sp_xs = splines.getSamples( func, 15*2, xmin=0.0, dx=xstep*0.5 )
        plt.plot( sp_xs, sp_ys, 'x'  )
        ys_sp_2 = splines.evalValTable( xs, sp_ys, xstep*0.5 ) ;
        plt.plot( xs, ys_sp_2 , label="sp2x" )
        
        sp_ys,sp_dys,sp_xs = splines.getSamplesDeriv( func, 15, xmin=0.0, dx=xstep )
        ys_sp_ = splines.evalValTableDeriv( xs, sp_ys, sp_dys, xstep ) ;
        plt.plot( xs, ys_sp_ , label="spDeriv" )
        
        sp_ys,sp_dys,sp_xs = splines.getSamplesDeriv( func, 15, xmin=0.0, dx=xstep, fddx=0.025 )
        ys_sp__ = splines.evalValTableDeriv( xs, sp_ys, sp_dys, xstep ) ;
        plt.plot( xs, ys_sp__ , label="sp_0.025" )
        
        ys_hi = BoysF_high( x2s, 0 );   plt.plot( xs, ys_hi , label="hi" )
        
        plt.ylim(0.0,1.2)
        plt.legend()
        plt.grid()
        
        plt.figure()
        plt.plot(xs, abs(ys_sp -ys_ref ), label="sp" )
        plt.plot(xs, abs(ys_sp_2 -ys_ref ), label="sp2" )
        plt.plot(xs, abs(ys_sp_-ys_ref ), label="Deriv")
        plt.plot(xs, abs(ys_sp__-ys_ref), label="Deriv")
        plt.plot(xs, abs(ys_hi-ys_ref), label="hi")
        plt.yscale('log')
        plt.legend()
    
    '''
    n = 10
    ys_invpow = invPowSeries(xs, [-1.0,0.0,0.0,5.0] ); plt.plot(xs,ys_invpow)

    ys_ref = BoysF_scipy(xs, n , bPref=False ); plt.plot(xs,ys_ref,'k')
    zs     = np.sqrt(1/(xs**(2*n+1)))
    
    pref = zs/ys_ref;    plt.plot(xs,pref,':')
    '''

    
    '''
    x0 = 30.0
    npow = 8
    c_hi = np.empty(npow)
    for n in range(npow):
        y = BoysF_scipy( x0, n, bPref=False )
        z = np.sqrt(1/(x0**(2*n+1)))
        print " n, y/z, y  ",  "%i %20.20f %20.20f  "  %(n, y/z, y)
        c_hi[n]  =  y/z
        #print " n, y/z, y  ",  "%i %20.20f %20.20f  "  %(n, y/z, y)
    print c_hi
    '''
    
    '''
    def BoysF_0(x):
        return BoysF_scipy(xs, 0, bPref=False )[0]
    
    cP,cQ = makePade( BoysF_0, 1., 2, 2 )
    print cP,cQ
    '''
    
    '''
    ifit_min   = 10
    xs_fit     = xs[ifit_min:] 
    basis = makePowSeriesBasis( xs_fit, nmax=4 )
    coefs = fitPowerSeries( xs_fit+1., ys_ref[ifit_min:], basis=basis )
    ys_fit = np.dot(basis.T,coefs); plt.plot(xs_fit,ys_fit )
    print "coefs ", coefs
    '''
    
    plt.grid()
    plt.ylim(0,1.5)
    plt.show()


