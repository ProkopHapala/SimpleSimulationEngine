# https://stackoverflow.com/questions/18662261/fastest-implementation-of-sine-cosine-and-square-root-in-c-doesnt-need-to-b
# https://en.wikipedia.org/wiki/Bhaskara_I's_sine_approximation_formula

# https://docs.scipy.org/doc/numpy-1.13.0/reference/generated/numpy.polynomial.chebyshev.Chebyshev.fit.html

import numpy as np

def polyAprox( func=None, xs=None, ys=None, xmin=0.0, xmax=1.0,  nps=100, ndeg=6 ):
    if xs is None:
        xs = np.linspace( xmin, xmax, nps )
    if ys is None:
        ys = func(xs)
    return np.polyfit( xs, ys, ndeg )

def PolyDerivMatrix( ndeg, nderiv ):
    '''
    P = Sum_k{  a_k x^k }
    '''
    D   = np.empty( (nderiv+1,ndeg), dtype=np.int)
    D[0,:] = 1
    ints = np.arange(0,ndeg,dtype=np.int)
    D[1,:] = ints
    for i in xrange(2,nderiv ):
        D[i,:] = D[i-1,:] * (ints-(i-1))
    return D

def HermiteSplineMatrix( nder ):
    '''
    construct polynominal given number of derivatives at 0 and 1
    '''
    ndeg = nder*2 
    print "HermiteSplineMatrix nder, ndeg ", nder, ndeg
    D = PolyDerivMatrix( ndeg, nder )
    A = np.zeros( (ndeg,ndeg) )
    for i in xrange(nder):
        A[i,i] = D[i,i]
    A[nder:,:] = D[:-1,:]
    return A

def numDeriv( xs, ys ):
    dys = (ys[1:]-ys[:-1])/(xs[1:]-xs[:-1])
    xs_ = (xs[1:]+xs[:-1])*0.5
    return dys, xs_

def getNthDeriv_poly( xs, i, D, coefs ):
    csi = D[i] * coefs
    csi = csi[i:] 
    return np.polyval( csi[::-1], xs )

def printHornerPolynom(coefs):
    print " %.16g" %coefs[0],
    for c in coefs[1:]:
        print "+x*(%.16g" %c,
    for c in coefs[1:]:
        print ")",
    print

def printHornerPolynom_EvenOdd(coefs):
    print "even: ", coefs[0],
    for c in coefs[2::2]:
        print "+xx*(%.16g" %c,
    for c in coefs[2::2]: print ")",
    print "\nodd: x*(", coefs[1],
    for c in coefs[3::2]:
        print "+xx*(%.16g" %c,
    for c in coefs[1::2]: print ")",
    #for c in coefs[1:]:
    #    print ")",
    print


def getHermitePoly( xs, H, yds ):
    coefs = np.dot( H, yds ) #  ;print "cC ", cC
    return np.polyval( coefs[::-1], xs )

def makePolynomBasis( xs, orders ):
    Bas = np.empty( ( len(orders), len(xs) ) )
    for i,k in enumerate( orders ):
        Bas[i,:] = xs**k
    return Bas

def polyFit( xs, y_ref, orders ):
    Bas    = makePolynomBasis( xs, orders )
    #print Bas.shape, xs.shape
    coefs_ = np.linalg.lstsq( Bas.T, y_ref )[0]
    coefs  = np.zeros( orders[-1]+1 )
    #print len(coefs_),len(coefs), orders
    for i,k in enumerate(orders):
        #print i,k
        coefs[k] = coefs_[i]
    return coefs

def polyFitFunc( xs, y_ref, coef0, orders, nps=100 ):
    coef0 = np.array( coef0 )
    y0s   = np.polyval( coef0[::-1], xs )
    coefs = polyFit( xs, y_ref-y0s, orders )
    coefs[:len(coef0)] = coef0[:]
    ys = np.polyval( coefs[::-1], xs )
    return coefs, ys

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    np.set_printoptions(precision=16, linewidth=200 )

    #phis   = np.linspace( -np.pi/3, np.pi/3, 1000 )
    #phis   = np.linspace( 0, np.pi/3, 1000 )
    phis   = np.linspace( np.pi/3, 2*np.pi/3, 1000 )
    xs     = np.cos(phis)
    ys     = np.sin(phis)
    y_ref  = phis/(np.pi/3) - 1.5
    print "y_ref[-1]-y_ref[0] ", y_ref[-1]-y_ref[0] 
    #args   = ys/xs  
    #args   = -xs/ys
    args   = -0.86602540378*xs/ys  

    plt.plot(phis, y_ref, '-',label="y_ref"    )
    plt.plot(phis, args , '-',label="args" )

    #coef0s = [0,1]
    coef0s = [0]
    #for maxOrder in [6,8,10,12,14,16,18,20,22]:
    for maxOrder in [6,8,10,12,14]:
        #coefs, ys = polyFitFunc( xs, y_ref, coef0s, range(5,maxOrder,2) )
        coefs, ys = polyFitFunc( args, y_ref, coef0s, range(1,maxOrder,2) )
        #coefs, ys = polyFitFunc( args, y_ref, coef0s, range(1,maxOrder) )
        #print maxOrder,":  "; printHornerPolynom(coefs)
        print maxOrder,":  "; printHornerPolynom_EvenOdd(coefs)
        y_err   = ys - y_ref
        #plt.plot(phis, ys       , ':',label=('y_%i' %maxOrder ) )
        plt.plot(phis,abs(y_err), '-',label=('err_%i' %maxOrder ) )

    plt.yscale('log')
    plt.legend()
    plt.grid()
    plt.title( "Polynominal Approx tan(x)" )
    plt.show()
    exit()



    '''
    # ========= PolyFit :  ata2(y,x)

    def atan2_pre(y,x, y_ref):
        abs_y = abs(y) + 1e-14;
        mask_x = x<0
        mask_y = y<0
        a         = (( x - abs_y ) / ( abs_y + x ))
        a[mask_x] = (( x + abs_y ) / ( abs_y - x ))[mask_x]
        #angle         = x*0 + 0.78539816339
        #angle[mask_x] =       2.35619449019
        #a[mask_y] *= -1
        y_ref[mask_y] *= -1
        y_ref[x<0 ] -=  2.35619449019
        y_ref[x>=0] -=  0.78539816339
        #y_ref[mask_y] *= -1
        return a, y_ref

    def atan2(y,x):
        a, angle = atan2_pre(y,x, x*0); angle*=-1
        aa     = a * a
        #a[y<0] *= -1
        angle += a * ( -1 + aa*( 0.331768825725 + aa*( -0.184940152398 + aa*( 0.091121250024 -0.0233480867489*aa ) ) ) )
        #angle[y<0] *= -1
        return angle 

    xmax   = np.pi*0.5
    phis   = np.linspace( -xmax, xmax, 1000 )
    xs = np.cos(phis)
    ys = np.sin(phis)
    y_ref  = np.arctan2(ys,xs)
    y      = atan2(ys,xs)
    args,  y_ref_mod = atan2_pre( ys, xs, y_ref )
    #plt.plot(phis,args*-1     , '-',label="a" )
    #plt.plot(phis,y_ref_mod, '-',label="y_ref_mod" )
    coef0s =  [0,-1.0]
    for maxOrder in [8,10,12,14,16]:
        #coefs, ys = polyFitFunc( xs, y_ref, coef0s, range(5,maxOrder,2) )
        coefs, ys = polyFitFunc( args, y_ref_mod, coef0s, range(3,maxOrder,2) )
        printHornerPolynom(coefs)
        y_err   = ys - y_ref_mod
        #ctg_err = 1/ys - 1/y_ref
        plt.plot(phis,abs(y_err), '-',label=('tan_%i' %maxOrder ) )
        #plt.plot(args,abs(ctg_err), ':',label=('cotg_%i' %maxOrder) )
    #plt.plot(phis, y_ref, '-',label="atan2_ref"    )
    #plt.plot(phis, y    , '-',label="atan2_approx" )
    plt.yscale('log')
    plt.legend()
    plt.grid()
    plt.title( "Polynominal Approx tan(x)" )
    plt.show()
    exit()
    '''

    '''
    # ========= PolyFit :  tan(x)
    xmax = np.pi/4
    xs     = np.linspace( -xmax, xmax, 100 )
    y_ref  = np.tan(xs)
    coef0s =  [0,1.0,0.0,1./3]
    for maxOrder in [8,10,12,14]:
        coefs, ys = polyFitFunc( xs, y_ref, coef0s, range(5,maxOrder,2) )
        printHornerPolynom(coefs)
        y_err   = ys - y_ref
        ctg_err = 1/ys - 1/y_ref
        plt.plot(xs,abs(y_err), '-',label=('tan_%i' %maxOrder ) )
        plt.plot(xs,abs(ctg_err), ':',label=('cotg_%i' %maxOrder) )
    
    plt.yscale('log')
    plt.legend()
    plt.grid()
    plt.title( "Polynominal Approx tan(x)" )
    plt.show()
    exit()
    '''

    '''
    # ========= PolyFit :   acos
    # acosf implementation:
    # https://github.com/bminor/glibc/blob/master/sysdeps/ieee754/flt-32/e_acosf.c
    #xs     = np.linspace( -1.0, 1.0, 100 )
    #xs__    = np.linspace( 0, 1.0, 100 )**2 - 1
    xs      = np.linspace( -0.75, 0.75, 100 )
    xs_     = np.linspace(  0.75, 1.000, 100 )
    y_ref   = np.arcsin(xs)
    y_ref_  = np.arcsin(xs_)
    coef0s =  [0,1.0]
    k = (y_ref_[-1]-y_ref_[0])/0.25
    coef0s_ =  [ y_ref_[0]-k*xs_[0] , k ]
    for maxOrder in [8,10,12,14]:
        coefs, ys = polyFitFunc( xs, y_ref, coef0s, range(3,maxOrder,2) )
        printHornerPolynom(coefs)
        y_err = ys - y_ref

        #plt.plot(xs,y_ref, '-',label='ref'   )
        #plt.plot(xs,ys   , ':',label='approx')

        plt.plot(xs,abs(y_err), '-',label='error'); plt.yscale('log')

        #coefs2, ys_ = polyFitFunc( xs_, y_ref_, coef0s_, range(3,6) )
        ys_ = np.polyval( coef0s_[::-1], xs_ )
        # sin(y) = x
        # y : y - sin(y)/sin'(y) = y - tan(y)
        #printHornerPolynom(coefs2)
        y_err_ = ys_ - y_ref_

        #plt.plot(xs_,y_ref_   , '-',label='ref')
        #plt.plot(xs_,ys_   , ':',label='approx')

        plt.plot(xs_,abs(y_err_), ':',label='error'); plt.yscale('log')

    plt.legend()
    plt.grid()
    plt.show()
    exit()
    '''



    '''
    # ========= PolyFit :   sin, cos
    dsc    = np.pi/4
    xs     = np.linspace( -1.0, 1.0, 100 )
    #y_ref  = np.cos(xs*dsc)
    y_ref  = np.sin(xs*dsc)
    #coef0s =  [1,0,-0.5*dsc*dsc,0]
    #coef0s =  [0,dsc,0,(-1./6)*dsc*dsc*dsc]
    coef0s =  [0,dsc,0 ]
    #coef0s =  [0,1,0,0]
    for maxOrder in [8,10,12]:
        coefs, ys = polyFitFunc( xs, y_ref, coef0s, range(3,maxOrder,2) )
        #ys = np.polyval( coef0s[::-1], xs )
        #print "coefs ", coefs
        #print " ", coefs
        printHornerPolynom(coefs)
        y_err = ys - y_ref

        #plt.plot(xs,y_ref, '-',label='ref'   )
        #plt.plot(xs,ys   , ':',label='approx')

        plt.plot(xs,abs(y_err), ':',label='error')
        plt.yscale('log')

    plt.legend()
    plt.grid()
    plt.show()

    exit()
    '''



    '''
    Hermite spline nth-degree Spline with n-de

    TaylorCos = [ 1., -1./2, 1./24,  -1./720,  1./40320,  -1./3628800,  1./479001600  ]
    TaylorSin = [ 1., -1./6, 1./120, -1./5040, 1./362880, -1./39916800, 1./6227020800 ]

    dsc=np.pi/4
    print "TaylorCos", printHornerPolynom( [ c*(dsc**(i*2  )) for i,c in enumerate(TaylorCos)  ] )
    print "TaylorSin", printHornerPolynom( [ c*(dsc**(i*2+1)) for i,c in enumerate(TaylorSin)  ] )
    exit()

    nderiv = 5

    A = HermiteSplineMatrix( nderiv ) #;print A
    H = np.linalg.inv(A)       #;print "inv(A)\n", H4_

    #yds = np.array( [  0.7,   +0.1316, -0.2669,   -0.116,        +0.65,  -0.1335,   -0.1456,   -0.0659,     ] )

    s = np.sqrt(0.5)

    #yds = np.array( [  1., 0., -1., 0., 1.,        0., -1., 0., +1, 0.     ] )
    #yds = np.array( [  1., 0., -1., 0., 1.,        s, -s, -s, s, s     ] )    # cos
    yds = np.array(  [  0., 1.,  0.,-1., 0.,        s,  s, -s,-s,+s     ] )    # sin
    
    ref_func = np.sin

    dsc  = np.pi/4
    dscs = np.array([ dsc**i for i in range(nderiv) ]*2)
    yds*=dscs

    coefs = np.dot( H, yds )

    print "coefs ", coefs

    D = PolyDerivMatrix( nderiv*2, nderiv )
    print "D\n", D

    # --- check plot with num derivs
    xs = np.linspace( 0, 1.0, 100 )
    ys_num = None
    for i in range(0, nderiv-1 ):
        ys = getNthDeriv_poly(xs, i, D, coefs )
        if ys_num is not None:
            ys_num,xs_num = numDeriv( xs_num,ys_num )
        else:
            #ys_num = ys.copy()
            ys_num = ref_func(xs*dsc)
            xs_num = xs.copy()
        plt.plot( xs, ys, '-', label='poly %i' %i )
        plt.plot( xs_num, ys_num, ':', label='num %i' %i )
        plt.plot( [0,1],[yds[i],yds[nderiv+i]],'+' )
    
    #plt.ylim(-2,2)
    plt.legend()
    plt.grid()
    plt.show()
    #exit()

    Hs = [  np.linalg.inv( HermiteSplineMatrix( ndv ) ) for ndv in range(1,nderiv+1) ]

    fnames = ['cos','sin']
    funcs  = [ np.cos, np.sin]
    yds=np.array([
        #[  1., 0.,-1., 0., 1.,       0.,-1., 0.,+1., 0.,   ],
        #[  0., 1., 0.,-1., 0.,       1., 0.,-1., 0., 1.,   ],
        [  1., 0., -1., 0., 1.,        s, -s, -s, s, s     ],    # cos
        [  0., 1.,  0.,-1., 0.,        s,  s, -s,-s,+s     ],    # sin
    ])

    for i in range(len(fnames)):
        yds[i]*=dscs

    xs = np.linspace( 0, dsc, 100 )
    for i,fname in enumerate(fnames):
        ys_ref = funcs[i](xs)
        plt.figure()
        #plt.plot(xs, ys_ref, '-', label=fname )
        for j,H in enumerate(Hs):
            inds = range(j+1)+range(nderiv,nderiv+j+1)
            #print "j, inds ", j, inds
            yds_j = yds[i, inds  ]
            #ys = getHermitePoly( xs/dsc, H, yds_j )
            coefs = np.dot( H, yds_j ) #;print "coefs[%s,%i] " %(fname,j), coefs
            printHornerPolynom(coefs)
            ys = np.polyval( coefs[::-1], xs/dsc )
            #plt.plot(xs, ys, ':', label=fname+("deg.%i" %(j+1))  )
            yerr = ys - ys_ref
            plt.plot(xs, abs(yerr), '-', label=fname+("_deg.%i" %(j+1))  )
        plt.yscale('log')
        plt.ylim(1e-12,1)
        plt.legend()
        plt.grid()

    plt.show()
    '''

    '''
    # ======== Cos

    import matplotlib.pyplot as plt

    xs = np.linspace( -np.pi/2, np.pi/2, 100 )
    ys_ref = np.cos(xs)

    ndegs = [4,6,8,10,12]

    for ndeg in ndegs:
        #coefs = polyAprox( xs=xs, ys= ndeg=6 )
        coefs   = np.polyfit( xs, ys_ref, ndeg ); print "degree ", ndeg, coefs
        ys_poly = np.polyval( coefs, xs )
        ys_err  = ys_poly-ys_ref

        #plt.figure()
        #plt.plot(xs,ys_ref , '-k')
        #plt.plot(xs,ys_poly, ':r')

        #plt.figure()
        plt.plot(xs,abs(ys_err) , '-', label="min deg.%i" %ndeg )
    
    c_taylor = [ 1./479001600,  -1./3628800,   1./40320,  -1./720,  1./24,  -1./2, 1. ] 

    x2s=xs**2
    ntay = len(c_taylor)
    for ndeg in range( ntay ):
        coefs = c_taylor[ntay-ndeg:]
        ys_poly = np.polyval( coefs, x2s )
        ys_err  = ys_poly-ys_ref
        plt.plot(xs,abs(ys_err) , '-', label="tay.%i" %ndeg )

    plt.yscale('log')
    plt.legend()

    plt.show()
    '''