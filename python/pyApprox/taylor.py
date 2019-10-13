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


def getHermitePoly( xs, H, yds ):
    coefs = np.dot( H, yds ) #  ;print "cC ", cC
    return np.polyval( coefs[::-1], xs )


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    np.set_printoptions(precision=16, linewidth=200 )

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