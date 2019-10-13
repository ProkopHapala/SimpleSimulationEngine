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

def PolyDerivMatrix(  ndeg, nderiv ):
    '''
    P = Sum_k{  a_k x^k }
    '''
    M   = np.empty( (nderiv,ndeg), dtype=np.int)
    M[0,:] = 1
    ints = np.arange(0,ndeg,dtype=np.int)
    M[1,:] = ints
    for i in xrange(2,nderiv ):
        M[i,:] = M[i-1,:] * (ints-(i-1))
    return M

def HermiteSplineN( derivs_0, derivs_1 ):
    '''
    construct polynominal given number of derivatives at 0 and 1
    '''
    n0   = len( derivs_0 )
    n1   = len( derivs_1 )
    ndeg = n0 + n1 
    D = PolyDerivMatrix( ndeg, max(n0,n1) )
    A = np.zeros( (ndeg,ndeg) )
    b = np.array( [d for d in derivs_0 ]+[ d for d in derivs_1 ] )
    for i in xrange(n0):
        A[i,i] = D[i,i]
    for i in xrange(n1):
        A[n0+i,:] = D[i,:]
    #print "A\n",A
    #print "b\n",b
    coefs = np.linalg.solve(A,b)
    return coefs


def HermiteSplineMatrix( nder ):
    '''
    construct polynominal given number of derivatives at 0 and 1
    '''
    ndeg = nder*2 
    D = PolyDerivMatrix( ndeg, nder )
    print "------D\n",D
    A = np.zeros( (ndeg,ndeg) )
    for i in xrange(nder):
        A[i,i]      = D[i,i]
    A[nder:,:] = D[:,:]
    #print "A\n",A
    #print "b\n",b
    #coefs = np.linalg.solve(A,b)
    return A

"""
HermiteSpline_C1    [y,dy]
// y0 dy0  y1 dy1
[[ 1.  0.  0.  0.]    // 1
 [ 0.  1.  0.  0.]    // x
 [-3. -2.  3. -1.]    // x^2
 [ 2.  1. -2.  1.]]   // x^3

HermiteSpline_C2    [y,dy,ddy]
//  y0   y'0  y''0    y1   y'1   y''1
[[  1.    0.    0.    0.    0.    0. ]   // 1
 [  0.    1.    0.    0.    0.    0. ]   // x
 [  0.    0.    0.5   0.    0.    0. ]   // x^2
 [-10.   -6.   -1.5  10.   -4.    0.5]   // x^3
 [ 15.    8.    1.5 -15.    7.   -1. ]   // x^4
 [ -6.   -3.   -0.5   6.   -3.    0.5]]  // x^5


[[ 1.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00 0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00]
 [ 0.00000000e+00  1.00000000e+00  0.00000000e+00  0.00000000e+00 0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00]
 [ 0.00000000e+00  0.00000000e+00  5.00000000e-01  0.00000000e+00 0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00]
 [ 0.00000000e+00  0.00000000e+00  0.00000000e+00  3.33333333e-02 0.00000000e+00  0.00000000e+00  0.00000000e+00  0.00000000e+00]
 [-2.80000000e+01 -1.70000000e+01 -4.50000000e+00 -1.26666667e-01  2.80000000e+01 -1.10000000e+01  1.50833333e+00 -8.33333333e-03]
 [ 6.30000000e+01  3.60000000e+01  8.50000000e+00  1.80000000e-01 -6.30000000e+01  2.70000000e+01 -4.02500000e+00  2.50000000e-02]
 [-4.90000000e+01 -2.70000000e+01 -6.00000000e+00 -1.13333333e-01  4.90000000e+01 -2.20000000e+01  3.52500000e+00 -2.50000000e-02]
 [ 1.30000000e+01  7.00000000e+00  1.50000000e+00  2.66666667e-02 -1.30000000e+01  6.00000000e+00 -1.00833333e+00  8.33333333e-03]]

[[  1.,   0.,  0.0,   0.,             0.,   0.,   0.,           0.,            ]
 [  0.,   1.,  0.0,   0.,             0.,   0.,   0.,           0.,            ]
 [  0.,   0.,  0.5,   0.,             0.,   0.,   0.,           0.,            ]
 [  0.,   0.,  0.0,   0.0333333333,   0.,   0.,   0.,           0.,            ]
 [-28., -17., -4.5,  -0.126666667,   28., -11.,   1.50833333,  -0.00833333333, ]
 [ 63.,  36.,  8.5,   0.18,         -63.,  27.,  -4.025,        0.0250000000,  ]
 [-49., -27., -6.0,  -0.113333333,   49., -22.,   3.525,       -0.0250000000,  ]
 [ 13.,   7.,  1.5,   0.0266666667, -13.,   6.,  -1.00833333,   0.00833333333, ]]

# https://www.johndcook.com/rational_approximation.html
//  y0    y'0  y''0   y'''0      y1    y'1   y''1         y'''1
[[  1.,   0.,  0.0,   0.0,        0.,   0.,   0.,           0,.   ]    // 1
 [  0.,   1.,  0.0,   0.0,        0.,   0.,   0.,           0.,   ]    // x
 [  0.,   0.,  0.5,   0.0,        0.,   0,.   0.,           0.,   ]    // x^2
 [  0.,   0.,  0.0,   1./30,      0.,   0.,   0.,           0.,   ]    // x^3
 [-28., -17., -4.5,  -19./150,   28., -11.,   181./120,   1./120, ]    // x^4
 [ 63.,  36.,  8.5,   0.18,     -63.,  27.,  -161./40,    1./40,  ]    // x^5
 [-49., -27., -6.0,  -17./150,   49., -22.,    61./40,   -1./40,  ]    // x^6
 [ 13.,   7.,  1.5,   2./75 ,   -13.,   6.,  -121./120,   1./120, ]]   // x^7


"""


'''
Cos Coefs [ -pi/2, pi/2]:
degree  4  [                                                                3.71276128e-02 -4.96257107e-01  9.99555541e-01 ]
degree  6  [                                               -1.27676413e-03  4.15061557e-02 -4.99927144e-01  9.99994936e-01 ]
degree  8  [                                2.32115274e-05 -1.38563556e-03  4.16639592e-02 -4.99999219e-01  9.99999964e-01 ]
degree  10 [                -2.60928931e-07 2.47625965e-05 -1.38884110e-03  4.16666403e-02 -4.99999995e-01  1.00000000e+00 ]
degree  12 [ 1.99300103e-09 -2.75262798e-07 2.48010828e-05 -1.38888847e-03  4.16666665e-02 -5.00000000e-01  1.00000000e+00 ]

Taylor Cos [ 1./479001600,  -1./3628800,   1./40320,  -1./720,  1./24,  -1./2, 1. ] 
Taylor Sin [ 1./6227020800, -1./39916800,  1./362880, -1./5040, 1./120, -1./6, 1. ] 
'''

'''

1      a0
x      -a0 + a1 + b0 - b1 + b10 + b2 - b3 + b4 - b5 + b6 - b7 + b8 - b9
x**2      -a1 + a2 + b1 - 10*b10 - 2*b2 + 3*b3 - 4*b4 + 5*b5 - 6*b6 + 7*b7 - 8*b8 + 9*b9
x**3      -a2 + a3 + 45*b10 + b2 - 3*b3 + 6*b4 - 10*b5 + 15*b6 - 21*b7 + 28*b8 - 36*b9
x**4      -a3 + a4 - 120*b10 + b3 - 4*b4 + 10*b5 - 20*b6 + 35*b7 - 56*b8 + 84*b9
x**5      -a4 + a5 + 210*b10 + b4 - 5*b5 + 15*b6 - 35*b7 + 70*b8 - 126*b9
x**6      -a5 + a6 - 252*b10 + b5 - 6*b6 + 21*b7 - 56*b8 + 126*b9
x**7      -a6 + a7 + 210*b10 + b6 - 7*b7 + 28*b8 - 84*b9
x**8      -a7 + a8 - 120*b10 + b7 - 8*b8 + 36*b9
x**9      -a8 + a9 + 45*b10 + b8 - 9*b9
x**10      a10 - a9 - 10*b10 + b9


1          a0
x         -a0  + a1  + b0 - b1 +   b2 -   b3 +   b4 -   b5  +    b6 -    b7 +    b8 -     b9 +     b10
x**2      -a1  + a2       + b1 - 2*b2 + 3*b3 - 4*b4 + 5*b5  - 6 *b6 + 7* b7 - 8* b8 + 9  *b9 - 10 *b10
x**3      -a2  + a3            +   b2 - 3*b3 + 6*b4 - 10*b5 + 15*b6 - 21*b7 + 28*b8 - 36 *b9 + 45 *b10
x**4      -a3  + a4                   +   b3 - 4*b4 + 10*b5 - 20*b6 + 35*b7 - 56*b8 + 84 *b9 - 120*b10
x**5      -a4  + a5                          +   b4 - 5*b5  + 15*b6 - 35*b7 + 70*b8 - 126*b9 + 210*b10
x**6      -a5  + a6                               +     b5  - 6 *b6 + 21*b7 - 56*b8 + 126*b9 - 252*b10
x**7      -a6  + a7                                         +    b6 - 7 *b7 + 28*b8 - 84 *b9 + 210*b10
x**8      -a7  + a8                                             +        b7 - 8* b8 + 36 *b9 - 120*b10
x**9      -a8  + a9                                                    +         b8 - 9  *b9 + 45 *b10
x**10      a10 - a9                                                           +           b9 - 10 *b10


[    
    [+1, -1, +1, -1, +1, - 1, + 1, - 1, + 1,   -1,  +1],
    [ 0, +1, -2, +3, -4, + 5, - 6, + 7, - 8,   +9, -10],
    [ 0,  0, +1, -3, +6, -10, +15, -21, +28,  -36, +45],
    [ 0,  0,  0, +1, -4, +10, -20, +35, -56,  +84,-120],
    [ 0,  0,  0,  0, +1, - 5, +15, -35, +70, -126,+210],
    [ 0,  0,  0,  0,  0, + 1, - 6, +21, -56, +126,-252],
    [ 0,  0,  0,  0,  0,   0, + 1, - 7, +28,  -84,+210],
    [ 0,  0,  0,  0,  0,   0,   0, + 1, - 8,  +36,-120],
    [ 0,  0,  0,  0,  0,   0,   0,   0, + 1,  - 9,+ 45],
    [ 0,  0,  0,  0,  0,   0,   0,   0,   0,   +1,- 10],
]

'''


def getHermitePoly( xs, H, yds ):
    coefs = np.dot( H, yds ) #  ;print "cC ", cC
    return np.polyval( coefs[::-1], xs )


if __name__ == "__main__":
    np.set_printoptions(precision=16, linewidth=200 )

    A = HermiteSplineMatrix( 4 )
    print A
    H4_ = np.linalg.inv(A)
    print "inv(A)\n", H4_

    H1 = np.array([
    #      y0   y1 
        [  1.,  0. ],   # 1
        [ -1.,  1. ],   # x
    ])

    H2 = np.array([
    #    y0 dy0  y1 dy1
        [ 1.,  0.,  0.,  0.],   # 1
        [ 0.,  1.,  0.,  0.],   # x
        [-3., -2.,  3., -1.],   # x^2
        [ 2.,  1., -2.,  1.],   # x^3
    ])

    H3 = np.array([
        #  y0   y'0  y''0    y1   y'1   y''1
        [  1.,    0.,    0.,    0.,    0.,    0.  ], # 1
        [  0.,    1.,    0.,    0.,    0.,    0.  ], # x
        [  0.,    0.,    0.5,   0.,    0.,    0.  ], # x^2
        [-10.,   -6.,   -1.5,  10.,   -4.,    0.5 ], # x^3
        [ 15.,    8.,    1.5, -15.,    7.,   -1.  ], # x^4
        [ -6.,   -3.,   -0.5,   6.,   -3.,    0.5 ]  # x^5
    ])

    yds = np.array( [  0.7,   -1.3, 0.2,   -0.116,        0.65,  -0.335,   -0.456,   +1.65,     ] )
    coefs = np.dot( H4_, yds )

    print "coefs ", coefs

    D = PolyDerivMatrix( 8, 4 )
    print "D\n", D
    for i,l in enumerate(D):
        print "deriv[%i]" %i, np.dot( l, coefs ), "   ref ", yds[4+i] 
    for i,l in enumerate(D):
        csi= l * coefs
        print "csi ", csi
    xs_ends = np.array( [0.,1.] )
    for i,l in enumerate(D):
        csi= l * coefs
        csi = csi[i:] 
        ys_end = np.polyval( csi[::-1], xs_ends )
        print "derivs[%i]" %i, ys_end 

    #exit()




    Hs = [ H1, H2, H3, H4_]

    dsc = 1.57079632679

    dscs = np.array([ 1.,dsc,dsc*dsc, dsc*dsc*dsc,        1.,dsc,dsc*dsc,dsc*dsc*dsc, ])

    fnames = ['sin','cos']
    funcs = [ np.cos, np.sin]
    yds=np.array([
        [  1.,   0., -1.,     0.,        0.,  -1.,   0.,   0./6,     ],
        [  0.,   1.,  0.,   -0./6,        1.,   0.,  -1.,    0.,      ]
    ])

    for i in range(len(fnames)):
        yds[i]*=dscs

    import matplotlib.pyplot as plt

    xs = np.linspace( 0, np.pi/2, 100 )
    for i,fname in enumerate(fnames):
        ys_ref = funcs[i](xs)
        plt.figure()
        #plt.plot(xs, ys_ref, '-', label=fname )
        for j,H in enumerate(Hs):
            inds = range(j+1)+range(4,4+j+1)
            print "j, inds ", j, inds
            yds_j = yds[i, inds  ]
            ys = getHermitePoly( xs/dsc, H, yds_j )
            #plt.plot(xs, ys, ':', label=fname+("deg.%i" %(j+1))  )
            yerr = ys - ys_ref
            plt.plot(xs, abs(yerr), '-', label=fname+("_deg.%i" %(j+1))  )
        plt.yscale('log')
        plt.ylim(1e-9,1)
        plt.legend()
        plt.grid()

    plt.show()

    '''
    import sympy as sy

    x, a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,  b0,b1,b2,b3,b4,b5,b6,b7,b8,b9,b10 = sy.symbols('x a0 a1 a2 a3 a4 a5 a6 a7 a8 a9 a10  b0 b1 b2 b3 b4 b5 b6 b7 b8 b9 b10')

    y=x-1

    f1 = a0 +x*( a1 +x*( a2 +x*( a3 +x*(a4 +x*( a5 +x*( a6 +x*( a7 +x*( a8 +x*( a9 +x*( a10 ) ) ) ) ) ) ) ) ) )
    f2 = b0 +y*( b1 +y*( b2 +y*( b3 +y*(b4 +y*( b5 +y*( b6 +y*( b7 +y*( b8 +y*( b9 +y*( b10 ) ) ) ) ) ) ) ) ) )

    f = (1-x)*f1 + x*f2

    f = sy.expand( f )
    f.simplify()
    dct = sy.collect( f,x, evaluate=False)
    ndegs = 11
    for i in range(ndegs):
        term = x**i
        print term,"    ", dct[term] 
    '''

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