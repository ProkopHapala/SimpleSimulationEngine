# https://stackoverflow.com/questions/18662261/fastest-implementation-of-sine-cosine-and-square-root-in-c-doesnt-need-to-b
# https://en.wikipedia.org/wiki/Bhaskara_I's_sine_approximation_formula

# https://docs.scipy.org/doc/numpy-1.13.0/reference/generated/numpy.polynomial.chebyshev.Chebyshev.fit.html

import numpy as np

'''

       p(x)     a0 + a1*x + a2*x^2 + a3*x^3 ...
f(x) = ---- = -----------------------------------
       q(x)     b0 + b1*x + b2*x^2 + b3*x^3 ...

df(x)/dai    =   (1/q(x)) * x^i
df(x)/dai    = -f(x)/q(x) * x^i

            p'(x)q(x) - p(x)q'(x)
df(x)/dx = -----------------------
            q(x)^2

Err  =  (Sum_i { f(xi) - yi })^2
dErr/dai =  2 * (Sum_i { f(xi) - yi }) df(x)/dai
dErr/dbi =  2 * (Sum_i { f(xi) - yi }) df(x)/dbi

'''

#def evalPoly

def evalRationalFunction( x, ps, qs, pis, qis ):
    p  = np.zeros( x.shape )
    q  = np.zeros( x.shape )
    xn = np.ones( x.shape )
    order = max( max(pis), max(qis) )+1
    ip=0; iq=0
    #print "order ", order
    for k in xrange(order):
        #print "k ", k
        if k in pis:
            #print "p[%i]" %k, ps[ip] * xn
            p += ps[ip] * xn 
            ip+=1
        if k in qis:
            #print "q[%i]" %k, qs[iq] * xn
            q += qs[iq] * xn 
            iq+=1
        xn*=xs
    return p,q

def rationalVarDeriv( x, dy, y, inv_q, pis, qis, fps, fqs ):
    #f    = p/q
    dEdp = dy*inv_q
    dEdq = dy*(-y*inv_q)
    xn = np.ones( x.shape )
    order = max( max(pis), max(qis) )+1
    ip=0; iq=0
    for k in xrange(order):
        if k in pis:
            fps[ip] = ( dEdp * xn ).sum()
            #print "fp[%i]" %k, fps[ip], dEdp * xn
            ip+=1
        if k in qis:
            fqs[iq] = ( dEdq * xn ).sum()
            #print "fp[%i]" %k, fqs[iq], dEdq * xn
            iq+=1
        xn*=xs
    return fps, fqs

def move( ps, fs, vs, dt=0.1, damp=1.0, invMass=1.0 ):
    dtv     = dt*invMass
    damping = min( 0.99, max( 0.5, 1-dtv*damp ) )
    #vs[:]  = vs*(1-damp) + fs*dtv
    vs[:]  = vs*damping + fs*dtv
    ps    += vs*dt
    #print "ps: ", ps
    #print "fs: ", fs
    #ps += fs * dt

def fitRational( x, y_ref, ps, qs, pis, qis,  nMaxIter=100, Fconv=1e-16, dt=0.1, damp=0.1, invMass=10.0, ws=None, Fconvs=None, Econvs=None, Fconvs_p=None, Fconvs_q=None ):
    F2conv = Fconv**2
    ps=np.array(ps)
    qs=np.array(qs)
    fps = np.zeros(len(ps))
    fqs = np.zeros(len(qs))
    vps = np.zeros(len(ps))
    vqs = np.zeros(len(qs))
    nx = len(x)
    if ws is None:
        ws = np.ones(nx) * (1./nx)
    for itr in xrange(nMaxIter):
        #print "=========== ITER ", itr
        p,q = evalRationalFunction( x, ps, qs, pis, qis )
        y  = p/q
        dy = y_ref - y
        E = (dy**2).sum()
        rationalVarDeriv( x, dy*ws, y, 1./q, pis, qis, fps,fqs )
        #fps*=-1; fqs*=-1
        F2err_p = (fps**2).sum()
        F2err_q = (fqs**2).sum()
        F2err = F2err_p + F2err_q
        if Econvs is not None: Econvs.append(E)
        if Fconvs is not None: Fconvs.append(F2err)
        if Fconvs_p is not None: Fconvs_p.append(F2err_p)
        if Fconvs_q is not None: Fconvs_q.append(F2err_q)
        print "iter %i F2err %g RMSE %g " %( itr, np.sqrt(F2err), np.sqrt(E) )
        #print "iter %i F2err %g " %(itr,np.sqrt(F2err) ), ps, fps  #, qs
        move( ps, fps, vps, dt=dt, damp=damp )
        move( qs, fqs, vqs, dt=dt, damp=damp, invMass=invMass )
        if F2err<F2conv:
            break
    #print Fconvs
    #print Econvs
    return ps, qs

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    np.set_printoptions(precision=16, linewidth=200 )

    # Tan = 1/cotg !!!! => enough to have -pi/4,pi/4


    # Tan approx
    #    -1/(x-1.5) -1/(x+1.5) 
    #   --( 2*x )/( x*x-(9./4) )
    #ps  = [ 0., -2., 0., 0.    ];        pis=[0,1,2,3]
    #qs  = [ -(np.pi/2)**2, 0., 1., 0  ]; qis=[0,1,2,3]

    #ps  = [ -2., 0.          ];  pis=[1,3];
    #qs  = [ -(np.pi/2)**2, 1. ]; qis=[0,2];

    #ps = [-2.3395812063314856,  0.1733543863012315, 0.0];   pis=[1,3,5];
    #qs = [-2.337150576473808,   0.9459337690826924];        qis=[0,2];

    #ps = [-2.3372519146471693,  0.1677758691235614,  0.0045067070833865, 0., 0. ];   pis=[1,3,5,7,9];
    #qs = [-2.3373184658660437,  0.9472075025917591, 0.0 ];  qis=[0,2,4];

    ps = [-2.3372901747581194e+00,  1.6791936299541979e-01,  4.4586164347535087e-03, -8.8472791163457613e-05,  8.1593603822823545e-05, 0.0 ];   pis=[1,3,5,7,9, 11 ]
    qs = [-2.3373135470026667e+00,  9.4718375532064869e-01,  5.0300549988835022e-05, 0.0 ];    qis=[0,2,4, 6]

    #ps  = [ 0., 0.5 ];   pis={0,1}
    #qs  = [ 1.  ];       qis={0}

    eps  = 0.0
    xmax = (np.pi/4)-eps

    eps  = 0.2
    xmax = (np.pi/2)-eps
    xs   = np.linspace( -xmax, xmax, 100 )
    #xs    = np.linspace( -np.pi/2+eps, np.pi/2-eps, 30 )
    #xs    = np.linspace( -np.pi/2+eps, np.pi/2-eps, 5 )
    #y_ref = 0.25 + 0.9*xs
    y_ref = np.tan(xs)


    p,q   = evalRationalFunction( xs, ps, qs, pis, qis )
    y_pq  = p/q
   
    y_pq_bak = y_pq.copy()

    #Fconvs = []
    Fconvs_p = []
    Fconvs_q = []
    Econvs = []

    #ps,qs = fitRational( xs, y_ref, ps, qs, pis, qis, nMaxIter=2000, Fconv=1e-16, dt=0.5, damp=2.0, invMass=0.2, Econvs=Econvs, Fconvs_p=Fconvs_p, Fconvs_q=Fconvs_q )
    ps,qs = fitRational( xs, y_ref, ps, qs, pis, qis, nMaxIter=2000, Fconv=1e-16, dt=0.1, damp=0.01, invMass=0.3, Econvs=Econvs, Fconvs_p=Fconvs_p, Fconvs_q=Fconvs_q )
    p,q   = evalRationalFunction( xs, ps, qs, pis, qis )
    y_pq  = p/q

    #print Fconvs
    #print Econvs

    plt.figure()
    plt.plot( np.sqrt(np.array(Fconvs_p)), "-", label='|Fp|' )
    plt.plot( np.sqrt(np.array(Fconvs_q)), "-", label='|Fq|' )
    plt.plot( np.sqrt(np.array(Econvs))  , "-", label='|E|'  )
    plt.yscale('log')
    plt.legend()
    plt.grid()

    print " Ps ", ps
    print " Qs ", qs

    plt.figure()
    plt.plot(xs, y_ref, "-", label='y_ref' )
    plt.plot(xs, y_pq,  ":", label='y_p/q' )
    plt.plot(xs, y_pq_bak,  ":", label='y_p/q_bak' )
    plt.legend()
    plt.grid()

    y_err = y_pq - y_ref
    plt.figure()
    plt.plot(xs, abs(y_err),  ":", label="y_err" )
    plt.yscale('log')

    plt.legend()
    plt.grid()
    plt.show()