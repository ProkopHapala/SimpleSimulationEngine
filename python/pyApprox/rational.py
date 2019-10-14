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

def rationalVarDeriv( x, dy, y, q, pis, qis, fps, fqs ):
    #f    = p/q
    dEdp = dy/q
    dEdq = dy*(-y/q)
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
            iq+=1
        xn*=xs
    return fps, fqs

def move( ps, fs, vs, dt=0.1, damp=0.9 ):
    #vs[:]  = vs*damp + fs*dt
    #ps    += vs*dt
    ps += fs * dt

def fitRational( x, y_ref, ps, qs, pis, qis,  nMaxIter=100, Fconv=1e-16, dt=0.1, damp=0.9, ws=None ):
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
        p,q = evalRationalFunction( x, ps, qs, pis, qis )
        y  = p/q
        dy = y_ref - y
        rmse = (dy**2).sum()
        rationalVarDeriv( x, dy*ws, y, q, pis, qis, fps,fqs )
        #fps*=-1; fqs*=-1
        F2err = (fps**2).sum() + (fqs**2).sum()
        print "iter %i F2err %g RMSE %g " %( itr, np.sqrt(F2err), rmse )
        #print "iter %i F2err %g " %(itr,np.sqrt(F2err) ), ps, fps  #, qs
        move( ps, fps, vps, dt=dt, damp=damp )
        #move( qs, fqs, vqs, dt=dt, damp=damp )
        if F2err<F2conv:
            break
    return ps, qs

if __name__ == "__main__":
    import matplotlib.pyplot as plt
    np.set_printoptions(precision=16, linewidth=200 )

    # Tan approx
    #    -1/(x-1.5) -1/(x+1.5) 
    #   --( 2*x )/( x*x-(9./4) )
    ps  = [ 0., -2., 0., 0.    ];        pis=[0,1,2,3]
    qs  = [ -(np.pi/2)**2, 0., 1., 0  ]; qis=[0,1,2,3]

    #ps  = [ 0., 0.5 ];   pis={0,1}
    #qs  = [ 1.  ];       qis={0}

    eps = 0.2
    #xs    = np.linspace( -np.pi/2+eps, np.pi/2-eps, 100 )
    xs    = np.linspace( -np.pi/2+eps, np.pi/2-eps, 30 )
    #xs    = np.linspace( -np.pi/2+eps, np.pi/2-eps, 10 )
    #y_ref = 0.25 + 0.9*xs
    y_ref = np.tan(xs)


    p,q   = evalRationalFunction( xs, ps, qs, pis, qis )
    y_pq  = p/q
   
    y_pq_bak = y_pq.copy()


    ps,qs = fitRational( xs, y_ref, ps, qs, pis, qis, nMaxIter=100, Fconv=1e-16, dt=0.2, damp=0.95 )
    p,q   = evalRationalFunction( xs, ps, qs, pis, qis )
    y_pq  = p/q

    print " Ps ", ps
    print " Qs ", qs

    plt.plot(xs, y_ref, "-", label='y_ref' )
    plt.plot(xs, y_pq,  ":", label='y_p/q' )
    plt.plot(xs, y_pq_bak,  ":", label='y_p/q_bak' )

    plt.legend()
    plt.grid()
    plt.show()