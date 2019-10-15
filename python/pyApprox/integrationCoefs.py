
import numpy as np

'''

integral is basically scalar product between weights and function values
For test function fi(x) with function values F_ik = fi(x_k) we can write integral

Si = Sum_k{ F_ik w_k }
where w_k are integration weights

assuming we know the true value of the integral   Ii we want to minimize error
Ei = (Si - Ii)  for all i using a single set of weights w_k

That is basically least squere fitting

(F^T,F) w_k = F I

w_k = inv(F^T,F) F I

Usually the integration grids posses considerable symmetry, therefore there is only few non-equvalent w_k.
This symmetry can be used to considerably reduce the size of the problem

We simply sum the function values at equvalent points

'''

def Gauss2D( xs, p ):
    r2 = (xs[0,:]-p[0])**2 + (xs[1,:]-p[1])**2
    return np.exp( r2 * -p[2] ) * p[3]

def Lorenz2D( xs, p ):
    #print p
    r2 = (xs[0,:]-p[0])**2 + (xs[1,:]-p[1])**2
    return 1/( 1 + r2*p[2] ) * p[3]

def testFunc( xs, params, funcComp=Lorenz2D ):
    #print "xs.sh", xs.shape[1:]
    ys = np.zeros( xs.shape[1:] )
    #for i,p in enumerate(params):
    for p in params:
        #print " param[%i] " %i, p
        ys += funcComp( xs, p )
    return ys

def sample( params, Xs, ys=None ):
    #Ys = []
    if ys is None:
        ys = np.empty( len(Xs) )
    for i,xs in enumerate(Xs):
        #ys[i] = func( xs ).sum()
        ys[i] = testFunc( xs, params ).sum()
        #Ys.append( func(xs) )
    return ys

def findIntegrationWeights( funcs, Xs, xs_ref, dw_ref, xs_back=None, dw_back=None, errors=None, ws0=None ):
    nf = len(funcs)
    nw = len(Xs)
    I_ref = np.empty(nf)
    Ys    = np.empty( (nf,nw) )
    for i,func in enumerate(funcs):
        print "sample func[%i]" %i
        params = func
        #print params
        #I_ref[i]  = func( xs_ref  ).sum()      # reference integral using high-resolution grid
        I_ref[i]  = testFunc( xs_ref, params ).sum() * dw_ref      # reference integral using high-resolution grid
        if xs_back is not None:
            #Iref[i] -= func( xs_back ).sum()  # subtract integral over background points (surrounding, embeding)
            I_ref[i] -= testFunc( xs_back, params ).sum() * dw_back  # subtract integral over background points (surrounding, embeding)
        sample( params, Xs, Ys[i,:] )               # sample function at our new integration grid
    ws = np.linalg.lstsq( Ys, I_ref )[0]
    #ws0 = np.array([ dw_back, 0.0, 0.0, 0.0 ])
    if errors is not None:
        nerr = errors.shape[1]
        for i in xrange(nf):
            #Iref     = I_ref[i]
            #Iapprox  = np.dot( Ys[i,:], ws )
            #dI       = Iapprox-Iref
            errors[i,0] = I_ref[i]
            errors[i,1] = np.dot( Ys[i,:], ws )
            if nerr>2:
                errors[i,2] = np.dot( Ys[i,:], ws0 )
            #print "Func[%i] " %i, dI/Iref, Iref, Iapprox, " | ", (Iapprox_-Iref)/Iref
    return ws

def testIntegrationWeights( ws, funcs, Xs, xs_ref, dw_ref, xs_back=None, dw_back=None, ws0=None ):
    nf = len(funcs)
    nw = len(Xs)
    I_ref = np.empty(nf)
    Ys    = np.empty( (nf,nw) )
    for i,func in enumerate(funcs):
        print "test sample func[%i]" %i
        params = func
        I_ref[i]  = testFunc( xs_ref, params ).sum() * dw_ref      # reference integral using high-resolution grid
        if xs_back is not None:
            I_ref[i] -= testFunc( xs_back, params ).sum() * dw_back  # subtract integral over background points (surrounding, embeding)
        sample( params, Xs, Ys[i,:] )               # sample function at our new integration grid
    nerr=2
    if ws0 is not None:
        nerr+=1
    errors = np.empty( (nf, nerr) )
    nerr = errors.shape[1]
    for i in xrange(nf):
        errors[i,0] = I_ref[i]
        errors[i,1] = np.dot( Ys[i,:], ws )
        if nerr>2:
            errors[i,2] = np.dot( Ys[i,:], ws0 )
    return errors

def genSamplePointRect( n, xmin, xmax, offset=0.5 ):
    dx = (xmax-xmin)/n
    xs = np.arange( xmin, xmax, dx ) + dx*offset
    Xs,Ys = np.meshgrid(xs,xs)
    return np.array( [ Xs.flat.copy(), Ys.flat.copy() ] )

def index2mask( n, inds ):
    mask = np.zeros(n,dtype=np.bool)
    mask[inds] = True
    return mask

def splitByMask( xs, mask ):
    mask_  = np.logical_not( mask )
    return xs[:,mask], xs[:,mask_]

def splitByInds( xs, inds ):
    mask = index2mask( xs.shape[1], inds )
    return splitByMask( xs, mask )

def splitByRect( xs, rect ):
    mask  = np.logical_and( np.logical_and( xs[0,:]>rect[0], xs[0,:]<rect[1] ), np.logical_and( xs[1,:]>rect[2], xs[1,:]<rect[3] ) )
    return splitByMask( xs, mask )

def getParams( nfunc, nfcomp, xmax=0.5, sigmaMin=0.1, sigmaMax=0.5 ):
    params = np.random.rand(nfunc,nfcomp,4 )
    params[:,:,0] = params[:,:,0]*xmax*2 - xmax
    params[:,:,1] = params[:,:,1]*xmax*2 - xmax
    #params[:,:,2] = 1/( 0.004 + 0.125*params[:,:,2]**2 )
    params[:,:,2] = 1/( sigmaMin*sigmaMin + (sigmaMax*sigmaMax-sigmaMin*sigmaMin)*params[:,:,2]**2 )
    params[:,:,3] = params[:,:,3]*2 - 1.0
    return params


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    np.set_printoptions(precision=16, linewidth=200 )

    xmax   = 1.0
    xsub   = 0.5
    extent = ( -xmax,xmax,  -xmax,xmax )

    nref  = 128
    nback = 16
    dx_ref  = 2*xmax/nref
    dx_back = 2*xmax/nback
    dw_ref  = dx_ref **2
    dw_back = dx_back**2

    print dw_ref, dw_back

    xs_ref  = genSamplePointRect( nref,  -xmax, xmax )
    xs_back = genSamplePointRect( nback, -xmax, xmax )

    #mask1 = np.logical_and( np.logical_and( xs_back[0,:]>-xsub, xs_back[0,:]<xsub ), np.logical_and( xs_back[1,:]>-xsub, xs_back[1,:]<xsub ) )
    #mask2 = np.logical_and( np.logical_and( xs_back[0,:]>-xsub, xs_back[0,:]<0.6 ), np.logical_and( xs_back[1,:]>-xsub, xs_back[1,:]<.6 ) )
    #mask  = np.logical_not( mask1 )
    #xs_sub1 = xs_back[:,mask1 ]
    #xs_sub2 = xs_back[:,mask2 ]; xs_sub2[0,:]-=1./16; xs_sub2[1,:]-=1./16
    #xs_back = xs_back[:,mask  ]


    xs_sub2, _   = splitByRect( xs_back, (-xsub,xsub-dx_back,-xsub,xsub-dx_back) )
    xs2          = xsub - dx_back * 0.5
    xs_sub2[0,:]+=dx_back*0.5; xs_sub2[1,:]+=dx_back*0.5 
    #xs_sub2, xs_border  = splitByRect( xs_sub2, (-xs2,xs2,-xs2,xs2 ) )
    #print "xs_border.shape ", xs_border.shape
    #xs_corner,xs_edge = splitByInds( xs_border, [0,8,23,-1] )
    xs_sub, xs_back = splitByRect( xs_back, (-xsub,xsub,-xsub,xsub) )

    xs2          = xsub - dx_back * 0.5
    xs_sub, xs_border  = splitByRect( xs_sub, (-xs2,xs2,-xs2,xs2 ) )
    xs_corner,xs_edge = splitByInds( xs_border, [0,7,20,-1] )

    Xs = [ xs_sub2, xs_sub, xs_corner,xs_edge ]

    nfunc  = 1000
    nfcomp = 100

    params = getParams( nfunc, nfcomp, xmax=0.4, sigmaMin=0.06, sigmaMax=0.5 )

    errors = np.empty((nfunc,3))
    #ws0    =  np.array([ dw_back, 0.0, 0.0, 0.0 ])
    ws0    =  np.array([ 0, dw_back, dw_back, dw_back ])

    ws = findIntegrationWeights( params, Xs, xs_ref, dw_ref, xs_back=xs_back, dw_back=dw_back, errors=errors, ws0=ws0 )
    ws = np.array(ws)

    params_test = getParams( 100, nfcomp, xmax=0.4, sigmaMin=0.5, sigmaMax=1.0 )

    errors =  testIntegrationWeights( ws, params_test, Xs, xs_ref, dw_ref, xs_back=xs_back, dw_back=dw_back, ws0=ws0 )
    print "ws", ws
    print "ws/dw_back", ws/dw_back

    params = params_test

    Ys = testFunc( xs_ref, params[0] )

    plt.imshow( Ys.reshape(nref,nref), extent=extent, cmap='gray' )

    #plt.plot( xs_ref[0,:], xs_ref[1,:], '.', markersize=1.0 )
    plt.plot( xs_back[0,:], xs_back[1,:], '+' )
    #plt.plot( xs_sub1[0,:], xs_sub1[1,:], '+' )
    #plt.plot( xs_sub2[0,:], xs_sub2[1,:], 'x' )
    for xs in Xs:
        plt.plot( xs[0,:], xs[1,:], '.',  markersize=2.0 )
    plt.axis('equal')

    #plt.show(); exit()



    plt.figure()
    err  = (errors[:,1]-errors[:,0])
    err0 = (errors[:,2]-errors[:,0])

    rmse  = np.sqrt( (err**2).sum() )
    rmse0 = np.sqrt( (err0**2).sum() )
    print " rmse, rmse0 ", rmse, rmse0

    err_rel  = err/errors[:,0]
    err_rel0 = err0/errors[:,0]
    plt.plot( abs(err ), '.', label='opt(abs)' , markersize=2.0 )
    plt.plot( abs(err0), '.', label='prev(abs)', markersize=2.0 )
    #plt.plot( abs(err_rel ), '.-', label='opt(rel)'  )
    #plt.plot( abs(err_rel0), '.-', label='prev(rel)' )
    plt.yscale( 'log' )
    plt.grid()
    plt.legend()
    plt.show()
