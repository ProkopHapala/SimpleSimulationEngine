
import numpy as np

def getSamples( f, n, xmin=0, dx=1.0, fddx=0.1 ):
    ddx=dx*fddx
    xs = np.arange( xmin-dx, xmin+dx*(n+2)+1e-8, dx )
    xs[0 ]=xs[ 1]+ddx
    xs[-1]=xs[-2]-ddx
    #print xs
    ys  = f(xs)
    dy0 = (ys[0 ] - ys[1 ])/ddx
    dy1 = (ys[-2] - ys[-1])/ddx
    ys[ 0] = ys[2 ] - dy0*2*dx
    ys[-1] = ys[-3] + dy1*2*dx
    xs[0 ]=xs[ 1]-dx
    xs[-1]=xs[-2]+dx
    #print "xs ", xs
    #print "ys ", ys
    return ys, xs

def getSamplesDeriv( f, n, xmin=0, dx=1.0, fddx=0.1 ):
    ddx=dx*fddx
    xs  = np.arange( xmin, xmin+dx*(n+1)+1e-8, dx )
    ys  = f(xs-ddx)
    ys_ = f(xs+ddx)
    dys = (ys_-ys)/(2*ddx)
    #ys  = 0.5*(ys+ys_)
    ys  = f(xs)
    return ys, dys, xs

def hermite(x, y0,y1, dy0,dy1 ):
    y01 = y0-y1;
    return (     y0
        +x*(           dy0
        +x*( -3*y01 -2*dy0 - dy1
        +x*(  2*y01 +  dy0 + dy1 ))));

def evalValTable( x, ys, xstep ):
    #x+=xstep*2
    invstep = 1/xstep
    s = x*invstep + 1
    i = s.astype(dtype=np.int)
    #print "i ", i
    y0  = ys[i  ]
    y1  = ys[i+1]
    dy0 = 0.5 * ( ys[i+1] - ys[i-1] )
    dy1 = 0.5 * ( ys[i+2] - ys[i  ] )
    return hermite( s-i, y0,y1, dy0,dy1 )

def evalValTableDeriv( x, ys, dys, xstep ):
    invstep = 1/xstep
    s = x*invstep
    i = s.astype(dtype=np.int)
    return hermite( s-i, ys[i],ys[i+1], dys[i]*xstep, dys[i+1]*xstep )


