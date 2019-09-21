
import numpy as np

np.set_printoptions(linewidth=200)

# 2D integrals with simpson rule
#  http://mathfaculty.fullerton.edu/mathews/n2003/SimpsonsRule2DMod.html

import quadratureCoefs as qc

def gauss( r, b=1.0):
    return np.exp( -b*r*r )

def slater( r, b=1.0):
    return np.exp( -b*r )

def GR2( r, b=1.0):
    return (1-r*r)**2

rfunc_default = slater

def getFuncs( xs, ys, f1=rfunc_default, f2=rfunc_default ):
    Xs,Ys = np.meshgrid(xs,ys)
    R2s   = Xs**2 + Ys**2
    Rs    = np.sqrt(R2s)
    Ws    = Ys * (2*np.pi)
    f1s   = f1(Rs)
    f2s   = f2(Rs)
    return f1s, f2s, Ws

def intRfNumpy_brute( f1, f2, Ws, out=None ):
    if out is None: out=np.empty( Ws.shape[1] )
    for i in xrange(len(out)):
        f2_    = np.roll(f2,i,axis=1)
        f1f2   = f1 * f2_ * Ws
        out[i] = f1f2.sum()
    return out

def intRfNumpy_fft( f1, f2, Ws, out=None ):
    if out is None: out=np.empty( Ws.shape[1] )
    f1w    = np.fft.fftn(f1,axes=(1,))
    f2w    = np.fft.fftn(f2,axes=(1,))
    f1f2w  = f1w*f2w
    f1f2   = np.fft.ifftn(f1f2w,axes=(1,))
    out[:] = (f1f2*Ws).sum(axis=0)
    return out

if __name__ == "__main__":
    import matplotlib.pyplot as plt

    ymax = 5.0

    xs    = np.arange(-10.,10.,0.1)
    ys    = np.arange(  0., ymax, 0.1 )
    #ys    = np.linspace(  0., ymax, 50 )
    
    func1 = slater
    func2 = slater
    
    dx    = xs[1]-xs[0]
    dy    = ys[1]-ys[0]

    f1,f2,Ws = getFuncs( xs, ys, f1=func1, f2=func2 )

    I_brute = intRfNumpy_brute( f1, f2, Ws*(dx*dy) ); plt.plot(xs,I_brute,label='brute')
    I_fft   = intRfNumpy_fft  ( f1, f2, Ws*(dx*dy) ); plt.plot(xs,I_fft  ,label='fft')

    xs_    = np.arange(-10.,10.,0.5)
    dx_    = xs_[1]-xs_[0]
    f1,f2,Ws = getFuncs( xs_, ys, f1=func1, f2=func2 )
    I_low = intRfNumpy_fft( f1, f2, Ws*(dx_*dy) );   plt.plot(xs_,I_low,':',label='fft_low')

    order = 6
    ys = np.array(qc.GaussLegendreNodes  [order])*ymax
    ws = np.array(qc.GaussLegendreWeights[order])

    f1,f2,Ws = getFuncs( xs, ys )
    Ws*=ws[:,None]
    I_CbG = intRfNumpy_fft( f1, f2, Ws*dx*ymax );   plt.plot(xs,I_CbG,':',label='fftCheby')
    
    
    f1,f2,Ws = getFuncs( xs_, ys, f1=func1, f2=func2 )
    I_CbG_low = intRfNumpy_fft( f1, f2, Ws*dx*ymax );   plt.plot(xs,I_CbG,':',label='fftCheby_low')

    ratio = I_CbG/I_brute

    #plt.plot(xs,(ratio-1.0)*100.0,label='error_ratio')
    #print ratio

    plt.legend()
    plt.grid()
    plt.show()

