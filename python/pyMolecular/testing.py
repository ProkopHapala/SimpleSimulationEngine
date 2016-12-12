#!/usr/bin/python
import numpy as np

def saw_sine( xs ):
    ixs  = xs.astype(int)
    dxs  = xs - ixs
    mask = ( np.bitwise_and( ixs , 1 ) == 1 ) 
    dxs[mask] = 1-dxs[mask]
    return dxs

def hash_saw( xs, k ):
    '''
       expected properties of this hash function are:
       for 'N' atoms displaced by 'dx'  the mean change of hash <dh> ~ dx * sqrt(N)
    '''
    return np.sum( saw_sine((xs+100.0)*k)-0.5 )
    #return np.sum( np.sin( 2*np.pi*xs ) )         # this should be almost the same, but the above is a bit faster to compute in C++
    
def mutateN( xs, N, dx ):
    inds = (np.random.rand( N ) * len(xs)).astype(int)
    xs[inds] += (2*dx) * ( np.random.rand( N )-0.5 )
    return xs
