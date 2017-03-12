#!/usr/bin/python

''' lam,x = jacobi(a,tol = 1.0e-9).
    Solution of std. eigenvalue problem [a]{x} = lam{x}
    by Jacobi's method. Returns eigenvalues in vector {lam}
    and the eigenvectors as columns of matrix [x].
'''

import numpy as np
from numpy import array,identity,diagonal
from math import sqrt

import matplotlib.pyplot as plt

# ================ functions

def maxElem(a):
    n = len(a)
    aMax = 0.0
    for i in range(n-1):
        for j in range(i+1,n):
            abs_aij = abs(a[i,j])
            if  abs_aij >= aMax:
                aMax = abs_aij
                imax = i; jmax = j
    return aMax,imax,jmax

def maxElem_On(a,mjs,mvs):
    n = len(a)
    imax=0
    vmax=mvs[imax]
    for i in range(1,n):
        vi=mvs[i]
        if vi>vmax:
            vmax=vi
            imax=i
    return vmax,imax,mjs[imax]

def maxInRow(a,i):
    n = len(a)
    jmax = i+1
    vmax = abs(a[i,jmax])
    for j in range(i+2,n):
        vij = abs(a[i,j])
        if vij >= vmax:
            vmax=vij
            jmax=j
    return jmax,vmax
    
def updateMaxs(a,imax,jmax,mjs,mvs):
    n = len(a)
    #print mjs,(imax,jmax)
    #mjs[imax],mvs[imax]=maxInRow(a,imax) # this row-wise step should not be required since column-wise procedure will take care of [imax] 
    for i in range(n-1):
        #if( (mj == jmax) or (mj==imax) ):
        mj   = mjs[i]
        #print  "%i,%i  %i,%i " %( i,mj, imax,jmax )
        if (i==imax)or(i==jmax)or(mj==imax)or(mj==jmax):
            #print i,"case 1"
            mjs[i],mvs[i]=maxInRow(a,i)
            continue
        vmj  = mvs[i]
        if (imax>i):
            v=abs(a[i,imax])
            if(v>vmj):
                vmj=v
                mj =imax
            #print i,"case imax",mj,vmj
        if (jmax>i):
            v=abs(a[i,jmax])
            if(v>vmj):
                vmj=v
                mj =jmax   
            #print i,"case jmax",mj,vmj    
        #mj,vmj=maxInRow(a,i)
        mjs[i]=mj;mvs[i]=vmj
        
        # ---- check update correctness ( On^2 )
        j,v=maxInRow(a,i)
        if (mjs[i]!=j) or (mvs[i]!=v):
            print i,(imax,jmax)," wrong: ",(mjs[i],mvs[i])," correct: ", (j,v) 
            exit(0)

def initMaxs(a):
    n = len(a)
    mvs=np.zeros(n)
    mjs=np.zeros(n,dtype=int)
    for i in range(n-1):
        mjs[i],mvs[i]=maxInRow(a,i)
    return mjs,mvs
    
def rotate(a,p,k,l): # Rotate to make a[k,l] = 0
    #  https://en.wikipedia.org/wiki/Jacobi_rotation
    n = len(a)
    aDiff = a[l,l] - a[k,k]
    if abs(a[k,l]) < abs(aDiff)*1.0e-36: t = a[k,l]/aDiff
    else:
        phi = aDiff/(2.0*a[k,l])
        t = 1.0/(abs(phi) + sqrt(phi**2 + 1.0))
        if phi < 0.0: t = -t
    c = 1.0/sqrt(t**2 + 1.0); s = t*c
    tau = s/(1.0 + c)
    temp = a[k,l]
    a[k,l] = 0.0
    a[k,k] = a[k,k] - t*temp
    a[l,l] = a[l,l] + t*temp
    for i in range(k):      # Case of i < k
        temp = a[i,k]
        a[i,k] = temp   - s*(a[i,l] + tau*temp)
        a[i,l] = a[i,l] + s*(temp   - tau*a[i,l])
    for i in range(k+1,l):  # Case of k < i < l
        temp = a[k,i]
        a[k,i] = temp   - s*(a[i,l] + tau*a[k,i])
        a[i,l] = a[i,l] + s*(temp   - tau*a[i,l])
    for i in range(l+1,n):  # Case of i > l
        temp = a[k,i]
        a[k,i] = temp   - s*(a[l,i] + tau*temp)
        a[l,i] = a[l,i] + s*(temp   - tau*a[l,i])
    for i in range(n):      # Update transformation matrix
        temp = p[i,k]
        p[i,k] = temp   - s*(p[i,l] + tau*p[i,k])
        p[i,l] = p[i,l] + s*(temp   - tau*p[i,l])    

def rotate_affected(a,k,l, val = 1.0 ): 
    a[k,l] = val
    a[k,k] = val
    a[l,l] = val
    for i in range(k):      # Case of i < k
        a[i,k] = val
        a[i,l] = val
    for i in range(k+1,l):  # Case of k < i < l
        a[k,i] = val
        a[i,l] = val
    for i in range(l+1,n):  # Case of i > l
        a[k,i] = val
        a[l,i] = val
 
def jacobi(a,tol = 1.0e-9, debug_func=None): # Jacobi method
    n = len(a)
    maxRot = 5*(n**2)       # Set limit on number of rotations
    p = identity(n)*1.0     # Initialize transformation matrix
    mjs,mvs = initMaxs(a)
    for i in range(maxRot): # Jacobi rotation loop 
        aMax,imax,jmax = maxElem_On(a,mjs,mvs)
        if(True):
            aMax_,imax_,jmax_ = maxElem(a)
            #print " On",(aMax,imax,jmax)," On2", (aMax_,imax_,jmax_)
            if (imax!=imax_) or (jmax!=jmax_):
                print "ERROR found wrong pivot "
                print " On",(aMax,imax,jmax)," On2", (aMax_,imax_,jmax_)
                exit(0)
        if aMax < tol:
            print "converged by %i rotations; error < %e " %(i,  tol)
            return diagonal(a),p
        rotate(a,p,imax,jmax)
        updateMaxs(a,imax,jmax,mjs,mvs)
        if debug_func is not None:
            #debug_func(a,p,imax,jmax)
            aff=np.zeros(a.shape)
            rotate_affected(aff,imax,jmax)
            debug_func(a,aff,imax,jmax)
    print 'Jacobi method did not converge'

def debug_plot(a,p,k,l):
    #plt.subplot(1,2,1); plt.imshow(a, interpolation="nearest")
    plt.subplot(1,2,1); plt.imshow( np.log10(np.abs(a)+1e-9), interpolation="nearest", cmap='gray')
    plt.subplot(1,2,2); plt.imshow( p, interpolation="nearest")
    #plt.subplot(1,2,2); plt.imshow( np.log(np.abs(p)), interpolation="nearest")
    #plt.subplot(1,2,2); plt.imshow( aff, interpolation="nearest")
    plt.show()
    plt.close()

# ================ Body
    
if __name__ == "__main__":
    n = 6;
    sqrtA    = np.random.rand( n,n ) - 0.5
    A        = np.dot( sqrtA, np.transpose(sqrtA) )
    esA, vsA = np.linalg.eig( A )
    print "eigenvalues np.eig: ",np.sort(esA) 
    #esA_, vsA_ = jacobi(A, tol = 1.0e-9, debug_func=debug_plot)    
    esA_, vsA_ = jacobi(A, tol = 1.0e-9, debug_func=None)   
    print "eigenvalues jacobi: ",np.sort(esA_)
    
    
    
