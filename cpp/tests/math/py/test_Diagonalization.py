#!/usr/bin/python

import numpy as np
import libTest_Lingebra as LA
import JacobiEigen as JE
import matplotlib.pyplot as plt

n = 6;
#np.random.seed(15464)
sqrtA    = np.random.rand( n,n ) - 0.5
A        = np.dot( sqrtA, np.transpose(sqrtA) )





tol = 1e-9



V     = np.zeros((n,n))
es    = np.zeros(n)
mjs   = np.zeros(n,dtype=np.int32) 
ijmax = np.zeros(2,dtype=np.int32) 
A_=A.copy(); V_=V.copy();
print A

esA, vsA = np.linalg.eig( A )
print "eigenvalues np.eig:     ",np.sort(esA) 
esA,vsA = JE.jacobi    ( A_, tol = tol ) 
print "eigenvalues py.jacobi:  ",np.sort(esA)  
esA,vsA = LA.eig_Jacobi( A,  tol=tol )
print "eigenvalues cpp.jacobi: ",np.sort(esA)


'''
vmax=LA.lib.eig_Jacobi_init(n,A,V,mjs,ijmax)
mjs_,mvs_ = JE.initMaxs(A_)
print "mjs  ", mjs , vmax
print "mjs_ ", mjs_, mvs_.max()
for itr in range(100):
    print "============ itr ", itr,   (ijmax[0],ijmax[1])
    #print "--- py ---"
    vmax_,imax_,jmax_ = JE.maxElem(A_)
    JE.rotate    (A_,V_,imax_,jmax_)
    JE.updateMaxs(A_,   imax_,jmax_,mjs_,mvs_)
    vmax_,imax_,jmax_ = JE.maxElem(A_)
    #print "imax,jmax",(imax_,jmax_)
    #print "--- cpp ---" 
    #ijmax[0]=imax_; ijmax[1]=jmax_   # to check just rotation - works fine
    vmax=LA.lib.eig_Jacobi_step(n,A,V,mjs,ijmax,vmax)
    #print "mjs  ", mjs
    #print "mjs_ ", mjs_
    print "cpp ", (ijmax[0],ijmax[1]), mjs , vmax
    print "py  ", (imax_,jmax_)      , mjs_, vmax_
    if not np.array_equal(mjs, mjs_):
        print "ERROR :  mjs not equal "
        break;
    if (vmax<tol):
        print "converged by %i Jacobi rotations" %itr 
        break;
    #print ijmax,vmax,vmax
'''


plt.figure(figsize=(15,5))    
plt.subplot(1,3,1); plt.imshow( np.log10(np.abs(A   )+tol), interpolation="nearest", cmap='gray');
plt.subplot(1,3,2); plt.imshow( np.log10(np.abs(A_  )+tol), interpolation="nearest", cmap='gray');
plt.subplot(1,3,3); plt.imshow( np.log10(np.abs(A-A_)+tol), interpolation="nearest", cmap='gray');
plt.show()

