#!/usr/bin/python

import numpy as np
import libTest_Lingebra as LA
import JacobiEigen as JE
import matplotlib.pyplot as plt
import time

'''
test tol=1e-9
Lenovo-ideapad-Y700-15ISK    i7-6700HQ CPU @ 2.60GHz
N       time np.eig(A) [s]     time LA.eig_Jacobi(A,tol=1e-9)   Nrot
200         0.259829                0.902322
200         0.256582                0.900702
200         0.260135                0.875203        78648  => 11128 ns/rotation => 55ns/index
100         0.035707                0.1576
100         0.05374                 0.252058
100         0.061157                0.180691        18983  => 9518 ns/rotation =>  95ns/index
'''

n = 200;
#np.random.seed(15464)
sqrtA    = np.random.rand( n,n ) - 0.5
A        = np.dot( sqrtA, np.transpose(sqrtA) )

tol = 1e-9

V     = np.zeros((n,n))
es    = np.zeros(n)
mjs   = np.zeros(n,dtype=np.int32) 
ijmax = np.zeros(2,dtype=np.int32) 
A_=A.copy(); V_=V.copy();
#print A

t1=time.clock()
esA, vsA = np.linalg.eig( A )
t2=time.clock()
print "eigenvalues np.eig:     ",np.sort(esA) 
print "time:     ",t2-t1

t1=time.clock()
esA,vsA = LA.eig_Jacobi( A,  tol=tol, nMaxIter=100000 )
t2=time.clock()
print "eigenvalues cpp.jacobi: ",np.sort(esA)
print "time:     ",t2-t1

'''
esA,vsA = JE.jacobi    ( A_, tol = tol ) 
print "eigenvalues py.jacobi:  ",np.sort(esA)  
'''

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

