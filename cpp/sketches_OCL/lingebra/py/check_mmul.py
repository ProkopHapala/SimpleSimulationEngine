#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

A   =np.genfromtxt("A.dat")
B   =np.genfromtxt("B.dat")
C   =np.genfromtxt("C.dat")
Cref=np.genfromtxt("Cref.dat")

Cnp=np.dot(A,B)

print "C    vs Cnp error^2 : ", np.sum( (C    - Cnp)**2 )
print "Cref vs Cnp error^2 : ", np.sum( (Cref - Cnp)**2 )

plt.figure(figsize=(8,12))

plt.subplot(3,2,1); plt.imshow(A,    interpolation='nearest'); plt.colorbar(); plt.title("A")
plt.subplot(3,2,2); plt.imshow(B,    interpolation='nearest'); plt.colorbar(); plt.title("B")

plt.subplot(3,2,3); plt.imshow(C,     interpolation='nearest'); plt.colorbar(); plt.title("C")
plt.subplot(3,2,4); plt.imshow(Cnp,   interpolation='nearest'); plt.colorbar(); plt.title("Cnp")
#plt.subplot(3,2,4); plt.imshow(Cref,  interpolation='nearest'); plt.colorbar(); plt.title("Cref")

#plt.subplot(3,2,3); plt.imshow(Cref, interpolation='nearest'); plt.colorbar(); plt.title("Cref")
#plt.subplot(3,2,4); plt.imshow(Cnp,  interpolation='nearest'); plt.colorbar(); plt.title("Cnp")

plt.subplot(3,2,5); plt.imshow(Cref - Cnp,  interpolation='nearest'); plt.colorbar(); plt.title("err Cref")
plt.subplot(3,2,6); plt.imshow(C    - Cnp,  interpolation='nearest'); plt.colorbar(); plt.title("err C")

plt.show()
