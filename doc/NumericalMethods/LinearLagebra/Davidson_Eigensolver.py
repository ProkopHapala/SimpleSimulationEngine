#!/bin/python

# from Here:  https://joshuagoings.com/2013/08/23/davidsons-method/
# Joshua Goings 

'''
### Davidson's Method

Many times in quantum chemistry you want the lowest few eigenvalues of very large matrices. 
For example, in my work with quantum dots, we want the first few excitation energies of our dots from a TDDFT calculation. 
Our smallest systems have 1000 basis functions with roughly 300 electrons, which requires creating and diagonalizing a 210,000 by 210,000 matrix (the dimension is number occupied orbitals times the number unoccupied orbitals). 
This is way too huge! For perspective, just storing this matrix on disk requires about 300 GB space.

Ernest Davidson ran into this problem decades ago working on diagonalizing the Configuration Interaction Hamiltonians, which are generally enormous. 
He came up with a method — the so-called Davidson method – which iteratively diagonalizes a subspace of the matrix instead of the whole thing, and gives you the first few lowest (or highest!) eigenvalues. 
It is much more cost-efficient, and actually doesn’t require you to create the whole matrix in the first place (it projects the matrix onto an appropriate subspace instead). 
The method turned out to work so well that quantum chemists adopted it and have used it ever since.

I wanted to try it out, so I implemented Davidson’s method on a Hermitian matrix (in Python, of course :)).

Let’s step through what I’ve done and see if it makes this method any clearer. The first bit simply creates our fake Hamiltonian.

sparsity = 0.0001
A = np.zeros((n,n))
for i in range(0,n):
    A[i,i] = i + 1 
A = A + sparsity*np.random.randn(n,n) 
A = (A.T + A)/2   

While it may look arbitrary, take a closer look at the structure. First, it is diagonally dominant. 
The diagonal is filled with increasing integers, while the off-diagonals are random numbers multiplied by a scaling factor to “mute” them somewhat. 
This adds sparsity. Davidson’s method really excels with sparse, diagonally dominant matrices. 
This is actually very similar to the Hamiltonians we encounter as quantum chemists. 
If you scale the sparsity down and approach zero, the Davidson method speeds up fast. 
(But don’t set it to exactly zero or it will crash — in this case the matrix is already diagonalized!)

Next we set up our subspace “trial vectors”:
  
k = 8 # number of initial guess vectors  
eig = 4 # number of eignvalues to solve  
t = np.eye(n,k) # set of k unit vectors as guess  
V = np.zeros((n,n)) # array of zeros to hold guess vec  
I = np.eye(n) # identity matrix same dimen as A  
Since we are choosing to find the first four eigenvalues, we need at least four guess vectors k. In practice, we choose maybe twice to three times that, because we want to increase the span of our guess space. In other words, it helps us hone in on the appropriate eigenvectors faster. But don’t make the guess too big! If it gets too large, we basically end up diagonalizing the whole matrix — which we don’t want to do, since that is the whole point of Davidson’s method. Speaking of the span of the guess vectors, it is important to make a good initial guess. Because the matrix is diagonally dominant, I chose a set of unit vectors as my guess, which is a good since our matrix is so close to being scalar multiples of the identity matrix.

Finally we get to the meat of the main routine:

...

1. We first check to see if this is our first iteration (e.g. m < k). 
   - If it is, we add our guess vectors to our set V, and set theta_old equal to one. 
      - theta_one equal to one is arbitrary, and is used to ensure we don’t “converge” on our first try (we can’t, since we have to compare at least two iterations to determine convergence). 
   - If this isn’t our first iteration, we set theta_old to our last iteration’s eigenvalues.
2. Next we ensure our vectors are orthonormal among each other with numpy’s QR decomposition routine. 
   - Vectors must be orthonormal or the routine will fail. Then we project our matrix A onto the subspace defined by our guess vectors, V^T A V
3. Then we diagonalize this matrix projection subspace, sort the eigenvalues to find the lowest few desired, and compute the residual r. 
   - The residual will be zero if the guess eigenvectors are the real eigenvectors. 
     - => If the residual is lower than some criterion, then we quit. 
     -Else, we orthonormalize the eigenvectors from that iteration and add them to our guess vectors! 
     - In this way, our subspace grows each time we iterate…the hope is that each added vector will be in the span of the real eigenvectors. 
     - Once that happens, we will have solved the problem of getting the lowest few eigenvalues.

Here is some sample output:

davidson = [0.99999921 2.00000133 3.00000042 3.99999768] ; 0.6596 seconds  
numpy    = [0.99999921 2.00000133 3.00000042 3.99999768] ; 1.7068 seconds
You can see it works! (And the eigenvalues are really similar to the integer values we put along the diagonal).

I’ve attached the full routine at the end. With a sparse enough matrix, I can beat numpy by about a second. Of course, comparisons aren’t really fair, since I make numpy compute all the eigenvalues. 
From what I know, this method has a lot of intricacies, and I am still learning many of them. If you know of any ways I can improve this code, let me know! Comments and messages are always appreciated.

'''

from __future__ import division
from __future__ import print_function
import math
import numpy as np
import time

''' Block Davidson, Joshua Goings (2013)

    Block Davidson method for finding the first few
	lowest eigenvalues of a large, diagonally dominant,
    sparse Hermitian matrix (e.g. Hamiltonian)
'''

n = 1200					# Dimension of matrix
tol = 1e-8				# Convergence tolerance
mmax = n//2				# Maximum number of iterations	

''' Create sparse, diagonally dominant matrix A with 
	diagonal containing 1,2,3,...n. The eigenvalues
    should be very close to these values. You can 
    change the sparsity. A smaller number for sparsity
    increases the diagonal dominance. Larger values
    (e.g. sparsity = 1) create a dense matrix
'''

sparsity = 0.0001
A = np.zeros((n,n))
for i in range(0,n):
    A[i,i] = i + 1 
A = A + sparsity*np.random.randn(n,n) 
A = (A.T + A)/2 


k = 8					# number of initial guess vectors 
eig = 4					# number of eignvalues to solve 
t = np.eye(n,k)			# set of k unit vectors as guess
V = np.zeros((n,n))		# array of zeros to hold guess vec
I = np.eye(n)			# identity matrix same dimen as A

# Begin block Davidson routine

start_davidson = time.time()

for m in range(k,mmax,k):
    if m <= k:
        for j in range(0,k):
            V[:,j] = t[:,j]/np.linalg.norm(t[:,j])
        theta_old = 1 
    elif m > k:
        theta_old = theta[:eig]
    V[:,:m],R = np.linalg.qr(V[:,:m])
    T = np.dot(V[:,:m].T,np.dot(A,V[:,:m]))
    THETA,S = np.linalg.eig(T)
    idx = THETA.argsort()
    theta = THETA[idx]
    s = S[:,idx]
    for j in range(0,k):
        w = np.dot((A - theta[j]*I),np.dot(V[:,:m],s[:,j])) 
        q = w/(theta[j]-A[j,j])
        V[:,(m+j)] = q
    norm = np.linalg.norm(theta[:eig] - theta_old)
    if norm < tol:
        break

end_davidson = time.time()

# End of block Davidson. Print results.

print("davidson = ", theta[:eig],";",
    end_davidson - start_davidson, "seconds")

# Begin Numpy diagonalization of A

start_numpy = time.time()

E,Vec = np.linalg.eig(A)
E = np.sort(E)

end_numpy = time.time()

# End of Numpy diagonalization. Print results.

print("numpy = ", E[:eig],";",
     end_numpy - start_numpy, "seconds") 