#!/usr/bin/python

'''

We want to find coefficient of expantion of function which is linear combination of powers of R^2 since R^2 = x^2 + y^2 + z^2 (+w^2 ... ) is easy to compute

f(x) = Sum_i  ai * x^(2*i)

We have 4 boundary conditions
f(x1)=y1,  f'(x1)=dy1
f(x2)=y1,  f'(x2)=dy2

for f(x) = a1 + a2*x**2 + a3*x**4 + a4*x**6 simpy gives following : 

coef 0:  -(  dy0*x0**3*x1**4 - dy0*x0*x1**6    + x0**6*(dy1*x1 - 2*y1)  +       x0**4*(-dy1*x1**3 + 6*x1**2*y1) - 6*x0**2*x1**4*y0 + 2*x1**6*y0)                                  /(2*(x0 - x1)**3*(x0 + x1)**3)
coef 1:   (2*dy0*x0**4*x1**3 - dy0*x0**2*x1**5 - dy0*x1**7              +   dy1*x0**7  + dy1*x0**5*x1**2 + x0**3*(-2*dy1*x1**4 - 12*x1**3*y0 + 12*x1**3*y1))                      /(2*x0*x1*(x0 - x1)**3*(x0 + x1)**3)
coef 2:  -(  dy0*x0**4*x1    + dy0*x0**2*x1**3 - 2*dy0*x1**5            + 2*dy1*x0**5 + x0**3*(-dy1*x1**2 - 6*x1*y0 + 6*x1*y1) + x0*(-dy1*x1**4 - 6*x1**3*y0 + 6*x1**3*y1))       /(2*x0*x1*(x0 - x1)**3*(x0 + x1)**3)
coef 3:   (  dy0*x0**2*x1    - dy0*x1**3       + dy1*x0**3              +       x0*(-dy1*x1**2 - 4*x1*y0 + 4*x1*y1))                                                              /(2*x0*x1*(x0 - x1)**3*(x0 + x1)**3)

'''

def analytical_coefs_sympy():
	import sympy
	from sympy.matrices import Matrix
	x0, x1, y0, y1, dy0, dy1, x = sympy.symbols('x0 x1 y0 y1 dy0 dy1 x')
	#a1,a2,a3,a4                = sympy.symbols('a0 a1 a2 a3')
	# ===== setup 
	bas1=1+0*x
	bas2=x**2
	bas3=x**4
	bas4=x**6
	# ===== Main
	dbas1=sympy.diff(bas1, x); dbas2=sympy.diff(bas2, x); dbas3=sympy.diff(bas3, x); dbas4=sympy.diff(bas4, x)
	#print ( " dbas1: ", dbas1); print ( " dbas2: ", dbas2); print ( " dbas3: ", dbas3); print ( " dbas4: ", dbas4)
	A = Matrix([
		[  bas1.subs(x,x0),  bas2.subs(x,x0),  bas3.subs(x,x0),  bas4.subs(x,x0) ],
		[  bas1.subs(x,x1),  bas2.subs(x,x1),  bas3.subs(x,x1),  bas4.subs(x,x1) ],
		[ dbas1.subs(x,x0), dbas2.subs(x,x0), dbas3.subs(x,x0), dbas4.subs(x,x0) ],
		[ dbas1.subs(x,x1), dbas2.subs(x,x1), dbas3.subs(x,x1), dbas4.subs(x,x1) ]])
	system = Matrix(4,1,[y0,y1,dy0,dy1])
	print ("Solving matrix ...")
	coefs  = A.LUsolve(system)
	print ("Symplyfying ...")
	for i,coef in enumerate(coefs):
		coef = coef.factor()
		coef = coef.collect([x0,x1])
		#print ( " coef %i:  " %i, coef )
		coef_subs = sympy.cse(coef)
		print ( " coef %i:  " %i, coef_subs )


#analytical_coefs_sympy()


import numpy as np
import matplotlib.pyplot as plt

def bas1(x):
	return 1 +0*x

def bas2(x):
	#return x
	return x**2

def bas3(x):
	#return x**2
	return x**4

def bas4(x):
	#return x**3
	#return x**6
	return 1/x**2

#------------ derivatives

def dbas1(x):
	return 0*x

def dbas2(x):
	#return 1 +0*x
	return 2*x

def dbas3(x):
	#return 2*x
	return 4*x**3

def dbas4(x):
	#return 3*x**2
	#return 6*x**5
	return -2/x**3

#============ Functions

def x2m(x):
	return 1/x**2

def m2x(m):
	return np.sqrt(1/m)

def getXs( x0, x1, n=10 ):
	m0 = x2m(x0);
	m1 = x2m(x1);
	ms = np.linspace(m0,m1,n+1); print( ms )
	xs = m2x(ms)
	return xs, ms[0]-ms[1], ms[-1] 

def numerical_coefs( xs, ys, dys, bass=[bas1,bas2,bas3,bas4], dbass=[dbas1,dbas2,dbas3,dbas4]):
	bas1s  =  bass[0](xs);  bas2s  = bass[1](xs);  bas3s =  bass[2](xs);  bas4s = bass[3](xs)
	dbas1s = dbass[0](xs); dbas2s = dbass[1](xs); dbas3s = dbass[2](xs); dbas4s = bass[3](xs)
	#print (bas1s)
	M      = np.zeros((4,4))
	b      = np.zeros(4)
	nseg   = len(xs)-1
	coefs  = np.zeros((nseg,4))
	for i in range(nseg):
		M[0,0]= bas1s[i  ]; M[0,1]= bas2s[i  ]; M[0,2]= bas3s[i  ]; M[0,3]= bas4s[i  ];   b[0]= ys[i  ];
		M[1,0]= bas1s[i+1]; M[1,1]= bas2s[i+1]; M[1,2]= bas3s[i+1]; M[1,3]= bas4s[i+1];   b[1]= ys[i+1];
		M[2,0]=dbas1s[i  ]; M[2,1]=dbas2s[i  ]; M[2,2]=dbas3s[i  ]; M[2,3]=dbas4s[i  ];   b[2]=dys[i  ];
		M[3,0]=dbas1s[i+1]; M[3,1]=dbas2s[i+1]; M[3,2]=dbas3s[i+1]; M[3,3]=dbas4s[i+1];   b[3]=dys[i+1]; 
		#M=M.transpose()
		coefs[i,:] = np.linalg.solve(M, b)
	return coefs

def eval_func( xs, m0, dm, coefs, bass=[bas1,bas2,bas3,bas4] ):
	#inds  = np.floor( xs/dx )
	ms    = x2m(xs);        #print ("ms:", ms )
	ms_   = (ms-m0)/dm;     #print ("ms_", ms_)
	inds  = len(coefs)-1-ms_.astype(int); #print (inds)
	#for ii in range(len(inds)):
	#	i = inds[ii]
	#	x = xs[ii]
	#	#print(  "%i %f %f %f %f" %(i,coefs[i,0],coefs[i,1],coefs[i,2],coefs[i,3]) )
	#	print(  "%i %f %f %f %f" %(i,bass[0](x),bass[1](x),bass[2](x),bass[3](x)) )
	#print ( ms_  )
	#print ( inds )
	coefi = coefs[inds,:]
	ys    = coefi[:,0]*bass[0](xs) + coefi[:,1]*bass[1](xs) + coefi[:,2]*bass[2](xs) + coefi[:,3]*bass[3](xs)
	return ys

sc = 10000
C6 = 1.0    * sc
C12= 5000.0 * sc

'''
F    =  r/|r|    * f(r)
F    =  r/|r|**2 * g(r)
g(r) = |r|*f(r)
'''

def func(x):
	#return C12/x**12 - C6/x**6         # LJ potential
	#return -12*C12/x**13 + 6*C6/x**7   # LJ force
	return -12*C12/x**12 + 6*C6/x**6    # LJ force * x	

def dfunc(x):
	#return -12*C12/x**13    + 6*C6/x**7  # of LJ potential
	#return  12*13*C12/x**14 - 6*7*C6/x**8  # of LJ force
	return   12*13*C12/x**14 - 6*7*C6/x**8  # of LJ force * x

#xs = np.array([2.0,1.0])
#xs = np.array([2.0,3.0])
#xs = np.array([4.0,4.2])
#xs = np.linspace(1,10,16+1)
xs,dm,m0=getXs( 4, 10, n=16 )
print("m0 %f dm %f " %(m0,dm) )

#x0 = xs[0]; dx = xs[1]-xs[0]
#print ("x0=%f dx=%f" %(x0,dx))
#print ("xs=",xs)

ys  =  func(xs)
dys = dfunc(xs)

coefs=numerical_coefs( xs, ys, dys )
#coefs=numerical_coefs_( xs, ys, dys )

np.savetxt("splineR2_coef.dat",coefs,comments='',header="%i %9.9e %9.9e" %(len(coefs), dm,m0))

#print( coefs )

xs_test  = np.linspace(4+0.001,10-0.001,1000)

y_test   = func(xs_test)
y_model  = eval_func( xs_test, m0, dm, coefs )
#y_model = eval_func_( xs_test, coefs )

plt.plot(xs, ys, 'o')
plt.plot(xs_test, y_test)
plt.plot(xs_test, y_model)

plt.figure()
plt.plot(xs_test, (y_model-y_test))

#np.savetxt("spline_debug.dat", np.array([xs_test,y_model]).transpose() )

plt.show()









