# -*- coding: utf-8 -*-
"""
Created on Sat Jul 20 23:42:30 2013

@author: asiJa
"""

import pylab

def railGun_Velocity( length=10, acceleration=9.81*250000 ):
	return pylab.sqrt(2*length*acceleration)

def railGun_length( velocity, acceleration=9.81*250000 ):	
	return velocity**2/(2*acceleration)







if __name__ == "__main__":	
	vs = pylab.arange( 1,12e+3,1e+3 )
	# projectile gun
	ls = railGun_length( vs,  acceleration=9.81*250000 )
	pylab.plot(vs/1000,ls,'m-o'); pylab.xlabel('v [km/s]'); pylab.ylabel('L [m]'); 
	#space launch - human crew
	pylab.figure()
	ls = railGun_length( vs,  acceleration=9.81*10 )
	pylab.plot(vs/1000,ls/1000,'m-o'); pylab.xlabel('v [km/s]'); pylab.ylabel('L [km]');   
	#pylab.yscale('log')
	pylab.show()