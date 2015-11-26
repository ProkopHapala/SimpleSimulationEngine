#!/usr/bin/python

from pylab import *

# from  http://hyperphysics.phy-astr.gsu.edu/hbase/relativ/relmom.html
#http://hyperphysics.phy-astr.gsu.edu/hbase/relativ/releng.html#c5
# pc = sqrt( E**2 - m0**2 c**4  )
# v/c = pc/E
# Ek = E - E0
# v = sqrt( E**2 - m0**2 c**4  ) * ( c / E )
# v = sqrt( 1 - (m0*c**2/E)**2 )  * c
# v = sqrt( 1 - ( m0*c**2/( m0*c**2 + Ek ) )**2 ) * c
# v = sqrt( ( ( m0*c**2 + Ek )**2 -  (m0*c**2)**2  ) / ( m0*c**2 + Ek )**2 ) * c
# v = sqrt( (   2*Ek*m0*c**2 + Ek**2  ) / ( m0*c**2 + Ek )**2 ) * c
# v = sqrt( Ek*(   2*m0*c**2 + Ek  ) / ( m0*c**2 + Ek )**2 ) * c
# v = sqrt( Ek*(   2*E0 + Ek  ) / ( E0 + Ek )**2 ) * c

# pc = sqrt(  E**2 - E0**2               )
# v  = sqrt( (E0+Ek)**2 - E0**2 ) / (E0+Ek)**2 ) * c
# v  = sqrt( 2*E0*Ek + Ek**2 ) / (E0+Ek)**2 ) * c

SPEED_OF_LIGHT = 299792458.0
eVtoSI         = 1.60217657e-19



def classical_Energy_to_Velocity( Ek, m0 = 1.67262178e-27 ):
	return sqrt( 2*Ek/m0 )

def relativistic_Energy_to_Velocity( Ek, m0 = 1.67262178e-27 ):
	E0 = m0 * SPEED_OF_LIGHT**2
	return SPEED_OF_LIGHT * sqrt( Ek*(2*E0 + Ek )/(E0+Ek)**2 )


#
Es       = 10**linspace( -3, 10, 14 ) 
vs_rel   = relativistic_Energy_to_Velocity( Es*eVtoSI  )
vs_class = classical_Energy_to_Velocity   ( Es*eVtoSI  )

for i in range( len(Es) ):
	print " %2.2e eV %2.2e m/s ( %2.2e c ) %2.2e m/s ( %2.2e c ) " %(Es[i], vs_class[i], vs_class[i]/SPEED_OF_LIGHT, vs_rel[i], vs_rel[i]/SPEED_OF_LIGHT )


Es      = 10**linspace( -3, 10, 100 ) 
vs_rel   = relativistic_Energy_to_Velocity( Es*eVtoSI  )
vs_class = classical_Energy_to_Velocity   ( Es*eVtoSI  )

plot( vs_rel/SPEED_OF_LIGHT                 , Es, label='relativistic' )
plot( vs_class/SPEED_OF_LIGHT, Es, label='classical'    )

yscale('log')
xscale('log')
xlabel(" v/c ")
ylabel(" Ek [ MeV ] ")
legend(loc=2)

grid()

show()
