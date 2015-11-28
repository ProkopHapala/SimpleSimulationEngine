#!/usr/bin/python

from pylab import *

import pySimE as sim

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

#SPEED_OF_LIGHT = 299792458.0
#eVtoSI         = 1.60217657e-19

#
Es       = 10**linspace( -3, 10, 14 ) 
vs_rel   = sim.relativistic_Energy_to_Velocity( Es*sim.const.eV  )
vs_class = sim.classical_Energy_to_Velocity   ( Es*sim.const.eV  )

for i in range( len(Es) ):
	print " %2.2e eV %2.2e m/s ( %2.2e c ) %2.2e m/s ( %2.2e c ) " %(Es[i], vs_class[i], vs_class[i]/sim.const.lightSpeed, vs_rel[i], vs_rel[i]/sim.const.lightSpeed )


Es      = 10**linspace( -3, 10, 100 ) 
vs_rel   = sim.relativistic_Energy_to_Velocity( Es*sim.const.eV )
vs_class = sim.classical_Energy_to_Velocity   ( Es*sim.const.eV )

plot( vs_rel/sim.const.lightSpeed                 , Es, label='relativistic' )
plot( vs_class/sim.const.lightSpeed, Es, label='classical'    )

yscale('log')
xscale('log')
xlabel(" v/c ")
ylabel(" Ek [ MeV ] ")
legend(loc=2)

grid()

show()
