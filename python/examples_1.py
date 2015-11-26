

import pySimE as sim
from   pylab import *

Ts = np.linspace( 300,3000 )
Ps = sim.thermalRadiativePower( temperature = Ts )
plot( Ps )

show()


'''
print " Compute equlibrium temperature of black 1-side surface at Earth orbit : "
P_1AU  = ks.solarSystem.Sun['Power'] / (4*pi*ks.const_AU**2)
print " Power density : ",P_1AU," [W] "
T_1AU  = ks.equilibriumTemperature(P_1AU) - ks.const_0celsius
print " Temperature ( 1 side )   : ",T_1AU," [°C] "
T_1AU_ = ks.equilibriumTemperature(P_1AU, absorptivity=0.5, emissivity=2.0 ) - ks.const_0celsius
print " Temperature ( 2 side solar panel)  : ",T_1AU_," [°C] "
'''
