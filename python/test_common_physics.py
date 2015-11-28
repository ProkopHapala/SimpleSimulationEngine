#!/usr/bin/python

import pySimE as sim
from   pylab import *

Ts = np.linspace( 300,3000 )
Ps = sim.thermalRadiativePower( temperature = Ts )
plot( Ps )

show()
