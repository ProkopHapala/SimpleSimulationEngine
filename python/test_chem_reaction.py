#!/usr/bin/python

import re
import numpy as np
import sys
from pySimE import chemistry as ch

#print ch.str2composition( sys.argv[1] )
#sides =  ch.parseReaction( 'Fe+O2=Fe2O3' )
#sides =  ch.parseReaction( 'C12H22O11+KNO3=H2O+CO2+K2CO3+N2' )
#print sides
#print ch.reaction2string( sides )
#print ch.balanceReactionString( 'Fe+O2=Fe2O3' )

print ch.balanceReactionString( 'C12H22O11+KNO3=H2O+CO2+K2CO3+N2' )

#print atomicBalance( reaction[0], reaction[1] )


