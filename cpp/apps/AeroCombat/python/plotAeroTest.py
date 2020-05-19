# -*- coding: utf-8 -*-
"""
Created on Mon Apr 30 04:32:39 2018

@author: prokop.hapala
"""

import numpy as np
import matplotlib.pyplot as plt


data = np.transpose( np.genfromtxt( "O:\\vbs3_beta\\AeroTest_1.log" ) )
print data

elevator = data[1]
Lift =  data[6];  plt.plot( elevator, Lift, label="Lift"  )
Drag = -data[7];  plt.plot( elevator, Drag, label="Drag"  )
AoA  = np.arctan2( data[9], data[10] );  plt.plot( elevator, AoA*180/np.pi, label="AoA"  )

plt.legend()
plt.grid()

plt.show( )