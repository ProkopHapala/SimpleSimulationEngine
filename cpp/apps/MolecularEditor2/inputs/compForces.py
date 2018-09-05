#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt( "Fgrid.log" )

for i,d in enumerate(data[0]) : print i," : ",d

plt.plot( data[:, 5] ,data[:,9 ], "r-", label="F_CPU" )
plt.plot( data[:,14] ,data[:,18], "b:", label="F_GPU" )

#plt.plot( data[:, 5-1] ,data[:,9 ], "r-", label="F_CPU" )
#plt.plot( data[:,14-1] ,data[:,18], "b:", label="F_GPU" )

#plt.plot( data[:, 5-2] ,data[:,9 ], "r-", label="F_CPU" )
#plt.plot( data[:,14-2] ,data[:,18], "b:", label="F_GPU" )

plt.legend()
plt.show()



