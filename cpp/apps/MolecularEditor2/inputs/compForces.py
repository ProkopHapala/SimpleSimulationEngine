#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt( "Fgrid.log" )

for i,d in enumerate(data[0]) : print i," : ",d

plt.figure(figsize=(5,15))

i1=535;
i2=630;
xcpu = data[:, 5]
ycpu = data[:, 9]
yl = (xcpu-xcpu[i1]) * (  (ycpu[i2]-ycpu[i1])/(xcpu[i2]-xcpu[i1])  ) + ycpu[i1]


plt.subplot(3,1,1)
plt.plot( data[:, 5] ,data[:,9 ], "r-", label="F_CPU" )
plt.plot( data[:,14] ,data[:,18], "b:", label="F_GPU" )
plt.legend()

plt.subplot(3,1,2)
plt.plot( data[:, 5] ,data[:,9 ]-yl, "r-", lw=0.2, label="F_CPU" )
plt.plot( data[:,14] ,data[:,18]-yl, "b-", lw=0.2, label="F_GPU" )
#plt.plot( xcpu ,yl, "g:", label="yl" )
#plt.plot( xcpu[[i1,i2]] ,yl[[i1,i2]], "go", label="yl" )
#vmin = data[:,9].min()*1.2
#print "vmin : ", vmin, -vmin
vmin = 0.0005
plt.ylim( vmin, -vmin )
plt.legend()

'''
plt.subplot(3,1,2)
#plt.plot( data[:, 5] ,data[:,9 ]-data[:,18], "r-", label="F_CPU-F_GPU" )
plt.plot( data[:, 5] ,(data[:,9 ]-data[:,18])/np.abs(data[:,18]), "g-", label="rel. err." )
plt.legend()
'''

plt.subplot(3,1,3)
plt.plot( data[:, 5] ,data[:,9 ]-data[:,18], "r-", label="F_CPU-F_GPU" )
#plt.plot( data[:, 5] ,(data[:,9 ]-data[:,18])/np.abs(data[:,18]), "g-", label="rel. err." )
plt.legend()

#plt.plot( data[:, 5-1] ,data[:,9 ], "r-", label="F_CPU" )
#plt.plot( data[:,14-1] ,data[:,18], "b:", label="F_GPU" )

#plt.plot( data[:, 5-2] ,data[:,9 ], "r-", label="F_CPU" )
#plt.plot( data[:,14-2] ,data[:,18], "b:", label="F_GPU" )

plt.show()
