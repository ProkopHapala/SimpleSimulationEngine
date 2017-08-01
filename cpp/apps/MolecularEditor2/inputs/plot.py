#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt

for elname in ['H','Xe']:
    data=np.transpose( np.genfromtxt('force_%s.dat' %elname) )
    #print data
    plt.plot(data[3],data[6], label=elname )

#vmin=data[6].min()
vmin=-0.03
plt.ylim(vmin*1.2,-vmin*1.2)
plt.grid()
plt.legend()
plt.show()


