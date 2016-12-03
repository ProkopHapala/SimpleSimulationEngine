#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import pyMolecular as mol

psamples = np.linspace(1.0, 10.0, 100 )[:,None] * np.array([1.0,1.0,1.0])[None,:]

ps = np.random.random((10,3)) - 0.5
Qs = np.random.random(len(ps)) - 0.5

Eaprox, Eref = mol.testMultipole( ps, Qs, psamples, order=2 )

plt.plot( psamples[:,0], Eref ,  label="exact");
plt.plot( psamples[:,0], Eaprox, label="Quadrupole" );

Eaprox, Eref = mol.testMultipole( ps, Qs, psamples, order=1 )
plt.plot( psamples[:,0], Eaprox, label="Dipole" );

Eaprox, Eref = mol.testMultipole( ps, Qs, psamples, order=0 )
plt.plot( psamples[:,0], Eaprox, label="Monopole" );

plt.legend()
plt.show()
