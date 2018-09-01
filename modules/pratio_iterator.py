# -*- coding: utf-8 -*-
"""
Created on Fri May 02 16:04:46 2014

@author: Julian
"""
from thrust import vdk
kappa = 1.3
epsilon = 49
pratio = [0]
iterations = range(100)
Gamma = vdk(kappa)
f = (Gamma**2/((2*kappa/(kappa-1))*epsilon**2) + pratio[-1]**((1+kappa)/kappa))**(kappa/2) - pratio[-1]
it = 0
while f<>0 and it<100:
    f = (Gamma**2/((2*kappa/(kappa-1))*epsilon**2) + pratio[-1]**((1+kappa)/kappa))**(kappa/2) - pratio[-1]
    fp = ((1+kappa)/kappa) * pratio[-1]**((1+kappa)/kappa -1) * (kappa/2) * (Gamma**2/((2*kappa/(kappa-1))*epsilon**2) + pratio[-1]**((1+kappa)/kappa))**(kappa/2-1) - 1

    pratio.append(pratio[-1]-f/fp)
    it = it +1
print 1/pratio[-1]

import matplotlib.pyplot as plt
plt.plot(pratio)