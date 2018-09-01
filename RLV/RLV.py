# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 20:14:02 2014

@author: Julian

calculating the physics of a reusable rocket with controlled powered return capability

Initial Vehicle input is retrieved from reqs.dat
"""
import numpy as np
import matplotlib.pyplot as plt

input_data = np.genfromtxt('reqs.dat', dtype='float', comments='#', delimiter=' ')

I_sp = input_data[0]
PL = input_data[1]
DV = input_data[2]
lamb = PL/(505846-PL)

c = I_sp * 9.81

MR = np.exp(DV/c)

eps = (1+lamb+lamb*MR)/MR

#rocket mass ratios
M_s = PL * eps/lamb
M_e = M_s + PL
M_0 = MR * M_e
M_f = M_0 - M_e
M_s, M_e, M_0, M_f = (np.round((M_s, M_e, M_0, M_f), decimals=0))

print
print 'Required rocket data (not reusable):'
print
print 'Payload Mass:        ',  PL,    '[kg]'
print 'Structural Mass:     ', M_s,    '[kg]'
print 'Empty Mass:          ', M_e,    '[kg]'
print 'Fuel Mass:           ', M_f,    '[kg]'
print 'Lift-Off Mass:       ', M_0,    '[kg]'

MR = np.exp(DV/c)

M_sr = M_s/ ((MR*(1-1/MR))*(1+1/(MR*(1-1/MR))))
M_e = M_sr + PL
M_fr = (MR+1)*M_sr + M_f
M_0 = M_e + M_sr + M_fr
M_sr, M_e, M_0, M_fr = (np.round((M_sr, M_e, M_0, M_fr), decimals=0))

print
print 'Required rocket data (reusable):'
print
print 'Payload Mass:        ',  PL,    '[kg]'
print 'Structural Mass:     ', M_sr,    '[kg]'
print 'Empty Mass:          ', M_e,    '[kg]'
print 'Fuel Mass:           ', M_fr,    '[kg]'
print 'Lift-Off Mass:       ', M_0,    '[kg]'

print
print (PL/lamb+PL)*np.exp(-2*DV/c)