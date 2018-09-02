#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  1 20:16:12 2018

@author: jvdgoltz
"""

import modules.trajectory as traj
import modules.guidance as guid
import matplotlib.pyplot as plt
import numpy as np

dt = 0.25
t_max = 12*60

t,m,x,y,T,vx,vy,q,th,ax,ay,M,D = traj.trajectory(traj.vehicle,traj.state,0,dt,t_max)


aconstr = T/m/np.sqrt(2)
x0 = np.array([0,0,0,0,0,-9.81])
targetalt = 200*1e3
targetvel = 7.9*1e3
#t1, X = guid.solve(x0,targetalt,targetvel,aconstr,dt,min(max(t),t_max))

plt.figure()
plt.plot(t,X[:,1])
plt.figure()
plt.plot(X[:,0],X[:,1])
plt.figure()
plt.plot(t,X[:,6])
plt.plot(t,X[:,7])
plt.plot(t,np.sqrt(X[:,6]**2+X[:,7]**2))
plt.plot(t,aconstr)
plt.legend(['x-acceleration','y-acceleration','total acceleration','max acceleration constraint'])


print('Downrange:',X[-1,0]/1000,'km')
print('Altitude:',X[-1,1]/1000,'km')
print('Horizontal Velocity:',X[-1,2]/1000,'km/s')
print('Vertical Velocity:',X[-1,3]/1000,'km/s')