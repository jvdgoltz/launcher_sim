#import all necessary funtions
from isatmos import dens, temp, press
from math import log, sqrt, atan2, sin, cos, pi				#TODO: if possible, import functions to do integrations
from numpy import genfromtxt
import matplotlib.pyplot as plt

#define variables and read out data source
table = genfromtxt("V2.dat", usecols=(1), dtype="float", comments="#", delimiter="=",)

oem=table[1]
fm=table[2]
Isp=table[3]
mflw=table[4]

x=0
y=0
t=0
v=0

g0=9.81
m0=oem+fm
m=m0

dt=input("Please enter the resolution of the simulation in seconds.")

ttab=[]
mtab=[]
vtab=[]
ytab=[]

while m>oem:
    t = t + dt
    ttab.append(t)
    
    m = m - dt * mflw
    mtab.append(m)
    
    v = log(m0/m)*Isp*g0
    vtab.append(v)

    y=y+v*dt
    ytab.append(y)
 
plt.plot(ttab, vtab)
plt.plot(ttab, mtab)
plt.plot(ttab, ytab)
plt.show()
