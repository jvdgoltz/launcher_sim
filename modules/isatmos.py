#with this python module it is possible to calculate temperature, density and pressure of the earth atmosphere at a certain altitude.
#The base of these calculations is the International Standard Atmosphere.
#The functions are called "temp(altitude)", "dens(altitude)" and "press(altitude)"
#It is possible to use them in other modules by importing them ("from isatmos import dens, temp, press")

#made by Julian von der Goltz
#last edit: 24-04-2013 10:15 by Julian


#import mathematical functions
from math import e

#define variables

t0=(273.15+15)
p0=101325
rho0=1.225

R=287.05

g=9.81

a0=-0.0065
a1=0
a2=0.001
a3=0.0028
a4=0
a5=-0.0028
a6=-0.002                                               #it's more efficient to save these values in a list!
                                                        #same for rho0,rho1...and t0,t1,...

def temp(h):
    
    if h>=0 and h<=11000:                               #L0 troposphere
        t=t0+a0*h                           
        
    elif h>11000 and h<=20000:                          #L1 tropopause isotherm
        t=t0+a0*11000
        
    elif h>20000 and h<=32000:                          #L2 stratosphere
        t2=t0+a0*11000                      
        t=t2+a2*(h-20000)
        
    elif h>32000 and h<=47000:                          #L3 stratosphere
        t3=t0+a0*11000+a2*(12000)
        t=t3+a3*(h-32000)

    elif h>47000 and h<=51000:                          #L4 stratopause isotherm
        t=t0+a0*11000+a2*(12000)+a3*(15000)

    elif h>51000 and h<=71000:                          #L4 mesosphere
        t5=t0+a0*11000+a2*(12000)+a3*(15000)
        t=t5+a5*(h-51000)

    elif h>71000 and h<=84854:                          #L5 mesosphere
        t6=t0+a0*11000+a2*(12000)+a3*(15000)+a5*(20000)
        t=t6+a6*(h-71000)

    else:
        t=0
    return t

def press(h):   
    if h>=0 and h<=84854:
        p=dens(h)*R*temp(h)                                 #Ideal Gas law to calculate pressure with density and temperature
    else:
        p=0
    return p

def dens(h):

    rho1=rho0*(temp(11000)/temp(00000))**(-g/(R*a0)-1)      #density at the beginning of each layer, respectively
    rho2=rho1*e**(-g/(R*temp(20000))*(20000-11000))
    rho3=rho2*(temp(32000)/temp(20000))**(-g/(R*a2)-1)
    rho4=rho3*(temp(47000)/temp(32000))**(-g/(R*a3)-1)
    rho5=rho4*e**(-g/(R*temp(51000))*(51000-47000))
    rho6=rho5*(temp(71000)/temp(51000))**(-g/(R*a5)-1)      

    if h>=0 and h<=11000:  
        rho=rho0*(temp(h)/temp(00000))**(-g/(R*a0)-1)
        
    elif h>11000 and h<=20000:                              #L1 isotherm!
        rho=rho1*e**(-g/(R*temp(h))*(h-11000))
        
    elif h>20000 and h<=32000:              
        rho=rho2*(temp(h)/temp(20000))**(-g/(R*a2)-1)
        
    elif h>32000 and h<=47000:                   #
        rho=rho3*(temp(h)/temp(32000))**(-g/(R*a3)-1)

    elif h>47000 and h<=51000:                              #L4 isotherm!
        rho=rho4*e**(-g/(R*temp(h))*(h-47000))

    elif h>51000 and h<=71000:                     #
        rho=rho5*(temp(h)/temp(51000))**(-g/(R*a5)-1)

    elif h>71000 and h<=84854:                #
        rho=rho6*(temp(h)/temp(71000))**(-g/(R*a6)-1)

    else:
        rho=0
        
    return rho

def atmos(h):
    prt=[press(h), dens(h), temp(h)]
    return prt
