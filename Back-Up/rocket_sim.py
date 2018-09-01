#this programm simulates the flight of multistage rocket from lift-off to orbit insertion
#input parameters are:
	#the empty mass of the rocket
	#the fuel mass of the rocket
	#the engine expansion ratio, exit area, chamber pressure, propellant molecular mass
	#the mass flow of the rocket engine
	#the drag coefficient of the rocket
	#the drag surface of the rocket
#it is suggested that the input variables are taken out from a file that will be loaded in the beginning

#the effects of wind and turbulences are neglected, as well as the rounded character of the earth and celestial effects

#the result of the simulation are the graphs:
	#y-position vs x-position
	#y-position vs time
	#velocity vs time
	#acceleration vs time
	#dynamic pressure vs time
     #Mach number vs time
	
#import all necessary funtions
from isatmos import dens
from math import sqrt, atan2, sin, cos, pi				
from xlrd import *
import matplotlib.pyplot as plt

#read out data source

#filename=raw_input("Please enter the data source file name.")
filename="Saturn5.xls"
table = open_workbook(filename)
sheet = table.sheet_by_name('Sheet')

#define rocket technical data lists
oem = []
fm = []
Isp = []
mflw = []
cd = []
S = []
c_alt = []
sepdur = []

dt = sheet.cell_value(3,1)
maxt = sheet.cell_value(4,1)
nstages = int(sheet.cell_value(1,1))
payld = sheet.cell_value(2,1)
targetalt = sheet.cell_value(14,1)
targetvel = sheet.cell_value(15,1)
for i in range(nstages):
    oem.append(sheet.cell_value(6,1+i))
    fm.append(sheet.cell_value(7,1+i))
    Isp.append(sheet.cell_value(8,1+i))
    mflw.append(sheet.cell_value(9,1+i))
    cd.append(sheet.cell_value(10,1+i))
    S.append(sheet.cell_value(11,1+i))
    c_alt.append(sheet.cell_value(12,1+i))
    sepdur.append(sheet.cell_value(13,i+1))

alpha=(pi*90/180)


x=0
y=0
t=0
v=0

theta=alpha
vx=v*cos(theta)
vy=v*sin(theta)
rho=1.225
g0=9.81                        
q=0.5*rho*v*v
D0=cd[0]*q*S[0]
Dx=-cd[0]*0.5*rho*vx*v*S[0]
Dy=-cd[0]*0.5*rho*vy*v*S[0]

m0=sum(oem)+sum(fm)+payld
m=m0
W=m*g0

T0=Isp[0]*mflw[0]*g0
Tx=T0*cos(theta)
Ty=T0*sin(theta)
Z=0
Rx=Dx+Tx
Ry=Dy+Ty-W
a0=sqrt(Rx**2+Ry**2)/m0

ttab=[0]
mtab=[m0]
vtab=[0]
atab=[a0]
xtab=[0]
ytab=[0]
thetab = [0]
qtab=[0]
for i in range(nstages):
    if i>0:
        fm[i-1]=0
        oem[i-1]=0
        
    while fm[i]>0 and y>=0 and Z<W:
            t = t + dt
            ttab.append(t)
            fm[i]=fm[i]-dt*mflw[i]
            m = sum(oem) + sum(fm) +payld
        
            mtab.append(m)
            
            v=sqrt(vx*vx+vy*vy)
            rho=float(dens(y))
            q=0.5*rho*v*v
            qtab.append(q)
            g=g0*(6371000/(6371000+y))**2
            W=m*g
            T=Isp[i]*mflw[i]*g           #Thrust
            Tx=T*cos(theta)
            Ty=T*sin(theta)
            
            Dx=-cd[i]*0.5*rho*vx*v*S[i]   #Drag
            Dy=-cd[i]*0.5*rho*vy*v*S[i]

            Z=m*vx*vx/(y+6371000)            
            
            Rx=Dx+Tx                #Resulting Force on the Rocket
            Ry=Dy+Ty-W+Z
            
            ax=Rx/m                 #acceleration due to resulting 
            ay=Ry/m
            a=sqrt(ay*ay+ax*ax)
            atab.append(a)
            vx=(vx+ax*dt)
            vy=(vy+ay*dt)
                      
            gamma=theta
            
            if y>=c_alt[i] and y<=targetalt:
                targetvel = sqrt(g*(6371000+targetalt))
                deltagamma = (vx-targetvel)*gamma/112000+(targetalt-y)*(pi/2-gamma)*0.9/2890260
                if deltagamma*dt+gamma >= -0 and deltagamma*dt+gamma <= pi/2:                
                    theta = deltagamma*dt+gamma
                
            if y>=targetalt and vx<=targetvel:
                targetvel = sqrt(g*(6371000+y))
                deltagamma = (vx-targetvel)*gamma/110000-vy*0.002
                if deltagamma*dt+gamma >= -0 and deltagamma*dt+gamma <= pi/2:                
                    theta = deltagamma*dt+gamma
                
            if y>=targetalt and vx>=targetvel:
                theta = -vy/abs(vy)*pi/2
            vtab.append(v)
            x=x+vx*dt
            y=y+vy*dt
            thetab.append(theta)
            xtab.append(x)
            ytab.append(y)
            tsep = ttab[-1]
    while fm[i]<=0 and t<=tsep+sepdur[i] and y>=0:  #  this loop is for the free flight phase after stage separation
            t = t + dt
            ttab.append(t)
        
            m = sum(oem) + sum(fm) + payld
        
            mtab.append(m)
            
            v=sqrt(vx*vx+vy*vy)
            rho=float(dens(y))
            q=0.5*rho*v*v
            qtab.append(q)
            
            g=g0*(6371000/(6371000+y))**2
            W=m*g
            try:
                Dx=-cd[i+1]*0.5*rho*vx*v*S[i+1]
                Dy=-cd[i+1]*0.5*rho*vy*v*S[i+1]
            except:
                Dx=-cd[i]*0.5*rho*vx*v*S[i]
                Dy=-cd[i]*0.5*rho*vy*v*S[i]
            
            Z=m*vx*vx/(y+6371000)
            
            Rx=Dx
            Ry=Dy-W+Z
            ax=Rx/m
            ay=Ry/m
            a=sqrt(ay*ay+ax*ax)
            atab.append(a)
            vx=(vx+ax*dt)
            vy=(vy+ay*dt)
            theta=(-g*cos(theta)/(sqrt(vx*vx+vy*vy)))*dt+theta
            thetab.append(theta)
            vtab.append(v)
            x=x+vx*dt
            y=y+vy*dt
            xtab.append(x)
            ytab.append(y)
plt.subplot(322)
plt.plot(xtab, ytab)
plt.title("flight profile: altitude vs ground range")

plt.subplot(321)
plt.plot(ttab, ytab)
plt.title("flight profile: altitude vs time")

plt.subplot(323)
plt.plot(ttab, vtab)
plt.title("velocity vs time")

plt.subplot(325)
plt.plot(ttab, atab)
plt.title("acceleration vs time")

plt.subplot(324)
plt.plot(ttab, qtab)
plt.title("dynamic pressure vs time")

plt.subplot(326)
plt.plot(ttab, thetab) #change!
plt.title("rocket mass vs time")

print targetvel
print vx
print vy
print y

plt.show()

print "done"

