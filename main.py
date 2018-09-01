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
from math import sqrt, atan2, sin, cos, pi				
from xlrd import *
import numpy as np
import matplotlib.pyplot as plt
import aerodynamics as aero
import thrust as th
import controls as ctl
from isatmos import press

#read out data source

#filename=raw_input("Please enter the data source file name.")
filename="Saturn5.xls"
table = open_workbook(filename)
sheet = table.sheet_by_name('Sheet')

#define rocket technical data lists
oem = []
fm = []
Tc = []
pc = []
epsilon = []
Mw = []
kappa = []
Ae = []
nengines = []
cd = []
S = []
t_ctl = []
sepdur = []

dt = sheet.cell_value(3,1)
maxt = sheet.cell_value(4,1)
nstages = int(sheet.cell_value(1,1))
payld = sheet.cell_value(2,1)
targetalt = sheet.cell_value(5,1)
targetvel = sheet.cell_value(6,1)
for i in range(nstages):
    oem.append(sheet.cell_value(8,1+i))
    fm.append(sheet.cell_value(9,1+i))
    Tc.append(sheet.cell_value(10,1+i))
    pc.append(sheet.cell_value(11,1+i))
    epsilon.append(sheet.cell_value(12,1+i))
    Mw.append(sheet.cell_value(13,1+i))
    kappa.append(sheet.cell_value(14,1+i))
    Ae.append(sheet.cell_value(15,1+i))
    nengines.append(sheet.cell_value(16,1+i))
    cd.append(sheet.cell_value(17,1+i))
    S.append(sheet.cell_value(18,1+i))
    t_ctl.append(sheet.cell_value(19,1+i))
    sepdur.append(sheet.cell_value(20,1+i))
    

At = np.divide(Ae,epsilon)
pratio = nstages * [0]
for p in range(nstages):
    pratio[p] = th.pratio(epsilon[p],kappa[p])

gamma=(pi/2)        #flight path / pitch angle
omega=0             #angular velocity
alpha=0             #angular acceleration
xmax = 1000000.

x=0
y=0
t=0
v=0.0001


vx=v*cos(gamma)
vy=v*sin(gamma)
rho=1.225
p0=101325.
g0=9.81                        
q=0.5*rho*v*v
D0=cd[0]*q*S[0]
Dx=-cd[0]*0.5*rho*vx*v*S[0]
Dy=-cd[0]*0.5*rho*vy*v*S[0]

m0=sum(oem)+sum(fm)+payld
m=m0
W=m*g0

mflw = th.mdot(At[0],pc[0],Tc[0],Mw[0],kappa[0])
ve = th.vexit(epsilon[0], Tc[0], Mw[0], kappa[0])
T0 = nengines[0]*th.FT(mflw,ve,At[0],epsilon[0],pc[0],pratio[0],p0)
Tx=T0*cos(gamma)
Ty=T0*sin(gamma)
Z=0
Rx=Dx+Tx
Ry=Dy+Ty-W
ax = Rx/m0
ay = Ry/m0
a0=sqrt(Rx**2+Ry**2)/m0

ttab=[0]
mtab=[m0]
Ttab=[T0]
Drag=[0]
vtab=[0]
atab=[a0]
xtab=[0]
ytab=[0]
y_ref=[0]
gammab = [0]
ga_ref = [0]
qtab=[0]
Mtab=[0]
startstage = []
for i in range(nstages):
    if i>0:
        fm[i-1]=0
        oem[i-1]=0
    startstage.append(t)  
    while fm[i]>0 and y>=0 and t<2000:
            t = t + dt
            ttab.append(t)
            if t>135 and i==0:
                nen = nengines[i]-1
            elif t>400 and i==1:
                nen = nengines[i]-1
            else:
                nen = nengines[i]
            mflw = nen*th.mdot(At[i],pc[i],Tc[i],Mw[i],kappa[i])
            fm[i]=fm[i]-dt*mflw
            m = sum(oem) + sum(fm) +payld
            mtab.append(m)


            Tk  = [Ttab[-1]]
            axk = [ax]
            ayk = [ay]
            vxk = [vx]
            vyk = [vy]
            xk  = [x]
            yk  = [y]
            
            for k in range(3):
                             
                g=g0*(6371000/(6371000+yk[-1]))**2
                W=m*g
            
                
                ve = th.vexit(epsilon[i], Tc[i], Mw[i], kappa[i])
                pa = press(yk[-1])
                
                if Z<W:
                    T = th.FT(mflw,ve,At[i],epsilon[i],pc[i],pratio[i],pa)
                else:
                    T = 0
                Tk.append(T)
                Tx=T*cos(gamma)
                Ty=T*sin(gamma)                 
                
                aerodyn = aero.drag(yk[-1],vxk[-1],vyk[-1],cd[i],S[i])
                Dx = aerodyn[0]
                Dy = aerodyn[1]
                M  = aerodyn[2]
                q  = aerodyn[3]
                
            
                Z=m*vxk[-1]*vxk[-1]/(yk[-1]+6371000)            
            
                Rx=Dx+Tx                #Resulting Force on the Rocket
                Ry=Dy+Ty-W+Z
            
                axk.append(Rx/m)                 #acceleration due to Resultant
                ayk.append(Ry/m)
                
                
                vxk.append(vx+axk[-1]*dt)
                vyk.append(vy+ayk[-1]*dt)
                
                xk.append(x+vxk[-1]*dt)
                yk.append(y+vyk[-1]*dt)

            T = Tk[-1]
            Ttab.append(T)
            ax = axk[-1]
            ay = ayk[-1]
            a=sqrt(ay*ay+ax*ax)
            atab.append(a)
            
            Drag.append(np.sqrt(Dx*Dx+Dy*Dy))            
            
            vx = vxk[-1]
            vy = vyk[-1]
            v=sqrt(vx*vx+vy*vy)
            vtab.append(v)
            
            x  =  xk[-1]
            y  =  yk[-1]
            xtab.append(x)
            ytab.append(y)
            
            Mtab.append(M)
            qtab.append(q)            
            y_ref.append(ctl.reference(x,y,xmax,targetalt)[1])
            ga_ref.append(ctl.reference(x,y,xmax,targetalt)[0])
            alpha = ctl.pitchrate(x,y,vx,vy,gamma,omega,targetvel,xmax,targetalt,t_ctl[i],W,Z,T,Dx,sum(fm),mflw,m,t,11*60+40)
            omega = omega + alpha*dt
            gamma = gamma + omega*dt
            
            gammab.append(gamma)
            
            tsep = ttab[-1]

    while fm[i]<=0  and y>=0 and t<=tsep+sepdur[i]:  #  this loop is for the free flight phase after stage separation
            t = t + dt
            ttab.append(t)
        
            m = sum(oem) + sum(fm) + payld
        
            mtab.append(m)
            
            
            
            g=g0*(6371000/(6371000+y))**2
            W=m*g
            try:
                aerodyn = aero.drag(y,vx,vy,cd[i+1],S[i+1])
            except:
                aerodyn = aero.drag(y,vx,vy,cd[i],S[i])
            Dx = aerodyn[0]
            Dy = aerodyn[1]
            M  = aerodyn[2]
            q  = aerodyn[3]
            Mtab.append(M)
            qtab.append(q)
            
            Z=m*vx*vx/(y+6371000)
            
            Rx=Dx
            Ry=Dy-W+Z
            Ttab.append(0)
            ax=Rx/m
            ay=Ry/m
            a=sqrt(ay*ay+ax*ax)
            atab.append(a)
            
            Drag.append(np.sqrt(Dx*Dx+Dy*Dy))                
            
            vx=(vx+ax*dt)
            vy=(vy+ay*dt)
            
            vtab.append(v)
            x=x+vx*dt
            y=y+vy*dt
            
            xtab.append(x)
            ytab.append(y)
            y_ref.append(ctl.reference(x,y,xmax,targetalt)[1])
            ga_ref.append(ctl.reference(x,y,xmax,targetalt)[0])
            alpha = ctl.pitchrate(x,y,vx,vy,gamma,omega,targetvel,xmax,targetalt,t_ctl[i],W,Z,T,Dx,sum(fm),mflw,m,t,11*60+40)
            omega = omega + alpha*dt
            gamma = gamma + omega*dt
            gammab.append(gamma)
ttab = np.array(ttab)/60
plt.subplot(331)
plt.plot(ttab, np.array(ytab)/1000)
plt.title("flight profile: altitude [km] vs time [min]")

plt.subplot(332)
plt.plot(np.array(xtab)/1000, np.array(ytab)/1000, 'b')
plt.plot(np.array(xtab)/1000, np.array(y_ref)/1000, 'r')
plt.title("flight profile: altitude [km] vs ground range [km]")

plt.subplot(338)
plt.plot(ttab, np.array(Ttab)/1000000)
plt.title("Thrust [MN] vs time [min]")

plt.subplot(334)
plt.plot(ttab, np.array(vtab)/1000)
plt.title("velocity [km/s] vs time [min]")

plt.subplot(336)
plt.plot(ttab, np.array(qtab)/1000)
plt.title("dynamic pressure [kPa] vs time [min]")

plt.subplot(335)
plt.plot(ttab, np.array(gammab)*180/np.pi,'b')
plt.plot(ttab, np.array(ga_ref)*180/np.pi,'r')
plt.title("pitch angle [deg] vs time [min]")

plt.subplot(337)
plt.plot(ttab, np.array(atab)/g0+1)
plt.title("acceleration [g0] vs time [min]")

plt.subplot(333)
plt.plot(ttab, np.array(Mtab)) 
plt.title("Mach number [-] vs time [min]")

plt.subplot(339)
plt.plot(ttab, np.array(Drag)/1000) 
plt.title("Drag [kN] vs time [min]")

print 'Orbit Insertion Velocity: ', round(targetvel), 'm/s'
print 'Horizontal Velocity: ', round(vx), 'm/s'
print 'Vertical Veloctity: ', round(vy), 'm/s'
print 'Orbit Height: ', round(y/1000), 'km'

plt.show()
# Required Total Energy as function of orbit altitude
mu = 398600.4418
Re = 6371.0
Vmax = targetvel/1000
hmax = targetalt/1000
rid = np.linspace(Re,Re+hmax)
Vid = np.linspace(0,Vmax)
[X, Y] = np.meshgrid(rid,Vid)

Etot = Y**2/2-mu/X
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

fig = plt.figure()
ax = fig.gca(projection='3d')

surf = ax.plot_surface(X, Y,Etot, rstride=1, cstride=1, cmap=cm.coolwarm,
        linewidth=0, antialiased=False)
r = Re+np.array(ytab)/1000
V = np.array(vtab)/1000
Etot = V**2/2-mu/r
ax.plot(r,V,zs=Etot, zdir='z', label='zs=0, zdir=z')
Etot = Vid**2/2-mu/rid
ax.plot(rid,Vid,zs=Etot)


plt.show()
print "[done]"

