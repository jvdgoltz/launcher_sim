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
import xlrd
import numpy as np
import matplotlib.pyplot as plt
import modules.aerodynamics as aero
import modules.thrust as th
import modules.controls as ctl
import modules.fdm as fdm
import modules.guidance as guid
from modules.isatmos import press

class vehicle:
    def __init__(self,fm,epsilon,Tc,Mw,kappa,At,pc,pratio,
        cd,S,nengines):
        self.fm = fm
        self.epsilon = epsilon
        self.Tc = Tc
        self.Mw = Mw
        self.kappa = kappa
        self.At = At
        self.pc = pc
        self.pratio = pratio
        self.cd = cd
        self.S = S
        self.nengines = nengines

#read out data source

#filename=raw_input("Please enter the data source file name.")
filename="launcher_data/Saturn5.xls"
table = xlrd.open_workbook(filename)
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

ttab=[0]
mtab=[sum(oem)+sum(fm)+payld]
Ttab=[0]
Dtab=[0]
vtab=[0]
atab=[0]
xtab=[0]
ytab=[0]
thtab = [0]
qtab=[0]
Mtab=[0]
m = sum(oem)+sum(fm)+payld
I = 1.
cg = 0.
t = 0
x = 0.
y = 0.
theta = 90.*np.pi/180
vx = 0.
vy = 0.
th_dot = 0.
ax = 0.
ay = 0.
th_ddot = 0.
g0 = 9.81
T = 0.
M = 0.
D = 0.
q = 0.
state = np.array([x,y,theta,vx,vy,th_dot,ax,ay,th_ddot,m,I,cg,T,M,D,q,fm[0]])
startstage = []
reference = targetalt
for i in range(nstages):
    print("Simulating stage:",i+1)
    if i>0:
        fm[i-1]=0
        oem[i-1]=0
    state[9] = sum(oem)+sum(fm)+payld
    state[16] = fm[i]
    stage = vehicle(fm[i],epsilon[i],Tc[i],Mw[i],kappa[i],At[i],pc[i],pratio[i],cd[i],S[i],nengines[i])
    startstage.append(t)
    sep = False
    while y>=0 and t<2000 and not sep:
        t = t + dt
        ttab.append(t)

        throttle, T_angle, sep = ctl.control(state,reference)
        
        #A,B = fdm.linearize(stage,state,throttle,T_angle,dt)
        state = fdm.fdm(stage,state,throttle,T_angle,dt)
        y = state[1]
        fm[i] = state[16]
        
        mtab.append(state[9])
        xtab.append(state[0])
        ytab.append(y)
        Ttab.append(state[12])
        vtab.append(np.sqrt(state[3]**2+state[4]**2))
        qtab.append(state[15])
        thtab.append(state[2])
        atab.append(np.sqrt(state[6]**2+state[7]**2))
        Mtab.append(state[13])
        Dtab.append(state[14])

ttab = np.array(ttab)/60
plt.subplot(331)
plt.plot(ttab, np.array(ytab)/1000)
plt.title("flight profile: altitude [km] vs time [min]")

plt.subplot(332)
plt.plot(np.array(xtab)/1000, np.array(ytab)/1000, 'b')
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
plt.plot(ttab, np.array(thtab)*180/np.pi,'b')
plt.title("pitch angle [deg] vs time [min]")

plt.subplot(337)
plt.plot(ttab, np.array(atab)/g0+1)
plt.title("acceleration [g0] vs time [min]")

plt.subplot(333)
plt.plot(ttab, np.array(Mtab)) 
plt.title("Mach number [-] vs time [min]")

plt.subplot(339)
plt.plot(ttab, np.array(Dtab)/1000) 
plt.title("Drag [kN] vs time [min]")

print('Orbit Insertion Velocity: ', round(targetvel), 'm/s')
print('Horizontal Velocity: ', round(state[3]), 'm/s')
print('Vertical Veloctity: ', round(state[4]), 'm/s')
print('Orbit Height: ', round(state[1]/1000), 'km')

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
print("done")