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
from modules.isatmos import press

class stage:
    def __init__(self,oem,fm,pld,epsilon,Tc,Mw,kappa,At,pc,pratio,
        cd,S,nengines):
        self.oem = oem
        self.fm = fm
        self.pld = pld
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

class rocket:
    def __init__(self,stagelist):
        self.nstages = len(stagelist)
        self.stagelist = stagelist
        self.oem = 0
        self.fm = 0
        self.pld = self.stagelist[-1].pld
        for i in range(nstages):
            self.oem = self.oem + self.stagelist[i].oem
            self.fm = self.fm + self.stagelist[i].fm
        self.m0 = self.oem+self.fm+self.pld
        
    def getMass(self):
        self.oem = 0
        self.fm = 0
        self.pld = self.stagelist[-1].pld
        for i in range(nstages):
            self.oem = self.oem + self.stagelist[i].oem
            self.fm = self.fm + self.stagelist[i].fm
        self.m = self.oem+self.fm+self.pld
        return self.m
        
            
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
stagelist = []
for i in range(nstages):
    oem = (sheet.cell_value(8,1+i))
    fm = (sheet.cell_value(9,1+i))
    Tc =(sheet.cell_value(10,1+i))
    pc =(sheet.cell_value(11,1+i))
    epsilon =(sheet.cell_value(12,1+i))
    Mw = (sheet.cell_value(13,1+i))
    kappa = (sheet.cell_value(14,1+i))
    Ae = (sheet.cell_value(15,1+i))
    At = np.divide(Ae,epsilon)
    nengines = (sheet.cell_value(16,1+i))
    cd = (sheet.cell_value(17,1+i))
    S = (sheet.cell_value(18,1+i))
    t_ctl = (sheet.cell_value(19,1+i))
    sepdur = (sheet.cell_value(20,1+i))
    pratio = th.pratio(epsilon,kappa)
    stagelist.append(stage(oem,fm,0,epsilon,Tc,Mw,kappa,At,pc,pratio,cd,S,nengines))
stagelist[-1].pld = payld    
vehicle = rocket(stagelist)

m = vehicle.m0
I = 1.
cg = 0.
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
state = np.array([x,y,theta,vx,vy,th_dot,ax,ay,th_ddot,m,I,cg,T,M,D,q,vehicle.fm])

reference = targetalt
def trajectory(vehicle,state,reference,dt=0.5,t_max=12*60):
    t = 0
    ttab=np.array([0])
    mtab=np.array([vehicle.m0])
    Ttab=np.array([0])
    Dtab=np.array([0])
    vxtab=np.array([0])
    vytab=np.array([0])
    axtab=np.array([0])
    aytab=np.array([0])
    xtab=np.array([0])
    ytab=np.array([0])
    thtab=np.array([0])
    qtab=np.array([0])
    Mtab=np.array([0])
    for i in range(vehicle.nstages):
        print("Simulating stage:",i+1)
        if i>0:
            vehicle.stagelist[i-1].fm=0
            vehicle.stagelist[i-1].oem=0
        state[9] = vehicle.getMass()
        state[16] = vehicle.stagelist[i].fm
        sep = False
        while state[1]>=0 and t<t_max and not sep:
            t = t + dt
            ttab = np.append(ttab,t)
    
            throttle, T_angle, sep = ctl.control(state,reference)
            
            #A,B = fdm.linearize(stage,state,throttle,T_angle,dt)
            state = fdm.fdm(vehicle.stagelist[i],state,throttle,T_angle,dt)
            vehicle.stagelist[i].fm = state[16]
            
            mtab = np.append(mtab,state[9])
            xtab = np.append(xtab,state[0])
            ytab = np.append(ytab,state[1])
            Ttab = np.append(Ttab,state[12])
            vxtab = np.append(vxtab,state[3])
            vytab = np.append(vytab,state[4])
            qtab = np.append(qtab,state[15])
            thtab = np.append(thtab,state[2])
            axtab = np.append(axtab,state[6])
            aytab = np.append(aytab,state[7])
            Mtab = np.append(Mtab,state[13])
            Dtab = np.append(Dtab,state[14])
    return ttab,mtab,xtab,ytab,Ttab,vxtab,vytab,qtab,thtab,axtab,aytab,Mtab,Dtab
