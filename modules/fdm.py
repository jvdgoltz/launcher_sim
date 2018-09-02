'''flight dynamics model'''
import modules.thrust as th
import modules.isatmos as isatmos
import modules.aerodynamics as aero
import numpy as np

g0 = 9.81

def fdm(vehicle,state,throttle,T_angle,dt):
    g=g0*(6371000/(6371000+state[1]))**2
    W=state[9]*g

    ve = th.vexit(vehicle.epsilon, vehicle.Tc, vehicle.Mw, vehicle.kappa)
    pa = isatmos.press(state[1])

    Z = state[9]*state[3]*state[3]/(state[1]+6371000)
    mflw = throttle*th.mdot(vehicle.At,vehicle.pc,vehicle.Tc,vehicle.Mw,vehicle.kappa)
    T = th.FT(mflw,ve,vehicle.At,vehicle.epsilon,vehicle.pc,vehicle.pratio,pa)

    Tx = T*np.cos(state[2]+T_angle)
    Ty = T*np.sin(state[2]+T_angle)

    Dx, Dy, M, q = aero.drag(state[1],state[3],state[4],vehicle.cd,vehicle.S) 
    D = np.sqrt(Dx**2+Dy**2)
    
    Rx = Dx+Tx        #Resulting Force on the Rocket
    Ry = Dy+Ty-W+Z
    
    # Update states
    m = state[9] - mflw*dt
    fm = state[16] - mflw*dt

    I = state[10]
    cg = state[11]

    ax = Rx/m                 #acceleration due to Resultant
    ay = Ry/m
    th_ddot = cg*T*np.sin(T_angle)/I
    
    vx = state[3]+ax*dt
    vy = state[4]+ay*dt
    th_dot = state[5] + th_ddot*dt
    
    x = state[0]+vx*dt
    y = state[1]+vy*dt
    theta = state[2] + th_dot*dt
    
    state = np.array([
            x,      #0
            y,      #1
            theta,  #2
            vx,     #3
            vy,     #4
            th_dot, #5
            ax,     #6
            ay,     #7
            th_ddot,#8
            m,      #9
            I,      #10
            cg,     #11
            T,      #12
            M,      #13
            D,      #14
            q,      #15
            fm])    #16
    
    return state

def linearize(vehicle,state,throttle,T_angle,dt,eps=0.1):
    A = np.zeros([len(state),len(state)])
    B = np.zeros([len(state),2])
    for i in range(len(state)):
        state_i1 = state
        state_i2 = state
        state_i2[i] = state[i]+eps
        state_i1[i] = state[i]-eps
        A[:,i] = (fdm(vehicle,state_i2,throttle,T_angle,dt)-fdm(vehicle,state_i1,throttle,T_angle,dt))/(2*eps)
    
    throttle2 = throttle+eps
    throttle1 = throttle-eps
    B[:,0] = (fdm(vehicle,state,throttle2,T_angle,dt)-fdm(vehicle,state,throttle1,T_angle,dt))/(2*eps)
    
    T_angle2 = T_angle+eps
    T_angle1 = T_angle-eps
    B[:,1] = (fdm(vehicle,state,throttle,T_angle2,dt)-fdm(vehicle,state,throttle,T_angle1,dt))/(2*eps)
        
    return A,B
