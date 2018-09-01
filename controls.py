import math
import numpy as np

def reference(x,y,xmax,targetalt):
    ymax = targetalt
    if x<xmax:
        y_ref = math.sqrt(abs(ymax**2 * (1 - (abs(x)- xmax)**2/xmax**2)))
        xs    = math.sqrt(abs(xmax**2 * (1 - (abs(y))**2/ymax**2)))
        xs = x
        yp = ymax**2*(xmax-abs(xs)) / (xmax**2 *np.sqrt( (ymax**2 * abs(xs) *(2*xmax-abs(xs)) / (xmax**2) )))
        ga_ref = np.arctan2(yp,1)
    else:
        y_ref = ymax
        ga_ref = 0
    return ga_ref, y_ref
    
def pitchrate(x,y,vx,vy,gamma,omega,targetvel,xmax,ymax,t_ctl,W,Z,T,Dx,Mf,mflw,M,t,t_fin):

    K1 = 2
    K2 = 1./math.pi
    K3 = 1

    if y>0000:

        ref = reference(x,y,xmax,ymax)
        ga_ref = ref[0]
        y_ref = ref[1]

        ay_req = -2*vy/(Mf/mflw) + (ymax-y)/(Mf/mflw)**2
        ax_req = (targetvel-vx)/(Mf/mflw)*0
        
        Ry_req = ay_req *M *0
        Rx_req = ax_req *M

        a = (Ry_req+W-Z)/T
        b = (Rx_req-Dx) /T
        
        K  = K1*(math.pi/2-ga_ref)/(ymax-y)
        deltaga = K  * (y_ref-y)
        ga_min = np.arcsin(a)
        ga_req = ga_ref + deltaga
        if ga_req<ga_min:
            ga_req=ga_min

        omega_ref   = K2 * (ga_req-gamma) 

    else:
        omega_ref = 0
    if omega>0:
        alpha = 0    
    else:
        if T<>0:
            alpha = K3 * (omega_ref - omega)
        else:
            alpha = 0
    
    return alpha