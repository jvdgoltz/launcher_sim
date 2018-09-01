import modules.isatmos as isa
import math
def drag(h,vx,vy,cd0,S):
    rho = float(isa.dens(h))
    Temp = float(isa.temp(h))
    v = math.sqrt(vx*vx+vy*vy)
    a = math.sqrt(1.4*isa.R*Temp)
    if a>0:
        M = v/a
    else:
        M = 0

    q = 0.5 * rho * v**2
    if M<1:
        cd = cd0 / ( math.sqrt(1-M**2) + (M**2 / (1+math.sqrt(1-M**2)))*cd0/2)
    else:
        cd = cd0 / math.sqrt(M**2-1)
    Dx=-cd*0.5*rho*v*vx*S
    Dy=-cd*0.5*rho*v*vy*S
    return Dx, Dy, M, q

