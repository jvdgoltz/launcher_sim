import numpy as np
import modules.fdm as fdm

def control(state,reference):
    #TODO add control law here:
    throttle = 5
    T_angle = 0
    if state[16] <= 0:
        sep = True
    else:
        sep = False 
    return throttle, T_angle, sep

