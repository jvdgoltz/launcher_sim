#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Sep  1 22:24:28 2018

@author: jvdgoltz
"""
import numpy as np
from cvxopt import matrix, solvers, sparse, spdiag, spmatrix

def solve(x0,targetalt,targetvel,amax,dt=0.5,t_max=12*60):
    r = np.array([0,targetalt,targetvel,0,0,0,0,9.81])
        
    t = np.arange(0,t_max+dt,dt)
        
    # Inequality Constraints
    G = spmatrix(np.zeros([1,6*len(t)]), range(6*len(t)), range(6*len(t)), (6*len(t),8*len(t)))
    h = np.zeros([6*len(t),1])
    for i in range(len(t)):
        G[2*i:2*i+2,8*i:8*i+2] = -np.eye(2)
        G[2*len(t)+4*i:2*len(t)+4*i+4,8*i+6:8*i+8] = np.array([[1,0],[-1,0],[0,1],[0,-1]])
        if i<len(amax):
            h[2*len(t)+4*i:2*len(t)+4*i+4,0] = np.ones(4)*amax[i]
        else:
            h[2*len(t)+4*i:2*len(t)+4*i+4,0] = np.zeros(4)
    h = matrix(h)
    
    # Equality Constraints
    A = np.array([[1,0,dt,0,0,0],
                  [0,1,0,dt,0,0],
                  [0,0,1,0,dt,0],
                  [0,0,0,1,0,dt],
                  [0,0,0,0,0,0],
                  [0,0,0,0,0,0]])
    B = np.array([[0,0],
                 [0,0],
                 [0,0],
                 [0,0],
                 [1,0],
                 [0,1]])
    Aeq = spmatrix(np.zeros([1,6*len(t)]),range(6*len(t)),range(6*len(t)),(6*len(t),8*len(t)))
    Aeq[:6,:6] = sparse(matrix(np.eye(6)))
    c = sparse(matrix(np.concatenate((-A,-B,np.eye(6)),axis=1)))
    for i in range(len(t)-1):    
        Aeq[6*i+6:6*i+12,8*i:8*i+14] = c
    beq = np.zeros([6*len(t),1])
    beq[:6,0] = x0
    for i in range(len(t)):
        beq[6*i+5,0] = -9.81
    beq = matrix(beq)

    # Objective
    Q = matrix(2*np.eye(8))
    Q[0,0] = 0
    QQ = sparse(Q)
    for i in range(len(t)-1):
        QQ = spdiag([QQ,Q])
    p = -r.T.dot(Q)
    pp = matrix(np.kron(np.ones([1,len(t)]), p).T)
    
    
    sol = solvers.qp(QQ, pp, G, h, Aeq, beq)
    
    x = np.array(sol['x']).reshape((-1,8))
    return t,x