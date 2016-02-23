import numpy as np

def acceleration(r, v, m, m_dot = lambda r,v,m: 0):
    accel = -v*m_dot(r,v,m)/m
    veloc = np.zeros(3)
    return np.concatenate((veloc,accel,m_dot(r,v,m)*np.ones(1)))
