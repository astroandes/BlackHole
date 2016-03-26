import numpy as np

def acceleration(r, v, m, parameters):
    G, M, R = parameters
    dr = R - r
    dr_norm = np.sum(dr*dr)**.5
    accel = G*M(dr_norm)*dr/dr_norm**3
    m_dot = np.zeros(1)
    veloc = np.zeros(3)
    return np.concatenate((veloc,accel,m_dot))
