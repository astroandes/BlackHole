import numpy as np

def acceleration(r, v, m, G = 1, M = lambda dr_norm: 1, R = np.zeros(3)):
    dr = R - r
    dr_norm = np.sum(dr*dr)**.5
    return G*M(dr_norm)*dr/dr_norm**3
