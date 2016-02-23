import numpy as np

def acceleration(r, v, m, v_h = 230160, q = 0.9, r_c = 0):
    x, y, z = r
    accel = -v_h*v_h*r/(x*x+y*y+(z/q)*(z/q)+r_c*r_c)
    m_dot = np.zeros(1)
    veloc = np.zeros(3)
    return np.concatenate((veloc,accel,m_dot))
