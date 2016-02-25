import numpy as np

def acceleration(r, v, m, v_h = 230.160, q = 0.9, r_c = 0):
    x, y, z = r
    accel = -v_h*v_h*r/(x*x+y*y+(z/q)*(z/q)+r_c*r_c)
    m_dot = np.zeros(1)
    veloc = np.zeros(3)
    return np.concatenate((veloc,accel,m_dot))


def density(r, v, m, v_h = 230.160, q = 0.9, r_c = 0):
    x, y, z = r
    rho = 2*q*q*v_h*v_h*(a*a*q*q*(2*q*q+1)-2*q*q*q*q*(x*x+y*y)+q*q*(x*x+y*y+2*z*z)-3*z*z)
    rho /= (a*a*q*q+q*q*(x*x+y*y)+z*z)**3
    return rho
