import numpy as np

def acceleration(r, v, m, G = 4490.0, a = 6.5e3, b = 0.26e3, M = 1e5):
    x, y, z = r
    accel = -G*M*r/((x*x+y*y+(a+np.sqrt(b*b+z*z))**2)**1.5)
    accel[2] *= (1+a/np.sqrt(b*b+z*z))
    m_dot = np.zeros(1)
    veloc = np.zeros(3)
    return np.concatenate((veloc,accel,m_dot))

def density(r, a = 6.5e3, b = 0.26e3, M = 1e5):
    x, y, z = r
    R2 = np.sum(r*r)
    rho = M*b*b*(a*R2+(a+3*(z*z+b*b)**.5)*(a+(z*z+b*b)**.5)**2)
    rho /=4*np.pi*((R2+(a+(z*z+b*b)**.5)**2)**2.5)*((z*z+b*b)**1.5)
    return rho
