import numpy as np

def acceleration(r, v, m, G = 1, a = 1, b = 1):
    x, y, z = r
    accel = -G*r/(x*x+y*y+(a+np.sqrt(b*b+z*z))**2)**1.5
    accel[2] *= (1+a/np.sqrt(b*b+z*z))
    m_dot = np.zeros(1)
    veloc = np.zeros(3)
    return np.concatenate((veloc,accel,m_dot))

def density(r, a = 1, b = 1):
    x, y, z = r
    rho = b*b*(a*a*a+5*a*a*np.sqrt(b*b+z*z)+a*(7*b*b+x*x+y*y+7*z*z)+3*(b*b+z*z)**(1.5))
    rho /= 4*np.pi*(b*b+z*z)**(1.5)*((a+np.sqrt(b*b+z*z))**2*x*x+y*y)**(2.5)
    return rho
