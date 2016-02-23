import numpy as np
from scipy.special import erf

def acceleration(r, v, m, G = 4490.0, logLambda = 1.0, rho = lambda r_norm : 1, sigma = 76.65):
    v_norm =  np.sqrt(np.sum(v*v))
    r_norm =  np.sqrt(np.sum(r*r))
    X =  v_norm/(np.sqrt(2)*sigma)
    accel = -4*np.pi*G*G*logLambda*m*v*rho(r_norm)*(erf(X)-2.*X*np.exp(-.5*X)/np.sqrt(np.pi))/v_norm**3
    m_dot = np.zeros(1)
    veloc = np.zeros(3)
    return np.concatenate((veloc,accel,m_dot))
