import numpy as np
from scipy.special import erf

def acceleration(s, t, G = 1, m = 1, logLambda = 1, rho = lambda r : 1, sigma = 1):
    ans = np.empty(6)
    ans[:3] = 0
    v =  np.sqrt(np.sum(s[3:]*s[3:]))
    r =  np.sqrt(np.sum(s[:3]*s[:3]))
    X =  v/(np.sqrt(2)*sigma)
    ans[3:] =  s[3:]*rho(r)*(erf(X)-2.*X*np.exp(-.5*X)/np.sqrt(np.pi))/(v**3)
    ans[3:] *= - 4*np.pi*G*G*logLambda*m
    return ans
