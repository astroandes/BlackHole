from scipy.special import __erf

def acceleration(s, t, G = 1, M = lambda r : 1, logLambda = 1, rho = lambda r : 1, sigma = 1):
    ans = np.empty(6)
    ans[:3] = 0
    X =  np.sum(s[3:]**2)**(1/2)/(2*sigma)
    ans[3:] = s[3:]*(__erf(X) - 2*X*np.exp(-X*X)/np.pi**(1/2))
    ans[3:] *=  -(4*np.pi*logLambda*G*rho(s[:3])*M(s[:3])/np.sum(s[3:]**2)**(3/2))
    return ans
