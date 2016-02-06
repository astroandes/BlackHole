import numpy as __np

def acceleration(s, t, G = 0, m = 0, center = __np.zeros(3)):
    ans = __np.empty(6)
    ans[:3] = 0
    ans[3:] = center-s[:3]
    ans[3:] /= __np.sum(ans[3:]**2)**(3/2)
    return ans

def potential(s, t, G = 0, m = 0, center = __np.zeros(3)):
    return G*m/__np.sum((center-s[:3])**2)**(1/2)
