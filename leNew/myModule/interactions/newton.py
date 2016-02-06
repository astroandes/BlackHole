import numpy as np

def acceleration(s, t, G = 0, m = 0, center = np.zeros(3)):
    ans = np.empty(6)
    ans[:3] = 0
    ans[3:] = center-s[:3]
    ans[3:] /= np.sum(ans[3:]**2)**(3/2)
    return ans

def potential(s, t, G = 0, m = 0, center = np.zeros(3)):
    return G*m/np.sum((center-s[:3])**2)**(1/2)
