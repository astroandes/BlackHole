import numpy as np
def speed(s):
    ans = np.empty(6)
    ans[:3] = s[3:]
    ans[3:] = 0
    return ans

class interaction:
    def __init__(self, name, acceleration, **constants):
        self.name = name
        self.acceleration = lambda s,t: acceleration(s, t, **constants)

class simulation:
    def __init__(self, integratorName):
        if integratorName == 'odeint':
            from scipy.integrate import odeint
            self.integrator = odeint
        else:
            raise ValueError('The integrator was not found')
        self.integratorName = integratorName
        self.interactions = []
        self.acceleration = speed
        self.path = None

    def updateAcceleration(self):
        self.acceleration = lambda s, t: sum([f.acceleration(s,t) for f in self.interactions]) + speed(s)

    def integrate(self, initial_state, initial_time, final_time, resolution):
        t = np.arange(initial_time, final_time, resolution)
        self.path = self.integrator(self.acceleration, initial_state, t)
