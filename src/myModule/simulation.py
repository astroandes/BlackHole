import numpy as np

class interaction:
    def __init__(self, name, acceleration, parameters):
        self.name = name
        self.acceleration = lambda r, v, m: acceleration(r, v, m, parameters)
        self.parameters = parameters

class simulation:
    def __init__(self):
        from scipy.integrate import odeint
        self.integrator = odeint
        self.interactions = []
        self.derivatives = lambda s,t: np.array([s[3],s[4],s[5],0,0,0,0])
        self.path = None

    def integrate(self, t_init, t_final, n_steps, s_init):
        t = np.linspace(t_init, t_final, n_steps)
        self.path = self.integrator(self.derivatives, s_init, t)

    def updateDerivatives(self):
        self.derivatives = lambda s,t: np.sum(np.array([f.acceleration(s[0:3],s[3:6],s[6]) for f in self.interactions]),axis = 0) + np.array([s[3],s[4],s[5],0,0,0,0])
