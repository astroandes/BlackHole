import numpy as np

def NullForce(s):
    ans = np.empty(6)
    ans[:3] = s[3:]
    ans[3:] = 0
    return ans

class Force:
    def __init__(self, name, forceType, force, *potential):
        if forceType == "Conservative":
            self.force = force
            self.potential = potential
        elif forceType == "Disipative":
            if potential == None:
                raise ValueError('Force objects of type \"Disipative\" cannot have a potential')
            self.force = force
            self.potential = None
        else:
            raise ValueError('Force objects should be of type \"Conservative\" or \"Disipative\"')
        self.name = name

class Simulation:
    def __init__(self, integratorName):
        if integratorName == 'odeint':
            from scipy.integrate import odeint
            self.integrator = odeint
        else:
            raise ValueError('The integrator was not found, please set it as a custom integrator')
        self.integratorName = integratorName
        self.forces = []
        self.path = None
        self.totalForce = None

    def updateTotalForce(self):
        self.totalForce = lambda s, t: sum([f.force(s,t) for f in self.forces]) + NullForce(s)

    def addForce(self, name, forceType, force, *potential):
            if forceType == 'Disipative' and self.integratorName == 'leapfrog':
                raise ValueError('Cannot use a symplectic integrator with a disipative force!')
            newForce = Force(name, forceType, force, potential)
            self.forces.append(newForce)
            self.updateTotalForce()

    def integrate(self, initial_state, initial_time, final_time, resolution):
        t = np.arange(initial_time, final_time, resolution)
        self.path = self.integrator(self.totalForce, initial_state, t)
