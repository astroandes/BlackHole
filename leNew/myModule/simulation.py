import numpy as np
def speed(s):
    ans = np.empty(6)
    ans[:3] = s[3:]
    ans[3:] = 0
    return ans

class interaction:
    """
    This class is an abstraction of an interaction inside a simulation.

    Attributes
    ----------
    name : str
        Human readable string describing the interaction.
    acceleration : function
        A numpy-compatible function which first argument is the 6D-vector of
        position + velocity, the other arguments must be keyword arguments.
        Unless it is desired, this function should not return any values on
        the first three components (which repesent velocity)
    """
    def __init__(self, name, acceleration, **constants):
        """
        Parameters
        ----------
        name : str
            Human readable string describing the interaction.
        acceleration : function
            A numpy-compatible function which first argument is the 6D-vector of
            position + velocity, the other arguments must be keyword arguments.
            Unless it is desired, this function should not return any values on
            the first three components (which repesent velocity)
        **constants
            Arbitrary keyword arguments of the acceleration function.
        """
        self.name = name
        self.acceleration = lambda s,t: acceleration(s, t, **constants)

class simulation:
    """
    This class is an abstraction of a numerical one-body simulation.

    Attributes
    ----------
    integrator : function
        Integrator function which takes as arguments a function and two numpy
        arrays, the first one is the initial state of the simulation, and the
        second one is the time array. It must return a numpy array object.
    integratorName : str
        Name of the integrator method. Currently it only recieves the \'odeint\'
        string which uses the scipy.integrate.odeint function.
    interactions : list
        A list of interaction objects.
    acceleration : function
        The total acceleration function which arguments are the 6D vector of
        position and velocity and time as a float.
    path : numpy array
        A numpy array containing all the values of 6D vector of position and velocity
        during time after using the integrate method.
    """
    def __init__(self, integratorName):
        """
        Parameters
        ----------
        integratorName: str
            Name of the integrator method. Currently it only recieves the \'odeint\'
            string which uses the scipy.integrate.odeint function.
        """
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
        """
        Method capable of updating the acceleration attribute. It must be called
        after appending an interaction to the interactions attribute.
        """
        self.acceleration = lambda s, t: sum([f.acceleration(s,t) for f in self.interactions]) + speed(s)

    def integrate(self, initial_state, initial_time, final_time, resolution):
        """
        Method capable of integrating the path of a particle under the influence
        of the acceleration method attribute using the integrator attribute.bIt saves the
        result into the path attribute.

        Parameters
        ----------
        initial_state : numpy array
            6D vector with the initial values of position and velocity in that
            order.
        initial_time : float
            Time in which the simulation starts. Tipically it is zero and its
            value is irrelevant if the acceleration attribute does not depends
            on time.
        final_time : float
            Time in which the simulation ends
        resolution : float
            The discrete time interval used to integrate numerically.
        """
        t = np.arange(initial_time, final_time, resolution)
        self.path = self.integrator(self.acceleration, initial_state, t)
