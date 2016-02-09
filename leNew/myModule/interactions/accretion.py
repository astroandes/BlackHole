import numpy as np

def acceleration(s, t, m = lambda s,t: 1, m_dot = lambda s,t: 0):
    """
    acceleration due to mass accretion suffered by a body of mass $m$, mass accretion $\dot{m}$
    and velocity $\\vec{v}$ this function returns
    $$
    -\\frac{\dot{m}}{m} \\vec{v}
    $$
    Parameters
    ----------
    s : numpy array
        6D vector of position and velocity
    t : float
        time value at an specific moment
    m : function
        Mass of the body suffering mass accretion.
    m_dot : function
        Mass accretion of the body suffering mass accretion.

    Returns
    -------
    numpy array
        A 6 component numpy array where the first three are zero and the other
        three are the values of the acceleration in x, y and z.
    """
    ans = np.empty(6)
    ans[:3] = 0
    ans[3:] = - m_dot(s,t) * s[3:] / m
