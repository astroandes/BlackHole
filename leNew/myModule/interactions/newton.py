import numpy as np

def acceleration(s, t, G = 1, m = lambda x: 1, center = np.zeros(3)):
    """
    Inverse square law for gravity. If a body is in a position $\\vec{r}$ with
    velocity $\\vec{v}$ under the influence of another body with mass $m(r)$
    and position $\\vec{r_{0}}$
    this function returns
    $$
    Gm\\frac{\\vec{r}_{0}-\\vec{r}}{\left\Vert \\vec{r}_{0}-\\vec{r}\\right\Vert ^{3}}
    $$
    Where $G$ is the Gravitational constant.

    Parameters
    ----------
    s : numpy array
        6D vector of position and velocity
    t : float
        time value at an specific moment
    G : float
        Gravitational constant
    m : function
        Mass of the body producing the force as a function of the radial
        distance to itself.
    center : numpy array
        3D vector with the position of the body producing the force

    Returns
    -------
    numpy array
        A 6 component numpy array where the first three are zero and the other
        three are the values of the acceleration in x, y and z.
    """
    ans = np.empty(6)
    ans[:3] = 0
    ans[3:] = center-s[:3]
    ans[3:] /= np.sum(ans[3:]**2)**(3/2)
    ans[3:] *= G*m(np.sqrt(np.sum(s[:3]**2)))
    return ans

def potential(s, t, G = 1, m = lambda x: 1, center = np.zeros(3)):
    """
    Inverse square law for gravity. If a body is in a position $\\vec{r}$ with
    velocity $\\vec{v}$ under the influence of another body with mass $m(r)$
    and position $\\vec{r_{0}}$
    this function returns
    $$
    -\\frac{Gm}{\left\Vert \\vec{r}_{0}-\\vec{r}\\right\Vert}
    $$
    Where $G$ is the Gravitational constant.

    Parameters
    ----------
    s : numpy array
        6D vector of position and velocity
    t : float
        time value at an specific moment
    G : float
        Gravitational constant
    m : function
        Mass of the body producing the force as a function of the radial
        distance to itself.
    center : numpy array
        3D vector with the position of the body producing the force

    Returns
    -------
    float
        The value of the gravitational potential (Energy per mass unit)
    """
    return G*m(np.sqrt(np.sum(s[:3]**2)))/np.sum((center-s[:3])**2)**(1/2)
