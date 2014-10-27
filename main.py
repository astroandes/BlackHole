import numpy as np,sympy as sp, pylab
from scipy.special import erf
from scipy.integrate import odeint as integrate

# Declares some global constants
N = 3 # Number of dimensions
G = 1.0 # Gravitational constant
n_iter = 50000 # Number of iterations
m = 1.0 # BH mass
R = 1.0 # BH radius

# Constants for the gravitational potential
rho_0 = 1.0

# Constants for the DM friction force
lnLambda= 1.0
M_tot = 0.01
rho = 1.0
sigma = 1.0
K = 4.0*np.pi*G*G*lnLambda*rho*M_tot

# Initial conditions
s_0 = np.array([1,0,0,-2,-2,0]) #[x,y,z,v_x,v_y,v_z]

t = np.linspace(0,1000.0,n_iter)

# Defines a mass distribution for the gas
def M(x,y,z):
    r= np.sqrt(x*x+y*y+z*z)
    return 2.0*sigma*sigma*r*(1.0-R/r*np.arctan(r/R))/G

# Defines a provisional potential
def phi(x,y,z):
    #   return 2.0*np.pi*G*rho_0*sp.log(x*x+y*y+z*z)
    return -G*m/sp.sqrt(x*x+y*y+z*z)

# Calculates gradient of an arbitrary potential using sympy
def getForce(Potential):
    x,y,z = sp.symbols('x y z')
    Sx = -sp.diff(Potential(x,y,z),x)
    Sy = -sp.diff(Potential(x,y,z),y)
    Sz = -sp.diff(Potential(x,y,z),z)
    print 'The Force is:\n',Sx,'\n',Sy,'\n',Sz,'\n'
    Fx = sp.lambdify((x,y,z),Sx,'numpy')
    Fy = sp.lambdify((x,y,z),Sy,'numpy')
    Fz = sp.lambdify((x,y,z),Sz,'numpy')
    return Fx,Fy,Fz

Fx,Fy,Fz = getForce(phi)

# Defines the time derivative of s: [v_x,v_y,v_z,a_x,a_y,a_z]
def f(s,t):
    x,y,z = s[:N]
    v = s[N:]
    V = np.sqrt(np.sum(v*v))
    b = V/(np.sqrt(2.0)*sigma)
    fx,fy,fz = -K*v*(erf(b)-2.0*b*np.exp(-0.5*b)/np.sqrt(np.pi))/(V*V*V)
    return np.array([v[0],v[1],v[2],(M(x,y,z)*Fx(x,y,z)+fx)/m,(M(x,y,z)*Fy(x,y,z)+fy)/m,(M(x,y,z)*Fz(x,y,z)+fz)/m])

s = integrate(f,s_0,t)

x,v = s[:,:N],s[:,N:]

#E = 0.5*m*np.sum(v*v,axis=1)+2.0*np.pi*G*rho_0*np.log(x[:,0]*x[:,0]+x[:,1]*x[:,1]+x[:,2]*x[:,2])
#print E
'''
pylab.plot(x[:,0],x[:,1])
pylab.xlabel('$x$')
pylab.ylabel('$y$')
pylab.axes().set_aspect('equal','datalim')
pylab.savefig('xy.png',dpi=200)
pylab.close()
'''
pylab.plot(t,np.sqrt(x[:,0]*x[:,0]+x[:,1]*x[:,1]))
pylab.xlabel('$t$',fontsize=16)
pylab.ylabel('$r=\sqrt{x^{2}+y^{2}+z^{2}}$',fontsize=16)
pylab.xscale('log')
pylab.yscale('log')
pylab.ylim([1,100])
pylab.savefig('r.png',dpi=200)
pylab.close()
