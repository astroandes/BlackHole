import numpy as np,sympy as sp, pylab
from scipy.special import erf
from scipy.integrate import odeint as integrate
## Units
# Mass:     10**6 Solar Masses
# Distance: Parsecs
# Time:     10**6 Years

## General Constants

N = 3                   # Number of dimensions
G = 4490.0              # Gravitational constant
n_iter = 50000          # Number of iterations

## System Constants

lnLamb = 1.0            # Coulomb logarithm
sigma = 76.65           # Velocity dispersion
M = 3.0                 # BlackHole mass
M_tot = 2.0*M           # Total mass
R = G*M/(sigma*sigma)   # BlackHole radius
K = 4.0*np.pi*G*G*lnLamb*M_tot

## Density and Mass function for the gas as numpy functions

def M(r):
    return 2.0*sigma*sigma*(r-R*np.arctan(r/R))/G

def rho(r):
    return sigma*sigma/(2.0*np.pi*G*(r*r+R*R))

# Defines a provisional potential
def phi(x,y,z):
    #   return 2.0*np.pi*G*rho_0*sp.log(x*x+y*y+z*z)
    return -G/sp.sqrt(x*x+y*y+z*z)

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
    r = np.sqrt(x*x+y*y+z*z)
    v = s[N:]
    V = np.sqrt(np.sum(v*v))
    b = V/(np.sqrt(2.0)*sigma)
    fx,fy,fz = -K*rho(r)*v*(erf(b)-2.0*b*np.exp(-0.5*b)/np.sqrt(np.pi))/(V*V*V)
    return np.array([v[0],v[1],v[2],(M(r)*Fx(x,y,z)+fx),(M(r)*Fy(x,y,z)+fy),(M(r)*Fz(x,y,z)+fz)])


# Initial conditions
vlct = np.array([200,300,400])*1.022
t = np.linspace(0,1000.0,n_iter)
for vel in vlct:
    print vel
    s_0 = np.array([1,0,0,0,0,vel]) #[x,y,z,v_x,v_y,v_z]
    s = integrate(f,s_0,t)
    #E = 0.5*m*np.sum(v*v,axis=1)+2.0*np.pi*G*rho_0*np.log(x[:,0]*x[:,0]+x[:,1]*x[:,1]+x[:,2]*x[:,2])
    x,v = s[:,:N],s[:,N:]
    pylab.plot(t*1E6,np.sqrt(x[:,0]*x[:,0]+x[:,1]*x[:,1]+x[:,2]*x[:,2]),label='$\mathrm{'+str(vel)+'\ Pc/Myr}$')

pylab.legend(loc=2)
pylab.xlabel('$\mathrm{Time\ (Years)}$',fontsize=16)
pylab.ylabel('$\mathrm{Radial\ Distance\ (Pc)}$',fontsize=16)
pylab.xscale('log')
pylab.yscale('log')
pylab.ylim(ymin=3)
pylab.xlim(xmin=2E4)
pylab.savefig('r.png',dpi=200)
pylab.close()
