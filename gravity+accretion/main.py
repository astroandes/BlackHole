import numpy as np, sympy as sp, pylab
from scipy.integrate import odeint as integrate
# Mass:     10**6 Solar Masses
# Distance: Parsecs
# Time:     10**6 Years

## universal constants
G = 4490.0
n_iter = 5000

## halo parameters
v_h = 220.0
q   = 1.0
r_c = 8.5e3

## gas parameters
m_g = 33.0e3
a = 34.86e3
b = 0.3e3
v_s = 347.5

## black hole parameters
m_0 = 3.0

m_tot = m_g + m_0

def getForce(Potential):
    x,y,z = sp.symbols('x y z')
    Sx = -sp.diff(Potential(x,y,z),x)
    Sy = -sp.diff(Potential(x,y,z),y)
    Sz = -sp.diff(Potential(x,y,z),z)
    Fx = sp.lambdify((x,y,z),Sx,'numpy')
    Fy = sp.lambdify((x,y,z),Sy,'numpy')
    Fz = sp.lambdify((x,y,z),Sz,'numpy')
    return Fx,Fy,Fz

def getDensity(Potential):
    x,y,z = sp.symbols('x y z')
    Sx = sp.diff(Potential(x,y,z),x,2)
    Sy = sp.diff(Potential(x,y,z),y,2)
    Sz = sp.diff(Potential(x,y,z),z,2)
    dens = sp.lambdify((x,y,z),(Sx+Sy+Sz)/(4*np.pi*G),'numpy')
    return dens

def halo_potential(x,y,z):
    return 0.5*v_h*v_h*sp.log(x*x+y*y+z*z/(q*q)+r_c*r_c)

def gas_potential(x,y,z):
    return -G/sp.sqrt(x*x+y*y+(a+sp.sqrt(z*z+b*b))**2)

Hx,Hy,Hz = getForce(halo_potential)
Gx,Gy,Gz = getForce(gas_potential)
rho = getDensity(gas_potential)

def mdot(x,y,z,u,v,w,m):
    if m < m_tot:
        return 4*np.pi*(m_tot-m)*rho(x,y,z)*G*G*m*m*(v_s*v_s+u*u+v*v+w*w)**-1.5
    return 0

def D(s,t):
    x,y,z = s[:3]
    u,v,w = s[3:-1]
    m = s[-1]
    dm = mdot(x,y,z,u,v,w,m)
    a_x = (Hx(x,y,z)+(m_tot-m)*Gx(x,y,z)-dm*u)/m
    a_y = (Hy(x,y,z)+(m_tot-m)*Gy(x,y,z)-dm*v)/m
    a_z = (Hz(x,y,z)+(m_tot-m)*Gz(x,y,z)-dm*w)/m
    return np.array([u,v,w,a_x,a_y,a_z,dm])

t = np.linspace(0,10000,n_iter)
s_0 = np.array([1,0,0,200/np.sqrt(3),200/np.sqrt(3),200/np.sqrt(3),m_0]) #[x,y,z,u,v,w,m]
s = integrate(D,s_0,t)
x,y,z = np.transpose(s[:,:3])
u,v,w = np.transpose(s[:,3:-1])
m = np.transpose(s[:,-1])
r = (x*x+y*y+z*z)**0.5

pylab.plot(t,r)
pylab.xlabel('$t$')
pylab.ylabel('$r(t)$')
pylab.savefig('r.png',dpi=200)
pylab.show()

r = 1.5e4
theta = np.linspace(0,2*np.pi,100)

pylab.plot(x,z)
pylab.plot(r*np.cos(theta),r*np.sin(theta),'--r')
pylab.xlabel('$x(t)$')
pylab.ylabel('$z(t)$')
pylab.savefig('x-z.png',dpi=200)
pylab.show()

pylab.plot(t,m)
pylab.xlabel('$t$')
pylab.ylabel('$m(t)$')
pylab.savefig('m.png',dpi=200)
pylab.show()
