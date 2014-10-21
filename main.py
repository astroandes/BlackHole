import numpy as np,sympy as sp, pylab
from scipy.special import erf

# Declares some global constants
N = 3
n_iter = 7000
m = 1.0
rho_0 = 1.0
G = 1.0
x = np.empty((n_iter,N))
v = np.empty((n_iter,N))
x[0] = np.array([0,1,0])
v[0] = (2**(0.5))*np.array([-1,0,0])

# Defines a provisional potential
def phi(x,y,z):
    return 2.0*np.pi*G*rho_0*sp.log(x*x+y*y+z*z)

# Euler integration method (provisional)
def integrate(F,x_0,v_0,dt,n_steps):
    x = x_0
    v = v_0
    for i  in range(n_steps):
        a = F(x,v)/m
        x += v*dt
        v += a*dt
    return x,v

# Calculates gradient of an arbitrary potential using sympy
def getForce(Potential):
    x,y,z = sp.symbols('x y z')
    print -sp.diff(Potential(x,y,z),x)
    Fx = sp.lambdify((x,y,z),-sp.diff(Potential(x,y,z),x),'numpy')
    Fy = sp.lambdify((x,y,z),-sp.diff(Potential(x,y,z),y),'numpy')
    Fz = sp.lambdify((x,y,z),-sp.diff(Potential(x,y,z),z),'numpy')
    return Fx,Fy,Fz

Fx,Fy,Fz = getForce(phi)

# Defines total forces: gravitational + dinamic friction
def f(r,v):
    x,y,z = r
    b = np.sqrt(np.sum(v*v))
    return np.array([Fx(x,y,z),Fy(x,y,z),Fz(x,y,z)])-1.0*(erf(b)-2.0*b*np.exp(-0.5*b*b))*v/b

for i in range(n_iter-1):
    x[i+1],v[i+1] = integrate(f,x[i],v[i],0.001,1)

print x
pylab.plot(x[:,0],x[:,1])
pylab.show()
