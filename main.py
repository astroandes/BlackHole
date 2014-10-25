import numpy as np,sympy as sp, pylab
from scipy.special import erf
from scipy.integrate import odeint as integrate

# Declares some global constants
N = 3
n_iter = 5000
m = 1.0
rho_0 = 1.0
G = 1.0
s_0 = np.array([0,1,0,-2**(0.5),0,0]) #[x,y,z,v_x,v_y,v_z]

# Defines a provisional potential
def phi(x,y,z):
    return 2.0*np.pi*G*rho_0*sp.log(x*x+y*y+z*z)

# Calculates gradient of an arbitrary potential using sympy
def getForce(Potential):
    x,y,z = sp.symbols('x y z')
    print -sp.diff(Potential(x,y,z),x)
    Fx = sp.lambdify((x,y,z),-sp.diff(Potential(x,y,z),x),'numpy')
    Fy = sp.lambdify((x,y,z),-sp.diff(Potential(x,y,z),y),'numpy')
    Fz = sp.lambdify((x,y,z),-sp.diff(Potential(x,y,z),z),'numpy')
    return Fx,Fy,Fz

Fx,Fy,Fz = getForce(phi)

# Defines the time derivative of s: [v_x,v_y,v_z,a_x,a_y,a_z]
def f(s,t):
    x,y,z = s[:N]
    v = s[N:]
    b = np.sqrt(np.sum(v*v))
    fx,fy,fz = -1.0*(erf(b)-2.0*b*np.exp(-0.5*b*b))*v/b
    return np.array([v[0],v[1],v[2],(Fx(x,y,z)+fx)/m,(Fy(x,y,z)+fy)/m,(Fz(x,y,z)+fz)/m])

t = np.linspace(0,20.0,n_iter)
s = integrate(f,s_0,t)

x,v = s[:,:N],s[:,N:]

#E = 0.5*m*np.sum(v*v,axis=1)+2.0*np.pi*G*rho_0*np.log(x[:,0]*x[:,0]+x[:,1]*x[:,1]+x[:,2]*x[:,2])
#print E

pylab.plot(x[:,0],x[:,1])
pylab.xlabel('$x$')
pylab.ylabel('$y$')
pylab.axes().set_aspect('equal','datalim')
pylab.savefig('xy.png',dpi=200)
pylab.close()

pylab.plot(t,np.sqrt(x[:,0]*x[:,0]+x[:,1]*x[:,1]))
pylab.xlabel('$t$')
pylab.ylabel('$r=\sqrt{x^{2}+y^{2}+z^{2}}$')
pylab.savefig('r.png',dpi=200)
pylab.close()
