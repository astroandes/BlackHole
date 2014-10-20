import numpy as np,sympy as sp, pylab

m = 1
G = 1
N = 3

def phi(x,y,z):
    return G*m*m*(x*x+y*y+z*z)**(-0.5)

def integrate(D,x_0,v_0,dt,n_steps):
    x = x_0
    v = v_0
    for i  in range(n_steps):
        a = D(x,v)
        x += v*dt
        v += a*dt
    return x,v

def getForce(Potential):
    x,y,z = sp.symbols('x y z')
    Fx = sp.lambdify((x,y,z),-sp.diff(Potential(x,y,z),x),'numpy')
    Fy = sp.lambdify((x,y,z),-sp.diff(Potential(x,y,z),y),'numpy')
    Fz = sp.lambdify((x,y,z),-sp.diff(Potential(x,y,z),z),'numpy')
    return Fx,Fy,Fz

Fx,Fy,Fz = getForce(phi)

def D(r,v):
    x,y,z = r
    return np.array([Fx(x,y,z),Fy(x,y,z),Fz(x,y,z)])/m
