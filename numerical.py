from cmath import cos, pi, sin, tan
import functions as f
import numpy as np
import math

# A. Newton's Method
def q1A():
    def function(x): return ((sin(x)/(x-1)**2)+5 *
                             (cos(0.5*(x**2))**2)+3*x+1.5).real
    def Derivative(x): return (-10*x*cos(x**2/2)*sin(x**2/2) -
                               ((2*sin(x))/(x-1)**3)+((cos(x))/(x-1)**2)+3).real
    approx = f.newton(-1, function, Derivative, 1e-6, True,100)
    

def q1B():
    def function(xy): 
        x,y = xy
        return [2*x*sin(1.5*x)-y,
                4*(x-6.7)**5-y-2]
    def Derivative(xy): 
        x,y = xy
        return [[3*x*cos(1.5*x)+2*sin(1.5*x), -1],
                [20*(x-6.7)**4,-1]]
    approx = f.newton_jacob([1.0, -2.0], function,
                            Derivative, 1e-6, False, 100)
    f.toTable(approx)
    
def q1C():
    def function(xyz): 
        x,y,z = xyz
        return [3*x-2*y-z-4,
                -x+3*y-z-8,
                x-z-2]
   
    def Derivative(xyz):
        x,y,z = xyz
        return [[3,-2,-1],
                [-1,3,-1],
                [1,0,-1]]
    
    approx = f.newton_jacob([1.0,1.0,1.0], function, Derivative, 1e-6, False, 100)
    print(approx)
    
def q1D():
    def function(x): return ((2*(x**4))-(6*(x**3))-(8*(x**2))+(24*x)).real
    def Derivative(x): return ((8*(x**3))-(18*(x**2))-(16*x)+24).real
    approx1 = f.newton(-1.5, function, Derivative, 1e-6, True, 100)
    approx2 = f.newton(1.5, function, Derivative, 1e-6, True, 100)
    approx3 = f.newton(2.7, function, Derivative, 1e-6, True, 100)
    approx4 = f.newton(-0.5, function, Derivative, 1e-6, True, 100)


# B. Euler's Method 
def q2A():
    def function(x,y): return (-4*sin((pi/6)*math.log(y))*(x-1)).real
    approx = f.euler(0,4,2,0.1,function,False)
    f.plotG2(approx,0.1)
    
# def q2B():
    # def function

    
def q2C():
    def function(x,y): return (y+x**2-math.e**(1.5*y)).real
    approx1 = f.euler(0,10,0,1e-5,function,False)
    approx2 = f.euler(0, 10, 0, 0.05, function, False)
    f.plotG2(approx1,1e-5)
    f.plotG2(approx2,0.05)
    
def q3A():
    def function(x): return tan(x-(math.pi/3))+4*cos(x**3)
    approx = f.midpoint(0, 2.5, 1e-5, function)
    f.plotG3(approx,1e-5)
    
def q3B():
    def function(x): return sin(2*x)+cos(3*x)
    approx1 = f.midpoint(0,4*math.pi,1e-5,function) 
    approx2 = f.midpoint(0, 4*math.pi, 1e-5, function)


q3A()