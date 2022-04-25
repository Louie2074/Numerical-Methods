import functions as f
import math
from cmath import cos,sin,tan,pi


'''Excellent version of Q2A Eulers method using Midpoint method of runge kutta'''
def q2AExcellent():
  def function(x, y): return (-4*sin((pi/6)*math.log(y))*(x-1)).real
  approx = f.midpointRungeKutta(0, 4,2,1e-5,function)
  f.plotG(approx, 0.1, "Eulers",False,None)

'''Excellent version of Q2C Eulers method using Midpoint method of runge kutta'''
def q2CExcellent():
    def function(x, y): return (y+x**2-math.e**(1.5*y)).real
    approx1 = f.midpointRungeKutta(0, 10, 0, 1e-5, function)
    approx2 = f.midpointRungeKutta(0, 10, 0, 0.05, function)
    f.plotG(approx1, 1e-5, "Eulers",False,None)
    f.plotG(approx2, 0.05, "Eulers",False,None)


'''Excellent version of Q3A Midpoint method using Trapezoid method of integration'''
def q3AExcellenceT():
    def function(x): return tan(x-(math.pi/3))+4*cos(x**3)
    approx = f.trapezoid(0, 2.5, 1e-5, function)
    f.plotG(approx, 1e-5, "Trapezoid", True, 2.5)


'''Excellent version of Q3B Midpoint method using Trapezoid method of integration'''
def q3BExcellenceT():
   def function(x): return sin(2*x)+cos(3*x)
   approx1 = f.trapezoid(0, 4*math.pi, 1e-5, function)
   approx2 = f.trapezoid(0, 4*math.pi, 1, function)
   f.plotG(approx1, 1e-5, "Trapezoid", True, 4*math.pi)

    
'''Excellent version of Q3A Midpoint method using Right Endpoint method of integration'''
def q3AExcellenceR():
    def function(x): return tan(x-(math.pi/3))+4*cos(x**3)
    approx = f.rightpoint(0, 2.5, 1e-5, function)
    f.plotG(approx, 1e-5, "Right Endpoint", True, 2.5)


'''Excellent version of Q3A Midpoint method using Right Endpoint method of integration'''
def q3BExcellenceR():
    def function(x): return sin(2*x)+cos(3*x)
    approx1 = f.rightpoint(0, 4*math.pi, 1e-5, function)
    approx2 = f.rightpoint(0, 4*math.pi, 1, function)
    f.plotG(approx1, 1e-5, "Right Endpoint", True, 4*math.pi)
