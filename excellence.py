import functions as f
import math
from cmath import cos,sin,tan

# Question 2 B. Euler's Method excellence













# Question 2 C. Midpoint Method excellence

# Trapezoid method part a
def q3CExcellenceAT():
    def function(x): return tan(x-(math.pi/3))+4*cos(x**3)
    approx = f.trapezoid(0, 2.5, 1e-5, function)

# Trapezoid method part b
def q3CExcellenceBT():
    def function(x): return sin(2*x)+cos(3*x)
    approx = f.trapezoid(0,4*math.pi,1e-5,function) 

# Right Endpoint method part a    
def q3CExcellenceAR():
    def function(x): return tan(x-(math.pi/3))+4*cos(x**3)
    approx = f.rightpoint(0, 2.5, 1e-5, function)

# Right Endpoint method part b
def q3CExcellenceBR():
    def function(x): return sin(2*x)+cos(3*x)
    approx = f.rightpoint(0,4*math.pi,1e-5,function) 
       
