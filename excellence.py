import functions as f
import math
import docx
from cmath import cos,sin,tan,pi


# Question 2 B. Euler's Method excellence using Midpoint Method as other Runge Katta

mydoc = docx.Document("./Doc1.docx")

def q2AExcellent():
    def function(x, y): return (-4*sin((pi/6)*math.log(y))*(x-1)).real
    approx = f.midpointRK(0, 4, 2, 0.1, function)
    f.plotG2(approx, 0.1,"Midpoint")



def q2CExcellent():
    def function(x, y): return (y+x**2-math.e**(1.5*y)).real
    approx1 = f.midpointRK(0, 10, 0, 1e-5, function)
    approx2 = f.midpointRK(0, 10, 0, 0.05, function)
    f.plotG2(approx1, 1e-5, "Midpoint")
    f.plotG2(approx2, 0.05, "Midpoint")




# Question 2 C. Midpoint Method excellence

# Trapezoid method part a
def q3CExcellenceAT():
    def function(x): return tan(x-(math.pi/3))+4*cos(x**3)
    approx = f.trapezoid(0, 2.5, 1e-5, function)
    f.plotG3(approx, 1e-5, 2.5,"Trapezoid")

# Trapezoid method part b
def q3CExcellenceBT():
    def function(x): return sin(2*x)+cos(3*x)
    approx1 = f.trapezoid(0,4*math.pi,1e-5,function)
    approx2 = f.trapezoid(0, 4*math.pi, 1, function)
    f.plotG3(approx1, 1e-5, 4*math.pi, "Trapezoid")
    

# Right Endpoint method part a    
def q3CExcellenceAR():
    def function(x): return tan(x-(math.pi/3))+4*cos(x**3)
    approx = f.rightpoint(0, 2.5, 1e-5, function)
    f.plotG3(approx, 1e-5, 2.5,"Right Endpoint")

# Right Endpoint method part b
def q3CExcellenceBR():
    def function(x): return sin(2*x)+cos(3*x)
    approx1 = f.trapezoid(0, 4*math.pi, 1e-5, function)
    approx2 = f.trapezoid(0, 4*math.pi, 1, function)
    f.plotG3(approx1, 1e-5, 4*math.pi, "Right Endpoint")
    
       
def runAllExcel():
    q2AExcellent()
    q2CExcellent()
    q3CExcellenceAT()
    q3CExcellenceBT()
    q3CExcellenceAR()
    q3CExcellenceBR()

runAllExcel()
