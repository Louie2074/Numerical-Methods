
from cProfile import run
import io
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import make_interp_spline, BSpline
import docx
from cmath import cos, pi, sin, tan
import numpy as np
import math

mydoc = docx.Document("./Doc1.docx")


def q1A():
    def function(x): return ((sin(x)/(x-1)**2)+5 *
                             (cos(0.5*(x**2))**2)+3*x+1.5).real
    def Derivative(x): return (-10*x*cos(x**2/2)*sin(x**2/2) -
                               ((2*sin(x))/(x-1)**3)+((cos(x))/(x-1)**2)+3).real
    approx = newton(-1, function, Derivative, 1e-6, True, 100)


def q1B():
    def function(xy):
        x, y = xy
        return [2*x*sin(1.5*x)-y,
                4*(x-6.7)**5-y-2]

    def Derivative(xy):
        x, y = xy
        return [[3*x*cos(1.5*x)+2*sin(1.5*x), -1],
                [20*(x-6.7)**4, -1]]
    approx = newton_jacob([1.0, -2.0], function,
                            Derivative, 1e-6, False, 100)
    tab = toTable(approx)
    print(tab)
    mydoc.add_paragraph(tab)
    mydoc.save("./Doc1.docx")


def q1C():
    def function(xyz):
        x, y, z = xyz
        return [3*x-2*y-z-4,
                -x+3*y-z-8,
                x-z-2]

    def Derivative(xyz):
        x, y, z = xyz
        return [[3, -2, -1],
                [-1, 3, -1],
                [1, 0, -1]]

    approx = newton_jacob([1.0, 1.0, 1.0], function,
                            Derivative, 1e-6, False, 100)


def q1D():
    def function(x): return ((2*(x**4))-(6*(x**3))-(8*(x**2))+(24*x)).real
    def Derivative(x): return ((8*(x**3))-(18*(x**2))-(16*x)+24).real
    approx1 = newton(-1.5, function, Derivative, 1e-6, True, 100)
    approx2 = newton(1.5, function, Derivative, 1e-6, True, 100)
    approx3 = newton(2.7, function, Derivative, 1e-6, True, 100)
    approx4 = newton(-0.5, function, Derivative, 1e-6, True, 100)


# B. Euler's Method
def q2A():
    def function(x, y): return (-4*sin((pi/6)*math.log(y))*(x-1)).real
    approx = euler(0, 4, 2, 0.1, function, False)
    plotG2(approx, 0.1, "Eulers")


def q2B():
    def function(x, y, z): return -sin(z*(2*x+y))
    def function2(x, y): return math.e**(-0.05*x)*cos(4*y)
    euler2(0, 7, -1, 0, 0.01, function, function2, False)


def q2C():
    def function(x, y): return (y+x**2-math.e**(1.5*y)).real
    approx1 = euler(0, 10, 0, 1e-5, function, False)
    approx2 = euler(0, 10, 0, 0.05, function, False)
    plotG2(approx1, 1e-5, "Eulers")
    plotG2(approx2, 0.05, "Eulers")


def q3A():
    def function(x): return tan(x-(math.pi/3))+4*cos(x**3)
    approx = midpoint(0, 2.5, 1e-5, function)
    plotG3(approx, 1e-5, 2.5,"Midpoint")


def q3B():
    def function(x): return sin(2*x)+cos(3*x)
    approx1 = midpoint(0, 4*math.pi, 1e-5, function)
    approx2 = midpoint(0, 4*math.pi, 1, function)
    plotG3(approx1, 1e-5, 4*math.pi, "Midpoint")
    plotG3(approx2, 1, 4*math.pi, "Midpoint")


def newton(x0, f, d, epsilon, plot, max_iter=100):
    xn = x0
    coords = {}
    for i in range(0, max_iter):
        fxn = f(xn)
        coords[i] = [xn.real, fxn.real]
        if plot == True:
            plotG(f, d, xn, i+1)
        if abs(fxn) < epsilon:
            print(
                f'Initial Guess: {x0}, Found solution after {i+1} iterations. The answer is {xn}')
            if plot == False:
                toTable(coords)
            return coords
        Dfxn = d(xn)
        if Dfxn == 0:
            print('Derivative is 0. No solution found.')
            return coords
        xn = xn - (fxn*(Dfxn**-1))

    print('No solution found.')
    return coords


def newton_jacob(guess, f, d, epsilon, plot, max_iter=100):
    xn = guess
    coords = {}
    for i in range(max_iter):
        # Solve J(xn)*( xn+1 - xn ) = -F(xn):
        J = np.array(d(xn))
        F = np.array(f(xn))

        diff = np.linalg.solve(J, -F)
        coords[i] = [xn[0].real, xn[1].real]
        xn = xn + diff

        # Stop condition:
        if np.linalg.norm(diff) < epsilon:
            print(f'{[i.real for i in xn]}')
            break

    else:  # only if the for loop end 'naturally'
        print('not converged')

    return coords


def euler(x0, n, y0, h, f, plot):
    # Calculating step size
    yn = y0
    coords = {}
    count = 0
    for i in np.arange(float(x0), float(n), h):
        coords[i] = yn.real
        new = yn+(h*f(i, yn))
        count += 1
        yn = new
    print(
        f'Initial Guess: {x0}, Found solution after {count} iterations. The answer is {yn}')
    return coords


def euler2(x0, n, y0,z0, h, f,f2, plot):
    # Calculating step size
    yn = y0
    zn = z0
    lastX = 0
    count = 0
    for i in np.arange(float(x0), float(n), h):
      
        
        newZ = zn+(h*f2(i, zn))
       
        newY = yn+(h*f(i, yn,zn))
        count += 1
        zn = newZ
        yn = newY
        lastX = i
    print(
        f'Initial Guess: {x0}, Found solution after {count} iterations. The answer is {lastX.real,yn.real,zn.real}')
    



def midpoint(start, stop, step, f):
    integral = 0
    coords = {}
    for i in np.arange(float(start), float(stop), step):
        coords[i] = f(i).real
        integral += f(((i+step)+i)/2)*step
    print("Integral is equal to: ", integral.real)
    return coords


def midpointRK(x0, n, y0, h, f):
    yn = y0
    count = 0
    coords = {}
    for i in np.arange(float(x0), float(n), h):
        coords[i] = yn.real
        k1 = f(i,yn)
        k2 = f(i+(h/2), yn+(k1*h/2))
        count += 1
        yn = yn+(k2*h)
    print(
        f'Initial Guess: {x0}, Found solution after {count} iterations. The answer is {yn}')
    return coords


def trapezoid(start, stop, step, f):
    integral = 0
    coords = {}
    for i in np.arange(float(start), float(stop), step):
        coords[i] = f(i)
        integral += ((f(i)+f(i+step))/2)*step
    print("Integral is equal to: ", integral.real)
    return coords

def rightpoint(start, stop, step, f):
    integral = 0
    coords = {}
    for i in np.arange(float(start), float(stop), step):
        coords[i] = f(i)
        integral += f(i)*step
    print("Integral is equal to: ", integral.real)
    return coords

def toTable(coords):
    str = ""
    str+=" "
    str+="{:<8} {:<15} {:<10}".format('Iteration', 'X', 'Y')
    for k, v in coords.items():
        x, y = v
        str+="\n"
        str+="{:<8} {:<15} {:<10}".format(k, x, y)
    str+=" "
    return str


def plotG2(dic, step,name):
    myList = dic.items()

    myList = sorted(myList)
    x, y = zip(*myList)

    plt.plot(x, y)
    plt.title(f"{name} method with step size of {step}")
    plt.xlabel("X")
    plt.ylabel("Y")
    figure = io.BytesIO()
    plt.savefig(figure, format='png')
    figure.seek(0)

    mydoc.add_picture(figure)
    mydoc.save('./Doc1.docx')
    plt.close()


def plotG3(dic, step,max, name):
    myList = dic.items()

    myList = sorted(myList)
    x, y = zip(*myList)

    plt.plot(x, y)
    plt.hlines(0, 0, max, color='red')
    plt.fill_between(x, y, color='green', alpha=0.5)
    plt.title(f"{name} method with step size of {step}")
    plt.xlabel("X")
    plt.ylabel("Y")
    figure = io.BytesIO()
    plt.savefig(figure, format='png')
    figure.seek(0)

    mydoc.add_picture(figure)
    mydoc.save('./Doc1.docx')
    plt.close()

def plotG(f, d, xn, iter):
    x = np.arange(-5, 5, 0.5)
    xnew = np.linspace(x.min(), x.max(), 300)

    y = []
    for i in x:
        try:
            y.append(f(i))
        except:
            y.append(10)
    spl = make_interp_spline(x, y, k=3)

    power_smooth = spl(xnew)
    fig = plt.figure(1)
    ax = fig.add_subplot(111)

    ax.spines['left'].set_position('zero')
    ax.spines['right'].set_color('none')
    ax.spines['bottom'].set_position('zero')
    ax.spines['top'].set_color('none')
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')
    ax.plot(xnew, power_smooth, linewidth=3, color='red', label="F(x)")

    ax.axline((xn, f(xn)), slope=d(xn), color='C0', label='Tangent')
    ax.plot(xn, f(xn), marker="o", color="green",
            label=f"Guess \n X: {f'{xn:.2f}'} \n Y: {f'{f(xn):.2f}'}")
    ax.set_xbound(-8, 8)
    ax.set_ybound(-8, 8)
    plt.title(f"Current Iteration: {iter}")
    plt.legend(loc="upper left")
    figure = io.BytesIO()
    plt.savefig(figure,format='png')
    figure.seek(0)

    mydoc.add_picture(figure)
    mydoc.save('./Doc1.docx')
    plt.close()


    
def runAll():
    q1A()
    q1B()
    q1C()
    q1D()
    q2A()
    q2B()
    q2C()
    q3A()
    q3B()

runAll()