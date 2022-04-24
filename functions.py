
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import make_interp_spline, BSpline


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


def midpoint(start, stop, step, f):
    integral = 0
    coords = {}
    for i in np.arange(float(start), float(stop), step):
        coords[i] = f(i)
        integral += f(((i+step)+i)/2)*step
    print("Integral is equal to: ", integral.real)
    return coords


def trapezoid(start, stop, step, f):
    integral = 0
    for i in np.arange(float(start), float(stop), step):
        integral += ((f(i)+f(i+step))/2)*step
    print("Integral is equal to: ", integral.real)


def rightpoint(start, stop, step, f):
    integral = 0
    for i in np.arange(float(start), float(stop), step):
        integral += f(i)*step
    print("Integral is equal to: ", integral.real)


def toTable(coords):
    print(" ")
    print("{:<8} {:<15} {:<10}".format(
        'Iteration', 'X', 'Y'))
    for k, v in coords.items():
        x, y = v
        print("{:<8} {:<15} {:<10}".format(k, x, y))
    print(" ")


def plotG2(dic, step):
    myList = dic.items()

    myList = sorted(myList)
    x, y = zip(*myList)

    plt.plot(x, y)
    plt.title(f"Eulers method with step size of {step}")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.show()


def plotG3(dic, step):
    myList = dic.items()

    myList = sorted(myList)
    x, y = zip(*myList)

    plt.plot(x, y)
    plt.hlines(0, 0, 2.5, color='red')
    plt.fill_between(x, y, color='green', alpha=0.5)
    plt.title(f"Eulers method with step size of {step}")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.show()

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

    # set bounds

    plt.show()
