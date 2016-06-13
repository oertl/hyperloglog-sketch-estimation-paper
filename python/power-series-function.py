#############################
# Copyright 2016 Otmar Ertl #
#############################

from numpy import ones, array
from math import expm1, cosh, tanh, sinh, log, sqrt, exp, pi, sin
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar

numPoints = 100000
xMin = -2
xMax = 2

def h2(x):
    return exp(-pow(2,x))*pow(2,x)

def h(x):

    a = h2(x)
    amin = 0.
    k = 0
    while True:
        k -= 1
        aminold = amin
        amin += h2(k+x)
        if (amin == aminold):
            break

    amax = 0.;
    k = 0
    while True:
        k += 1
        amaxold = amax
        amax += h2(k+x)
        if (amax == amaxold):
            break

    return log(2.)*(a + amin + amax)



def f1(x):
    return h(x) - 1.

def f2(x):
    return 1. - h(x)

iter = range(0,numPoints)
xValues = list([xMin + i*((xMax - xMin)/(numPoints-1)) for i in iter])
yValues = list([f1(xValues[i]) for i in iter])
plt.plot(xValues, yValues)


res1 = minimize_scalar(f1, bounds = (0,1), method='bounded')
res2 = minimize_scalar(f2, bounds = (0,1), method='bounded')

print(res1)
print(res2)

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

fig1, ax1 = plt.subplots(subplot_kw=dict(xlim=[xMin,xMax], ylim=[-1.1e-5,1.1e-5]))
fig1.set_size_inches(6, 4)
ax1.plot(xValues, yValues, color='black')
ax1.grid(b=True, which='major', color='black', linestyle='--')
ax1.ticklabel_format(axis='y', style = 'sci', useOffset=False,  scilimits=(-1e-6,1e-6))
ax1.set_xlabel(r"$x$")
ax1.set_ylabel(r"$\xi(x)-1$")
fig1.savefig('../paper/power-series-function-minus-1.svg', format='svg', dpi=600)

