#############################
# Copyright 2016 Otmar Ertl #
#############################

from numpy import ones, array
from math import expm1, cosh, tanh, sinh, log
import matplotlib.pyplot as plt

numPoints = 100000
xMin = 1e-3
xMax = 110

def helper(x):
    if x == 0:
        return 0
    else:
        return 1 - x/expm1(x)

def helperApprox(x):
    if x < 1e-2:
        x2 = x*x
        return 0.5*x*(1. - x/6.*(1.-x2/60.*(1.-x2/42.)))
    elif x > 45:
        return 1
    else:
        x_half = x*0.5
        q_half = helperApprox(x_half)
        t = 2*(1-q_half)
        return (x_half+q_half*t)/(x_half+t)


iter = range(0,numPoints)
xValues = list([xMin + i*((xMax - xMin)/(numPoints-1)) for i in iter])
yValues = list([helperApprox(xValues[i]) for i in iter])

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

fig1, ax1 = plt.subplots(subplot_kw=dict(xlim=[0,20], ylim=[0,1.2]))
fig1.set_size_inches(6, 4)
ax1.plot(xValues, yValues, color='black')
ax1.grid(b=True, which='major', color='black', linestyle='--')
ax1.set_xlabel(r"$x$")
ax1.set_ylabel(r"$h(x)$")
fig1.savefig('../paper/helper.svg', format='svg', dpi=600)

