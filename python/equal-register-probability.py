#############################
# Copyright 2016 Otmar Ertl #
#############################

from numpy import ones, array
from math import expm1, cosh, tanh, sinh, log
import matplotlib.pyplot as plt

numPoints = 100000
xMin = 0
xMax = 1

alphaInf = 1/(2.*log(2.))

def upperBound(d):
    return 1. + 2 * alphaInf * log(1 - d/2 + d*d/16)

def lowerBound(d):
    return 1. + 2 * alphaInf * log(1 - d/2)


iter = range(0,numPoints)
xValues = list([xMin + i*((xMax - xMin)/(numPoints-1)) for i in iter])
yValuesUpperBound = list([upperBound(xValues[i]) for i in iter])
yValuesLowerBound = list([lowerBound(xValues[i]) for i in iter])

plt.rc('text', usetex=True)
plt.rc('font', family='serif')

fig1, ax1 = plt.subplots(subplot_kw=dict(xlim=[0,1], ylim=[0,1]))
fig1.set_size_inches(6, 4)
ax1.plot(xValues, yValuesUpperBound, color='black')
ax1.plot(xValues, yValuesLowerBound, color='black')
ax1.fill_between(xValues, yValuesLowerBound, yValuesUpperBound, facecolor='#dddddd', edgecolor='none')
ax1.grid(b=True, which='major', color='black', linestyle='--')
ax1.set_xlabel(r"Jaccard distance $D$")
ax1.set_ylabel(r"probability of equal registers")
fig1.savefig('../paper/equal-register-probability.svg', format='svg', dpi=600)

