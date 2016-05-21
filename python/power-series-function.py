#############################
# Copyright 2016 Otmar Ertl #
#############################

from numpy import ones, array
from math import expm1, cosh, tanh, sinh, log, sqrt, exp, pi, sin
import matplotlib.pyplot as plt
from scipy.optimize import minimize_scalar

numPoints = 100000
xMin = 0.000
xMax = 1

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
plt.show()

