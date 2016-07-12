#############################
# Copyright 2016 Otmar Ertl #
#############################

from math import expm1, ceil
from scipy.optimize import minimize

epsBound = 0.1

def h(x):
    if x == 0:
        return 0
    else:
        return 1 - x/expm1(x)

def f1(x):
    return -((h(2*x)*(1-2*h(2*x)))/(x+h(2*x)*(1-h(2*x)))) - h(2*x)/(x+1-h(2*x))

def f2(x):
    return - h(2*x)*h(2*x)/(x+h(2*x)*(1-h(2*x)))

def f3(x):
    return - h(2*x)/(x+1-h(2*x))

res1 = minimize(f1, 1)
res2 = minimize(f2, 1)
res3 = minimize(f3, 1)

factor1 = ceil(-res1.fun*1000)/1000.
factor2 = ceil(-res2.fun*1000)/1000.
factor3 = ceil(-res3.fun*1000)/1000.
factor4 = ceil((factor1 + epsBound*factor2)/(1. - factor3*epsBound)*1000.)/1000.

out1 = open('../paper/max_errror_propagation_factor1.txt', 'w')
out1.write(str(factor1))
out2 = open('../paper/max_errror_propagation_factor2.txt', 'w')
out2.write(str(factor2))
out3 = open('../paper/max_errror_propagation_factor3.txt', 'w')
out3.write(str(factor3))
out4 = open('../paper/max_errror_propagation_factor4.txt', 'w')
out4.write(str(factor4))
