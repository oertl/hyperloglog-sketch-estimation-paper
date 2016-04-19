#############################
# Copyright 2016 Otmar Ertl #
#############################

from math import expm1, ceil
from scipy.optimize import minimize

def h(x):
    if x == 0:
        return 0
    else:
        return 1 - x/expm1(x)
        
def f(x):
    return -((2*h(x)*(1-2*h(x)))/(x+2*h(x)*(1-h(x)))) - 2*h(x)/(x+2*(1-h(x)))
    
    
res = minimize(f, 1)

out = open('../paper/max_errror_propagation_factor.txt', 'w')
out.write(str(ceil(-res.fun*1000)/1000.))

