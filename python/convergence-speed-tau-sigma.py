#############################
# Copyright 2016 Otmar Ertl #
#############################

from math import exp, sqrt, log

def tau(x, eps):
    numIterations = 0
    if x == 0. or x == 1.:
        return (0., numIterations)
    y = 1
    z = 1 - x
    while True:
        numIterations += 1
        x = sqrt(x)
        zPrime = z
        y *= 0.5
        z -= pow(1-x,2)*y
        if (zPrime - zPrime * eps <= z):
            return (z, numIterations)

def sigma(x, eps):
    numIterations = 0
    if (x == 1.):
        return (float('inf'), numIterations)
    y = 1
    z = x
    while True:
        numIterations += 1
        x *= x
        zPrime = z
        z += x * y
        y += y
        if (zPrime + zPrime * eps >= z):
            return (z, numIterations)

for p in range(1, 31):
    s  = sigma(1. - pow(2.,-p), 0)
    t  = tau(pow(2.,-p), 0)
    print(str(p) + " " + str(s) + " " + str(t))
