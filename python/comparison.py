#############################
# Copyright 2016 Otmar Ertl #
#############################

from numpy import ones, array
from math import expm1, cosh, tanh, sinh, log, erf, sqrt
import matplotlib.pyplot as plt
import csv
import sys
import numpy

dataDir = '../data/'

def readData(filename):
    csv.field_size_limit(sys.maxsize)
    with open(dataDir + filename , 'r') as file:
        reader = csv.reader(file, skipinitialspace=True)
        return numpy.array([[float(e) for e in r] for r in reader][1:])

config = [
["12_20_max_likelihood_estimates.dat", "blue", False, "maximum likelihood estimator"],
#["12_20_flajolet_mid_range_estimates.dat", "gray", False, "raw estimator"],
["12_20_flajolet_mid_range_estimates.dat", "red", True, "raw estimator (bias corrected)"],
["12_20_flajolet_small_range_estimates.dat", "black", False, "linear counting estimator"],
["12_20_corrected_raw_estimates.dat", "green", False, "corrected raw estimator"]
]

cardinalities = readData('cardinalities.dat').flatten()

def applyBiasCorrection(data, cardinalities):
    datatmp = numpy.copy(data)
    for i in range(0, 100):
        data_mean = numpy.mean(datatmp, axis=1) - cardinalities
        datatmp -= numpy.interp(datatmp, cardinalities, data_mean)
    return datatmp

def calculate(filename, doBiasCorrection):

    data = readData(filename)
    if (doBiasCorrection):
        data = applyBiasCorrection(data, cardinalities)
    relative_error = data/cardinalities[:,None]-1.
    datamean = numpy.mean(relative_error, axis=1)
    datastd = numpy.std(relative_error, axis=1)
    return (datamean, datastd)

fig1, ax1 = plt.subplots(subplot_kw=dict(xlim=[cardinalities[0],cardinalities[-1]], ylim=[0, 0.02], xscale='log', yscale='linear'))

fig1.set_size_inches(10, 5)
ax1.yaxis.grid(b=True, which='major', color='gray', linestyle='--',dash_capstyle='round')
ax1.xaxis.grid(b=True, which='major', color='gray', linestyle='--',dash_capstyle='round')
ax1.set_xlabel("cardinality")
ax1.set_ylabel("stdev")

fig2, ax2 = plt.subplots(subplot_kw=dict(xlim=[cardinalities[0],cardinalities[-1]], ylim=[-0.02,0.02], xscale='log', yscale='linear'))

fig2.set_size_inches(10, 5)
ax2.yaxis.grid(b=True, which='major', color='gray', linestyle='--',dash_capstyle='round')
ax2.xaxis.grid(b=True, which='major', color='gray', linestyle='--',dash_capstyle='round')
ax2.set_xlabel("cardinality")
ax2.set_ylabel("bias")


plt.rc('text', usetex=True)
plt.rc('font', family='serif')

for c in config:
    (mean, stddev) = calculate(c[0], c[2])
    ax1.plot(cardinalities, stddev, color=c[1], linewidth=0.75, label=c[3])
    ax2.plot(cardinalities, mean, color=c[1], linewidth=0.75, label=c[3])

ax1.legend(loc=4,prop={'size':11},ncol=1)
ax2.legend(loc=4,prop={'size':11},ncol=1)

fig1.savefig('../paper/stdev_comparison.svg', format='svg', dpi=600)
fig1.savefig('../paper/stdev_comparison.png', format='png', dpi=100)
plt.close(fig1)
fig2.savefig('../paper/mean_comparison.svg', format='svg', dpi=600)
fig2.savefig('../paper/mean_comparison.png', format='png', dpi=100)
plt.close(fig2)
