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

def createFigure(p, q, dataSource, figureName):

    filename = str(p) + '_' + str(q)

    m = pow(2, p)

    cardinalities = readData('cardinalities.dat').flatten()
    data = readData(filename + '_' + dataSource + '.dat')

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')

    yMax = 6./pow(2., p*0.5)
    fig1, ax1 = plt.subplots(subplot_kw=dict(xlim=[cardinalities[0],cardinalities[-1]], ylim=[-yMax,yMax], xscale='log', yscale='linear'))

    fig1.set_size_inches(10, 5)
    ax1.yaxis.grid(b=True, which='major', color='gray', linestyle='--',dash_capstyle='round')
    ax1.xaxis.grid(b=True, which='major', color='gray', linestyle='--',dash_capstyle='round')
    ax1.set_xlabel("cardinality")
    ax1.set_ylabel("relative error")

    pm3s = (1.-erf(3./sqrt(2.)))/2.
    pp3s = (1.+erf(3./sqrt(2.)))/2.
    pm2s = (1.-erf(2./sqrt(2.)))/2.
    pp2s = (1.+erf(2./sqrt(2.)))/2.
    pm1s = (1.-erf(1./sqrt(2.)))/2.
    pp1s = (1.+erf(1./sqrt(2.)))/2.

    datamin = numpy.divide(data.min(axis=1), cardinalities)-1.
    datam3s = numpy.divide(numpy.percentile(data, pm3s*100., axis=1), cardinalities)-1.
    datap3s = numpy.divide(numpy.percentile(data, pp3s*100., axis=1), cardinalities)-1.
    datam2s = numpy.divide(numpy.percentile(data, pm2s*100., axis=1), cardinalities)-1.
    datap2s = numpy.divide(numpy.percentile(data, pp2s*100., axis=1), cardinalities)-1.
    datam1s = numpy.divide(numpy.percentile(data, pm1s*100., axis=1), cardinalities)-1.
    datap1s = numpy.divide(numpy.percentile(data, pp1s*100., axis=1), cardinalities)-1.
    datamax = numpy.divide(data.max(axis=1), cardinalities)-1.
    datamean = numpy.divide(numpy.mean(data, axis=1), cardinalities)-1.
    datastd = numpy.divide(numpy.std(data, axis=1), cardinalities)
    datastdpos = datamean + datastd
    datastdneg = datamean - datastd
    data50 = numpy.divide(numpy.percentile(data, 50, axis=1), cardinalities)-1.

    ax1.plot(cardinalities, datamean, color='black', linewidth=0.75, label='mean', dashes=[1, 2], dash_capstyle='round')
    ax1.plot(cardinalities, datastdpos, color='black', linewidth=0.75, linestyle = 'dotted', label='stdev', dashes = [1, 4], dash_capstyle='round')
    ax1.plot(cardinalities, datastdneg, color='black', linewidth=0.75, linestyle = 'dotted', dashes = [1, 4], dash_capstyle='round')
    ax1.plot(cardinalities, data50, color='black', linewidth=0.75, linestyle = 'solid', label='median')
    ax1.fill_between(cardinalities, datam3s, datap3s, facecolor='#dddddd', edgecolor='none')
    ax1.fill_between(cardinalities, datam2s, datap2s, facecolor='#bbbbbb', edgecolor='none')
    ax1.fill_between(cardinalities, datam1s, datap1s, facecolor='#999999', edgecolor='none')
    label1 = "mid" + "{:10.2f}".format(100.*(pp1s-pm1s)) + "\%"
    label2 = "mid" + "{:10.2f}".format(100.*(pp2s-pm2s)) + "\%"
    label3 = "mid" + "{:10.2f}".format(100.*(pp3s-pm3s)) + "\%"
    plt.plot([], [], color='#999999', linewidth=10, label=label1)
    plt.plot([], [], color='#bbbbbb', linewidth=10, label=label2)
    plt.plot([], [], color='#dddddd', linewidth=10, label=label3)
    ax1.legend(loc=2,prop={'size':11},ncol=1)
    fig1.savefig('../paper/' + figureName + '.svg', format='svg', dpi=600)
    fig1.savefig('../paper/' + figureName + '.png', format='png', dpi=100)
    plt.close(fig1)

    numpy.savetxt(dataDir + figureName + ".min_error", datamin, delimiter=',')
    numpy.savetxt(dataDir + figureName + ".max_error", datamax, delimiter=',')

createFigure(12, 20, 'flajolet_estimates', 'original_estimate')
createFigure(12, 20, 'flajolet_mid_range_estimates', 'raw_estimate')
createFigure(12, 20, 'flajolet_small_range_estimates', 'small_range_estimate')

createFigure(8, 24,  'max_likelihood_estimates', 'max_likelihood_estimate_8_24')
createFigure(12, 20, 'max_likelihood_estimates', 'max_likelihood_estimate_12_20')
createFigure(16, 16, 'max_likelihood_estimates', 'max_likelihood_estimate_16_16')
createFigure(22, 10, 'max_likelihood_estimates', 'max_likelihood_estimate_22_10')
createFigure(12, 52, 'max_likelihood_estimates', 'max_likelihood_estimate_12_52')
createFigure(12, 14, 'max_likelihood_estimates', 'max_likelihood_estimate_12_14')
createFigure(12,  0, 'max_likelihood_estimates', 'max_likelihood_estimate_12_0')

createFigure(8, 24,  'corrected_raw_estimates', 'raw_corrected_estimate_8_24')
createFigure(12, 20, 'corrected_raw_estimates', 'raw_corrected_estimate_12_20')
createFigure(16, 16, 'corrected_raw_estimates', 'raw_corrected_estimate_16_16')
createFigure(22, 10, 'corrected_raw_estimates', 'raw_corrected_estimate_22_10')
createFigure(12, 52, 'corrected_raw_estimates', 'raw_corrected_estimate_12_52')
createFigure(12, 14, 'corrected_raw_estimates', 'raw_corrected_estimate_12_14')
createFigure(12, 0,  'corrected_raw_estimates', 'raw_corrected_estimate_12_0')
