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
		reader = csv.reader(file, skipinitialspace=True,delimiter=' ')
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
	ax1.yaxis.grid(b=True, which='major', color='gray', linestyle='--')
	ax1.xaxis.grid(b=True, which='major', color='gray', linestyle='--')
	ax1.set_xlabel("cardinality")
	ax1.set_ylabel("relative error")


	pm3s = (1.-erf(3./sqrt(2.)))/2.
	pp3s = (1.+erf(3./sqrt(2.)))/2.
	pm2s = (1.-erf(2./sqrt(2.)))/2.
	pp2s = (1.+erf(2./sqrt(2.)))/2.
	pm1s = (1.-erf(1./sqrt(2.)))/2.
	pp1s = (1.+erf(1./sqrt(2.)))/2.

	datam3s = numpy.divide(numpy.percentile(data, pm3s*100., axis=1), cardinalities)-1.
	datap3s = numpy.divide(numpy.percentile(data, pp3s*100., axis=1), cardinalities)-1.
	datam2s = numpy.divide(numpy.percentile(data, pm2s*100., axis=1), cardinalities)-1.
	datap2s = numpy.divide(numpy.percentile(data, pp2s*100., axis=1), cardinalities)-1.
	datam1s = numpy.divide(numpy.percentile(data, pm1s*100., axis=1), cardinalities)-1.
	datap1s = numpy.divide(numpy.percentile(data, pp1s*100., axis=1), cardinalities)-1.
	data05 = numpy.divide(numpy.percentile(data, 5, axis=1), cardinalities)-1.
	data50 = numpy.divide(numpy.percentile(data, 50, axis=1), cardinalities)-1.
	data95 = numpy.divide(numpy.percentile(data, 95, axis=1), cardinalities)-1.
	ax1.plot(cardinalities, data50, color='black', linewidth=1.0, label='median')
	ax1.fill_between(cardinalities, datam3s, datap3s, facecolor='#dddddd', edgecolor='none')
	ax1.fill_between(cardinalities, datam2s, datap2s, facecolor='#bbbbbb', edgecolor='none')
	ax1.fill_between(cardinalities, datam1s, datap1s, facecolor='#999999', edgecolor='none')
	label1 = "mid" + "{:10.2f}".format(100.*(pp1s-pm1s)) + "\%"
	label2 = "mid" + "{:10.2f}".format(100.*(pp2s-pm2s)) + "\%"
	label3 = "mid" + "{:10.2f}".format(100.*(pp3s-pm3s)) + "\%"
	plt.plot([], [], color='#999999', linewidth=10, label=label1)
	plt.plot([], [], color='#bbbbbb', linewidth=10, label=label2)
	plt.plot([], [], color='#dddddd', linewidth=10, label=label3)
	ax1.legend(loc=2,prop={'size':12})
	fig1.savefig('../paper/' + figureName + '.svg', format='svg', dpi=600)
	plt.close(fig1)

createFigure(12, 20, 'flajolet_estimates', 'original_estimate')
createFigure(12, 20, 'flajolet_mid_range_estimates', 'raw_estimate')
createFigure(12, 20, 'flajolet_small_range_estimates', 'small_range_estimate')
createFigure(8, 24, 'max_likelihood_estimates', 'max_likelihood_estimate_8_24')
createFigure(12, 20, 'max_likelihood_estimates', 'max_likelihood_estimate_12_20')
createFigure(16, 16, 'max_likelihood_estimates', 'max_likelihood_estimate_16_16')
createFigure(20, 12, 'max_likelihood_estimates', 'max_likelihood_estimate_20_12')
createFigure(12, 52, 'max_likelihood_estimates', 'max_likelihood_estimate_12_52')
createFigure(12, 14, 'max_likelihood_estimates', 'max_likelihood_estimate_12_14')
