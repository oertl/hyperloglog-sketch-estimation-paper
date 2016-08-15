#############################
# Copyright 2016 Otmar Ertl #
#############################

from math import expm1, cosh, tanh, sinh, log
import matplotlib.pyplot as plt
import csv
import sys
import numpy

numPoints = 100000
xMin = 1e-3
xMax = 110
filename_12_20 = '../data/benchmark_results_12_20.dat'
filename_12_52 = '../data/benchmark_results_12_52.dat'

def readData(filename):
    csv.field_size_limit(sys.maxsize)
    with open(filename , 'r') as file:
        reader = csv.reader(file, skipinitialspace=True)
        return numpy.array([[float(e) for e in r] for r in reader])


plt.rc('text', usetex=True)
plt.rc('font', family='serif')


data_12_20 = readData(filename_12_20)
data_12_52 = readData(filename_12_52)

cardinalities_12_20 = data_12_20[:,0]
max_likelihood_avg_exec_time_nanos_12_20 = data_12_20[:,1]
max_likelihood_avg_outer_loop_count_12_20 = data_12_20[:,2]
max_likelihood_avg_inner_loop1_count_12_20 = data_12_20[:,3]
max_likelihood_avg_inner_loop2_count_12_20 = data_12_20[:,4]
max_likelihood_avg_log_eval_count_12_20 = data_12_20[:,5]
corrected_raw_1_avg_exec_time_nanos_12_20 = data_12_20[:,6]
corrected_raw_small_iter_12_20 = data_12_20[:,7]
corrected_raw_small_iter_12_20 = data_12_20[:,8]
corrected_raw_2_avg_exec_time_nanos_12_20 = data_12_20[:,9]

cardinalities_12_52 = data_12_52[:,0]
max_likelihood_avg_exec_time_nanos_12_52 = data_12_52[:,1]
max_likelihood_avg_outer_loop_count_12_52 = data_12_52[:,2]
max_likelihood_avg_inner_loop1_count_12_52 = data_12_52[:,3]
max_likelihood_avg_inner_loop2_count_12_52 = data_12_52[:,4]
max_likelihood_avg_log_eval_count_12_52 = data_12_52[:,5]
corrected_raw_1_avg_exec_time_nanos_12_52 = data_12_52[:,6]
corrected_raw_small_iter_12_52 = data_12_52[:,7]
corrected_raw_small_iter_12_52 = data_12_52[:,8]
corrected_raw_2_avg_exec_time_nanos_12_52 = data_12_52[:,9]

xMin = 1
xMax = min(cardinalities_12_20[-1], cardinalities_12_52[-1])

# max likelihood exec times
fig_time, ax_time = plt.subplots(subplot_kw=dict(xlim=[xMin,xMax], ylim=[0,800], xscale='log', yscale='linear'))
fig_time.set_size_inches(10, 5)
ax_time.yaxis.grid(b=True, which='major', color='gray', linestyle='--')
ax_time.xaxis.grid(b=True, which='major', color='gray', linestyle='--')
ax_time.set_xlabel("cardinality")
ax_time.set_ylabel("time (ns)")

ax_time.plot(cardinalities_12_20, max_likelihood_avg_exec_time_nanos_12_20, color='green', linestyle='-',linewidth=1.0, label='p=12, q=20')
ax_time.plot(cardinalities_12_52, max_likelihood_avg_exec_time_nanos_12_52, color='blue', linestyle='-',linewidth=1.0, label='p=12, q=52')

ax_time.legend(loc=2,prop={'size':12})

fig_time.savefig('../paper/max_likelihood_avg_exec_time.png', format='png', dpi=100)
fig_time.savefig('../paper/max_likelihood_avg_exec_time.svg', format='svg', dpi=600)
plt.close(fig_time)

# corrected raw exec times
fig_time, ax_time = plt.subplots(subplot_kw=dict(xlim=[xMin,xMax], ylim=[0,400], xscale='log', yscale='linear'))
fig_time.set_size_inches(10, 5)
ax_time.yaxis.grid(b=True, which='major', color='gray', linestyle='--')
ax_time.xaxis.grid(b=True, which='major', color='gray', linestyle='--')
ax_time.set_xlabel("cardinality")
ax_time.set_ylabel("time (ns)")

ax_time.plot(cardinalities_12_20, corrected_raw_1_avg_exec_time_nanos_12_20, color='green', linestyle='-',linewidth=1.0, label='p=12, q=20, on-demand')
ax_time.plot(cardinalities_12_52, corrected_raw_1_avg_exec_time_nanos_12_52, color='blue', linestyle='-',linewidth=1.0, label='p=12, q=52, on-demand')
ax_time.plot(cardinalities_12_20, corrected_raw_2_avg_exec_time_nanos_12_20, color='red', linestyle='-',linewidth=1.0, label='p=12, q=20, precalculated')
ax_time.plot(cardinalities_12_52, corrected_raw_2_avg_exec_time_nanos_12_52, color='black', linestyle='-',linewidth=1.0, label='p=12, q=52, precalculated')

ax_time.legend(loc=2,prop={'size':12})

fig_time.savefig('../paper/corrected_raw_avg_exec_time.png', format='png', dpi=100)
fig_time.savefig('../paper/corrected_raw_avg_exec_time.svg', format='svg', dpi=600)
plt.close(fig_time)


# max likelihood outer loop iterations
fig_out, ax_out = plt.subplots(subplot_kw=dict(xlim=[xMin,xMax], ylim=[0,4], xscale='log', yscale='linear'))
fig_out.set_size_inches(10, 5)
ax_out.yaxis.grid(b=True, which='major', color='gray', linestyle='--')
ax_out.xaxis.grid(b=True, which='major', color='gray', linestyle='--')
ax_out.set_xlabel("cardinality")
ax_out.set_ylabel("outer loop iterations")

ax_out.plot(cardinalities_12_20, max_likelihood_avg_outer_loop_count_12_20, color='green', linestyle='-',linewidth=1.0, label='p=12, q=20')
ax_out.plot(cardinalities_12_52, max_likelihood_avg_outer_loop_count_12_52, color='blue', linestyle='-',linewidth=1.0, label='p=12, q=52')

ax_out.legend(loc=2,prop={'size':12})

fig_out.savefig('../paper/max_likelihood_avg_outer_loop_iterations.svg', format='svg', dpi=600)
plt.close(fig_out)

# max likelihood inner loop 1 iterations
fig_out, ax_out = plt.subplots(subplot_kw=dict(xlim=[xMin,xMax], ylim=[0,10], xscale='log', yscale='linear'))
fig_out.set_size_inches(10, 5)
ax_out.yaxis.grid(b=True, which='major', color='gray', linestyle='--')
ax_out.xaxis.grid(b=True, which='major', color='gray', linestyle='--')
ax_out.set_xlabel("cardinality")
ax_out.set_ylabel("inner loop 1 iterations")

ax_out.plot(cardinalities_12_20, max_likelihood_avg_inner_loop1_count_12_20, color='green', linestyle='-',linewidth=1.0, label='p=12, q=20')
ax_out.plot(cardinalities_12_52, max_likelihood_avg_inner_loop1_count_12_52, color='blue', linestyle='-',linewidth=1.0, label='p=12, q=52')

ax_out.legend(loc=2,prop={'size':12})

fig_out.savefig('../paper/max_likelihood_avg_inner_loop_1_iterations.svg', format='svg', dpi=600)
plt.close(fig_out)

# max likelihood inner loop 2 iterations
fig_out, ax_out = plt.subplots(subplot_kw=dict(xlim=[xMin,xMax], ylim=[0,60], xscale='log', yscale='linear'))
fig_out.set_size_inches(10, 5)
ax_out.yaxis.grid(b=True, which='major', color='gray', linestyle='--')
ax_out.xaxis.grid(b=True, which='major', color='gray', linestyle='--')
ax_out.set_xlabel("cardinality")
ax_out.set_ylabel("inner loop 2 iterations")

ax_out.plot(cardinalities_12_20, max_likelihood_avg_inner_loop2_count_12_20, color='green', linestyle='-',linewidth=1.0, label='p=12, q=20')
ax_out.plot(cardinalities_12_52, max_likelihood_avg_inner_loop2_count_12_52, color='blue', linestyle='-',linewidth=1.0, label='p=12, q=52')

ax_out.legend(loc=2,prop={'size':12})

fig_out.savefig('../paper/max_likelihood_avg_inner_loop_2_iterations.svg', format='svg', dpi=600)
plt.close(fig_out)

# max likelihood num log evaluations
fig_out, ax_out = plt.subplots(subplot_kw=dict(xlim=[xMin,xMax], ylim=[0,2], xscale='log', yscale='linear'))
fig_out.set_size_inches(10, 5)
ax_out.yaxis.grid(b=True, which='major', color='gray', linestyle='--')
ax_out.xaxis.grid(b=True, which='major', color='gray', linestyle='--')
ax_out.set_xlabel("cardinality")
ax_out.set_ylabel("number of logarithm evaluations")

ax_out.plot(cardinalities_12_20, max_likelihood_avg_log_eval_count_12_20, color='green', linestyle='-',linewidth=1.0, label='p=12, q=20')
ax_out.plot(cardinalities_12_52, max_likelihood_avg_log_eval_count_12_52, color='blue', linestyle='-',linewidth=1.0, label='p=12, q=52')

ax_out.legend(loc=2,prop={'size':12})

fig_out.savefig('../paper/max_likelihood_avg_number_log_evaluations.svg', format='svg', dpi=600)
plt.close(fig_out)
