#!/bin/python

import numpy as np
import matplotlib.pyplot as plt
import os, sys
import figStyle as fs

# Prepare figure

fs.prep(plt, 'small')

# Load data

data = np.loadtxt('postProcessing/graph/100/line_q1_q2_q3_q5_q4_q6_q7_q8_q9_q10_p_k_epsilon.xy')

y = data[:,0]

q1 = data[:,1]
q2 = data[:,2]
q3 = data[:,3]
q4 = data[:,4]
q5 = data[:,5]
q6 = data[:,6]
q7 = data[:,7]
q8 = data[:,8]
q9 = data[:,9]
q10 = data[:,10]

# Plot

fig = plt.figure('features')

plt.plot(y, q1, label=r'$q_1$')
plt.plot(y, q2, label=r'$q_2$')
plt.plot(y, q3, label=r'$q_3$')
plt.plot(y, q4, label=r'$q_4$')
plt.plot(y, q5, label=r'$q_5$')
plt.plot(y, q6, label=r'$q_6$')
plt.plot(y, q7, label=r'$q_7$')
plt.plot(y, q8, label=r'$q_8$')
plt.plot(y, q9, label=r'$q_9$')
plt.plot(y, q10, label=r'$q_{10}$')

# Style/save

plt.xlabel(r'$y$')
plt.ylabel('Features')

fs.post(fig, 'small', plt.legend())

plt.savefig('features.pdf')
