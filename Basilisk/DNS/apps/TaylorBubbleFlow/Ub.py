#!/usr/bin/python

import sys
import numpy as np
import matplotlib as mpl

import matplotlib.pyplot as plt

import figStyle as fs

minVal = 0
maxVal = 2

fs.prep(plt, 'square')

fig = plt.figure()

t = np.loadtxt('t.txt')
Ub = np.loadtxt('Ub.txt')

Ub = np.minimum(np.maximum(Ub,minVal), maxVal)

plt.plot(t, Ub)

plt.xlabel(r'$t$')
plt.ylabel(r'$U_b$')

fs.post(fig, 'square')

plt.savefig('Ub.pdf')
