#!/usr/bin/python

import sys
import numpy as np
import matplotlib as mpl

import matplotlib.pyplot as plt

import figStyle as fs

fs.prep(plt, 'square')

fig = plt.figure()

t = np.loadtxt('t.txt')
xb = np.loadtxt('xb.txt')

plt.plot(t, xb)

plt.xlabel(r'$t$')
plt.ylabel(r'$x_b$')

fs.post(fig, 'square')

plt.savefig('xb.pdf')
