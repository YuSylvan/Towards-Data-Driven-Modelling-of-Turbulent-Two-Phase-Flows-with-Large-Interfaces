import numpy as np
import matplotlib.pyplot as plt
import os, sys
import figStyle as fs

sys.path.append(os.environ['HOME'] + "/git/NRGBasilisk/apps/pipeFlow")

from TBFAverage import TBFAverage
from TBFData import TBFData
from KasagiData import KasagiData

z1 = 12
z2 = 16

casePaths = ['.']

startTimeNums = [0]
endTimeNums = [6]

# Load data

Kasagi = KasagiData(os.environ['HOME'] + '/git/NRGBasilisk/apps/pipeFlow/Kasagi.txt')

if not os.path.isfile('fl.npz'):
    TBFAverage(casePaths, startTimeNums, endTimeNums, 'fl').write('fl.npz')

if not os.path.isfile('fg.npz'):
    TBFAverage(casePaths, startTimeNums, endTimeNums, 'fg').write('fg.npz')

L = TBFData('fl.npz')
G = TBFData('fg.npz')

U = L.Re/L.Re_tau

# Plot data

fs.prep(plt, 'square')

##

fig = plt.figure('plot')

L.setRho(1.0)
G.setRho(0.1)

L.setRange(z1,z2)
G.setRange(z1,z2)

PHASE = L

# plt.plot(PHASE.r(), PHASE.fua(), label='mean u')

plt.plot(PHASE.r(), PHASE.production()*U**3, label='prod')
plt.plot(PHASE.r(), PHASE.transport()*U**3, label='trans')
plt.plot(PHASE.r(), PHASE.pressureDiffusion(1)*U**3, label='pdiff 1')
# plt.plot(PHASE.r(), PHASE.pressureDiffusion(2)*U**3, label='pdiff 2')
plt.plot(PHASE.r(), PHASE.viscousDiffusion(1)*U**3, label='visc 1')
# plt.plot(PHASE.r(), PHASE.viscousDiffusion(2)*U**3, label='visc 2')
plt.plot(PHASE.r(), -PHASE.dissipation()*U**3, label='diss')
plt.plot(PHASE.r(), PHASE.convection()*U**3, label='conv')
# plt.plot(PHASE.r(), PHASE.interface1()*U**3, label='int 1')
# plt.plot(PHASE.r(), PHASE.interface2()*U**3, label='int 2')
# plt.plot(PHASE.r(), PHASE.surfaceTension()*U**3, label='acc. s')
# plt.plot(PHASE.r(), PHASE.budgetSum(1,1)*U**3, label='sum')
plt.plot(PHASE.r(), PHASE.budgetSum(1,1)*U**3, label='sum with pdiff 1')
# plt.plot(PHASE.r(), PHASE.budgetSum(2,2)*U**3, label='sum with pdiff 2')

# plt.plot(Kasagi.r_c(), Kasagi.ua_c(), 'k')
# plt.plot(Kasagi.r_c(), Kasagi.production_c(), 'k')
# plt.plot(Kasagi.r_c(), Kasagi.transport_c(), 'k')
# plt.plot(Kasagi.r_c(), Kasagi.pressureDiffusion_c(), 'k')
# plt.plot(Kasagi.r_c(), Kasagi.viscousDiffusion_c(), 'k')
# plt.plot(Kasagi.r_c(), Kasagi.dissipation_c(), 'k')

# Plot layout

fig = plt.figure('plot')

plt.xlabel(r'$r^+$')
fs.post(fig, 'square', plt.legend())
plt.savefig('plot.pdf')

