import numpy as np
import matplotlib.pyplot as plt
import os, sys
import pandas as pd
import figStyle as fs

sys.path.append(os.environ['HOME'] + "/twophase/2phase/apps/channelFlow")

from TBFAverage import TBFAverage

from TBFData import TBFData

startt, endt = 8,9
startTimeNums = [startt]
endTimeNums = [endt]

casePaths = ['.']

PFl1 = TBFData('fl'+str(startt)+'to'+str(endt)+'.npz')
PFl1.setRho(10.0)  
PFg1 = TBFData('fg'+str(startt)+'to'+str(endt)+'.npz')
PFg1.setRho(1.0)

startt, endt = 0,2
startTimeNums = [startt]
endTimeNums = [endt]

casePaths = ['.']

PFl2 = TBFData('fl'+str(startt)+'to'+str(endt)+'.npz')
PFl2.setRho(10.0)  
PFg2 = TBFData('fg'+str(startt)+'to'+str(endt)+'.npz')
PFg2.setRho(1.0)

u_openfoam=np.array(pd.read_csv('./data/E05.csv')['U:0'])
x1 = np.linspace(-0.05,0.05,len(u_openfoam))


plt.plot( u_openfoam,x1, label='openFoam')
plt.plot( PFl1.fux()+PFg1.fux(),PFl1.y(), label='9-10')
plt.plot( PFl2.fux()+PFg2.fux(),PFl2.y(), label='11-12')

plt.title('mean u_x')
plt.xlabel('u_x(m/s)')
plt.ylabel('coordy(m)')
plt.legend()
plt.savefig('com ux')