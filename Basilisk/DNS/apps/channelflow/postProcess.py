import numpy as np
import matplotlib.pyplot as plt
import os, sys
import pandas as pd
import figStyle as fs

sys.path.append(os.environ['HOME'] + "/twophase/2phase/apps/channelFlow/aaa")

from TBFAverage import TBFAverage

from TBFData import TBFData

startt, endt = 0,2
startTimeNums = [startt]
endTimeNums = [endt]

casePaths = ['.']
if not os.path.isfile('fl'+str(startt)+'to'+str(endt)+'.npz'):
    TBFAverage(casePaths, startTimeNums, endTimeNums, 'fl').write('fl'+str(startt)+'to'+str(endt)+'.npz')
    
# PFl = TBFData('fl'+str(startt)+'to'+str(endt)+'.npz')
# PFl.setRho(10.0)

if not os.path.isfile('fg'+str(startt)+'to'+str(endt)+'.npz'):
    TBFAverage(casePaths, startTimeNums, endTimeNums, 'fg').write('fg'+str(startt)+'to'+str(endt)+'.npz')
    
# PFg = TBFData('fg'+str(startt)+'to'+str(endt)+'.npz')
# PFg.setRho(1.0)