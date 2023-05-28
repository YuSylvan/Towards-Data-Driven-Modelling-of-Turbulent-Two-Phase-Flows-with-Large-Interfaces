# Class to read TBF cases and ensemble-average the data

import numpy as np
import re, os, sys
import pandas as pd

class TBFAverage:

    # Static data

    # Number of data files, named <phase name>-<file num>.txt

    NFiles = 3

    # Read case

    def readCase(self):

        # Read times

        timeNames = [ f.name for f in os.scandir(self.casePath) if (f.is_dir() and re.match('[0-9]+', f.name)) ]
        timeNums = np.array([ float(t) for t in timeNames ])

        # Sort

        timeNames = [ timeNames[i] for i in timeNums.argsort() ]
        timeNums = [ timeNums[i] for i in timeNums.argsort() ]

        # Select required time range and store

        timeNames = timeNames[self.startTimeNum:self.endTimeNum+1]
        timeNums = timeNums[self.startTimeNum:self.endTimeNum+1]

        self.timeNames = timeNames
        self.timeNums = timeNums

        Nt = len(timeNums)
        print(timeNames)
        print(timeNums)

        # Load raw data and ensemble-average

        rawData = self.readRawData(0)
        

        y = rawData[:,0]

        self.y = np.unique(y)

        data = rawData[:,1:]

        for i in range(1,Nt):

            rawData = self.readRawData(i)

            data += rawData[:,1:]

        data = data/float(Nt)
        

        # Reduce spatial data

        Ny = len(self.y)
        

        if len(y) != Ny:

            print('Mismatch between Nr, Nz and the number of rows in the data')
            sys.exit(1)

        Nf = data.shape[1]

        data = data.reshape((Ny, Nf))

        if hasattr(self, 'data'):

            if self.data.shape != (Ny, Nf):

                print('Mismatch in data size across cases')
                sys.exit(1)

            else:

                NtNew = self.Nt + Nt

                self.data = (self.Nt*self.data + Nt*data)/NtNew
                self.Nt = NtNew

        else:

            self.data = data
            self.Nt = Nt

    # Read cases

    def readCases(self):

        for (self.casePath, self.startTimeNum, self.endTimeNum) in \
            zip(self.casePaths, self.startTimeNums, self.endTimeNums):

            self.readCase()

    # Construct the path of the data file

    def dataFilePath(self, i):

        return self.casePath+'/'+self.timeNames[i]

    # Read the data file and return the raw data

    def readRawData(self, i):

        print('Reading ' + self.phaseName + ' fields from', self.dataFilePath(i))

        filePath = self.dataFilePath(i) + '/' + self.phaseName + 'averages.txt'

        if os.path.isfile(filePath):

            # Legacy format in which all data is contained in one file

            return pd.read_csv(filePath, delim_whitespace=True, skiprows=None, header=None).values[:,:]

        else:

            # New format in which fields are stored to separte files. The first
            # file (void fraction) also contains the coordinates of the points
            # in the first two collumns.

            for j in range (1,self.NFiles+1):

                filePath = self.dataFilePath(i) + '/' + self.phaseName + '-' + str(j) + '.txt'

                print(filePath)

                datai = pd.read_csv(filePath, delim_whitespace=True, skiprows=None, header=None).values[:,:]

                if 'data' not in locals():
                    data = datai
                else:
                    data = np.hstack((data,datai))

            return data

    # Write the averaged data

    def write(self, fileName):

        np.savez(fileName, y=self.y, data=self.data)

    # Initialize and read data

    def __init__(self, casePaths, startTimeNums, endTimeNums, phaseName):

        self.casePaths = casePaths

        self.startTimeNums = startTimeNums
        self.endTimeNums = endTimeNums

        self.phaseName = phaseName

        if  len(casePaths) != len(startTimeNums) or \
            len(casePaths) != len(endTimeNums):

            print('Provided case path, start and end time number arrays do not have equal length')
            sys.exit(1)

        self.readCases()
