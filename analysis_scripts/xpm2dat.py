#!/usr/bin/env python 
import numpy as np 
import sys
import os

def Usage() : 
    print "Usage: %s < PixMap Matrix (.xpm) >" %sys.argv[0] 

##############################################
#  QUICK NON-COMPREHENSIVE FILE CHECKS       #
##############################################
try : 
    fileName = sys.argv[1] 
except IndexError : 
    print "ERROR: No file name supplied" 
    Usage() 
    sys.exit()

try : 
    base, ext = os.path.basename(fileName).split('.') 
except ValueError : 
    print "ERROR: Unable to parse fileName" 
    print fileName 
    sys.exit()

if not ext == 'xpm' : 
    print "ERROR: File must be an PixMap compatible matrix (.xpm file)" 
    sys.exit() 

if not os.path.isfile(fileName) : 
    print "ERROR: File %s not found"%fileName 
    sys.exit() 

outFile = "%s.dat"%base
outFile = os.path.join(os.getcwd(),outFile) 

if os.path.isfile(outFile) : 
    print "ERROR: %s already exists. Cowardly refusing to overwrite"%outFile
    sys.exit() 

##############################################
#  READ HEADER OF XPM TO GET CHAR->NUM DICT  #
##############################################
dataDict = {} 
dataStr = [] 
xRangeMin, yRangeMin = 1e10, 1e10
xRangeMax, yRangeMax = -1, -1
with open(fileName) as f : 
    for line in f : 
        if line.startswith('/*') : 
            if line.split()[1] == "x-axis:" : 
                xRange = np.array(line.split()[2:-1],dtype=float) 
                if np.min(xRange) < xRangeMin : xRangeMin = np.min(xRange) 
                if np.max(xRange) > xRangeMax : xRangeMax = np.max(xRange) 
            if line.split()[1] == "y-axis:" : 
                yRange = np.array(line.split()[2:-1],dtype=float) 
                if np.min(yRange) < yRangeMin : yRangeMin = np.min(yRange) 
                if np.max(yRange) > yRangeMax : yRangeMax = np.max(yRange) 
        elif line.startswith('static char') : 
            line = next(f) 
            dataLines = int(line.split()[1]) 
            numKeys = int(line.split()[2]) 
            numChars= int(line.split()[3][:-2]) ##slice off trailing ", from key length field

            for i in range(numKeys) : 
                line = next(f)   
                key, value = line.split()[0], line.split()[5]

                ##Strip leading " on key, and leading and trailing " on value
                key = key[1:] 
                value = value[1:-1]
                
                #print key, value 
                
                dataDict[key] = float(value ) 
        else : 
            dataStr.append(line) 

##############################################
#  CONVERT CHAR MATRIX TO NUMERIC MATRIX     #
##############################################
dataNums = []
for line in dataStr : 
    ##Strip leading " and trailing 
    line = line[1:] 
    numLine = [] 
    for char in [line[i:i+numChars] for i in range(0,len(line), numChars)] : ##list of every numChars
        #print char
        if not char.startswith('"') : 
            numLine.append(dataDict[char]) # * 10 ) ##nm -> A for GoodTuring model 
        else : 
            dataNums.append(numLine) 
            break 
dataNums = np.array(dataNums,dtype=float) 
dataNums = np.flip(dataNums, axis=0) ##xpm has time zero in bottom left of matrix. We need time zero to be top, and increase for each line

#for i in np.diagonal(dataNums) : 
#    if not i == 0 : 
#        print "Error: Non-zero element on diagonal %f" %i

##############################################
#  PRINT NUMERICA MATRIX TO FILE             #
##############################################
np.savetxt(outFile, dataNums,fmt="%0.3f", 
        header="x-range: %f\t%f\ny-range: %f\t%f"%(xRangeMin, xRangeMax, yRangeMin, yRangeMax)  ) 

