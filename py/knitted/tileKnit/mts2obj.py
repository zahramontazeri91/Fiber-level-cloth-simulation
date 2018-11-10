# -*- coding: utf-8 -*-
"""
Created on Sat Nov 10 10:02:23 2018

@author: zahra
"""
import numpy as np

mtsfile = 'slipstitchrib_tiled_5_5.txt'
objfile = 'slipstitchrib_tiled_5_5.obj'

with open(mtsfile, 'r') as fin:
    with open(objfile, 'w') as fout:
        yarnNum = int ( fin.readline() )
        fin.readline() 
        for y in range (0, yarnNum):
            vrtxNum = int ( fin.readline() )
            for v in range (0, vrtxNum):
                line = fin.readline()
                val = np.array([float(x) for x in line.strip().split()]) 
                fout.writelines ('v %.8f %.8f %.8f \n' % (val[0]*0.5, val[1], val[2] ) )
            print(yarnNum, vrtxNum)
        