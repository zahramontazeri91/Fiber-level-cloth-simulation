# -*- coding: utf-8 -*-
"""
Created on Thu May 03 15:54:49 2018

@author: zahra
"""

import numpy as np
import struct

#each voxel should be larger than z-step-size (0.01) to gaurantee between two points we don't skip any box
#curves  =
def fillVol(pnts, radius, res_x, res_y, res_z):
    max_x, max_y, max_z = np.amax(pnts, axis=0)
    min_x, min_y, min_z = np.amin(pnts, axis=0)
    
    assert (min_x != max_x and min_y != max_y and min_z != max_z)
    
    sz = pnts.shape[0]
    vol = np.zeros((res_x, res_y, res_z))
    
    for i in range (0, sz):
        
        len_x = float(max_x - min_x)
        len_y = float(max_y - min_y)
        len_z = float(max_z - min_z)
        
        idx_x =  int ((float(pnts[i][0] - min_x ) / len_x )* float(res_x) )
        idx_y =  int ((float(pnts[i][1] - min_y ) / len_y )* float(res_y) )
        idx_z =  int ((float(pnts[i][2] - min_z ) / len_z )* float(res_z) )
        
        if (idx_x == res_x):
            idx_x = idx_x - 1
        if (idx_y == res_y):
            idx_y = idx_y - 1
        if (idx_z == res_z):
            idx_z = idx_z - 1
            
        vol[idx_x][idx_y][idx_z] = 1
        
        # go d distance in all 6 directions and add neighbor voxels if needed
        bottom_x = min_x + idx_x * (len_x/float(res_x))
        top_x = bottom_x + (len_x/float(res_x))
        
        bottom_y = min_y + idx_y * (len_y/float(res_y))
        top_y = bottom_y + (len_y/float(res_y))

        bottom_z = min_z + idx_z * (len_z/float(res_z))
        top_z = bottom_z + (len_z/float(res_z))
        
#        if ( (top_x - pnts[i][0]) < radius ):
#            vol[idx_x + 1][idx_y][idx_z] = 1  
#        if ( (pnts[i][0] - bottom_x) < radius ):
#            vol[idx_x - 1][idx_y][idx_z] = 1
#            
#        if ( (top_y - pnts[i][1]) < radius ):
#            vol[idx_x][idx_y+1][idx_z] = 1  
#        if ( (pnts[i][1] - bottom_y) < radius ):
#            vol[idx_x][idx_y-1][idx_z] = 1
#            
#        if ( (top_z - pnts[i][2]) < radius ):
#            vol[idx_x][idx_y][idx_z+1] = 1  
#        if ( (pnts[i][2] - bottom_z) < radius ):
#            vol[idx_x][idx_y][idx_z-1] = 1

        
    return vol


#def writeVol(vol, fn_vol):
radius = 0.1
res_x, res_y, res_z = 4,4,4
pnts = np.array([[1.2,0,0],[2,0,0],[3,2,0], [4,0,1]])
vol = fillVol(pnts, radius, res_x, res_y, res_z)



fout = open('test.vol', 'wb')
versian = 3
#name = struct.pack('<cccB','V','O','L', versian  )
fout.write("\x56\x4f\x4c\x03\x01\x00\x00\x00")
name = struct.pack('<HI', 2,1)
fout.write(name)
fout.flush()

fout.close()