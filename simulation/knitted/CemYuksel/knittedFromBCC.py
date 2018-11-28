# -*- coding: utf-8 -*-
"""
Make sure to run with python 2.7 (python3.5 has different binary-WR functions)
Created on Thu May 24 14:19:45 2018
This gets a BCC file from http://www.cemyuksel.com/research/yarnmodels/ and returns the centerlines in obj file
@author: zahra
"""

# This gets a BCC file from http://www.cemyuksel.com/research/yarnmodels/ and returns the centerlines in obj file


import io, struct

#fn_bcc = 'openwork_trellis_pattern.bcc'
#fn_obj = 'openwork_trellis_pattern.obj'
#fn_txt = 'openwork_trellis_pattern.txt'

fn_bcc = 'alien_sweater.bcc'
fn_obj = 'alien_sweater.obj'
fn_txt = 'alien_sweater.txt'

with open(fn_txt, "w") as fout_txt:
    with open(fn_obj, "w") as fout:
        with io.open(fn_bcc, mode='rb') as fin:
            assert fin.read(3) == 'BCC'
            assert fin.read(1) == '\x44'
            assert fin.read(2) == 'C0'
            assert fin.read(1) == '\x03'
            up = ord(struct.unpack('c', fin.read(1))[0])
            curve = struct.unpack('2I', fin.read(8))[0]
            pnts = struct.unpack('2I', fin.read(8))[0]
            
            for i in range (0,40):
                info = struct.unpack('c', fin.read(1))[0]
                print(info)
                
            print(curve, pnts)
            vrtxList = []
            fout_txt.writelines('%d\n' % (curve) )
            cnt = 0
            for c in range (0, curve):
                vrtx = struct.unpack('i', fin.read(4))[0]
                written_vrtx = abs(vrtx)
#                print('curve ', c, 'vrtx ', vrtx)
                v=0
                curve_x=[]
                curve_y=[]
                curve_z=[]
                while v < abs(vrtx):
                    v += 1
                    x = struct.unpack('f', fin.read(4))[0]
                    y = struct.unpack('f', fin.read(4))[0]
                    z = struct.unpack('f', fin.read(4))[0]
                    
                    #skip short curves
                    if (abs(vrtx) < 1000):
                        written_vrtx = 0
                        continue
                    #pick the middle part of sweater
                    if (x > 13 ):
                        written_vrtx -= 1
                        continue
                    if (x < -7 ):
                        written_vrtx -= 1
                        continue
                    if (z > 0 ):
                        written_vrtx -= 1
                        continue
                    if (abs(y) > 10 ):
                        written_vrtx -= 1
                        continue
                    
                    curve_x.append(x)
                    curve_y.append(y)
                    curve_z.append(z)
                    
                vrtxList.append(written_vrtx)
                print('written_vrtx ', written_vrtx, cnt)
                j = 0
                fout_txt.writelines('%d\n' % (written_vrtx) )
                while j < written_vrtx: 
                    fout.writelines('v %.8f %.8f %.8f\n' % (curve_x[j], curve_y[j], curve_z[j]) )
                    fout_txt.writelines('%.8f %.8f %.8f\n' % (curve_x[j], curve_y[j], curve_z[j]) )   
                    j +=1
                
                #write each curve separately
                fixedN = 369 #minimum vrtx for all curves
                if (written_vrtx > fixedN): 
                    fn_obj_curve = 'frame_0000000fiber_' + str(cnt).zfill(2) + '.obj'
                    fn_fe_curve = 'frame_0000000fiber_' + str(cnt).zfill(2) + '.fe'
                    with open(fn_obj_curve, "w") as fout_obj_curve:
                        with open(fn_fe_curve, "w") as fout_fe_curve:
                            j = 0
                            while j < fixedN: 
                                fout_obj_curve.writelines('v %.8f %.8f %.8f 0.0\n' % (curve_x[j], curve_y[j], curve_z[j]) ) #don't forget twist
                                if (j != 0):
                                    fout_fe_curve.writelines('0 1 0 0 0 1 1 0 0\n' ) 
                                j += 1
                    cnt+=1
                 
            #write the edges (should be out of for because vertex indices accumulated)
            v_total=1 
            for c in range (0, len(vrtxList) ):   
#                print(c, vrtxList[c])
                v = 0
                while v < abs(vrtxList[c]): 
                    if v==abs(vrtxList[c])-1: #skip connecting last vrtx of curve_i with first vrtx of curve_i+1
                        v_total += 1
                        v += 1
                        continue
                    fout.writelines('l %d %d \n' % (v_total, v_total+1) )
                    v_total += 1
                    v += 1

            # write the fe (deformation gradient defined for edges not vertices)
#            for c in range (0, len(vrtxList) ):   
##                print(c, vrtxList[c])
#                v = 0
#                while v < abs(vrtxList[c]):
                    
            
print('done!')
        
        
        