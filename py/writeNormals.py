# -*- coding: utf-8 -*-
"""
write normals of the curve
"""

normfile = open('../genYarn_frame1_norms.txt','w')

vrtx_num = 1480
normfile.writelines('%d \n' % vrtx_num )  
for v in range (0, vrtx_num):
    normfile.writelines('%.6f %.6f %.6f \n' % (0.0, 1.0, 0.0) )
normfile.close()