# -*- coding: utf-8 -*-
"""
Created on Thu May 24 14:19:45 2018

@author: zahra
"""

# This gets a BCC file from http://www.cemyuksel.com/research/yarnmodels/ and returns the centerlines in obj file

import numpy as np
import io, struct, time, sys

fn_bcc = 'openwork_trellis_pattern.bcc'
fn_obj = 'openwork_trellis_pattern.obj'
fn_txt = 'openwork_trellis_pattern.txt'
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
             
            fout_txt.writelines('%d\n' % (vrtx) )
            for c in range (0, curve):
                vrtx = struct.unpack('I', fin.read(4))[0]
                for v in range (0, vrtx):
                    x = struct.unpack('f', fin.read(4))[0]
                    y = struct.unpack('f', fin.read(4))[0]
                    z = struct.unpack('f', fin.read(4))[0]
                    fout.writelines('v %.8f %.8f %.8f\n' % (x,y,z) )
                    fout_txt.writelines('%.8f %.8f %.8f\n' % (x,y,z) )
                for v in range (1, vrtx):
                    fout.writelines('l %d %d \n' % (v, v+1) )
            
            
        
        
        