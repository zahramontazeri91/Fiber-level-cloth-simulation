# -*- coding: utf-8 -*-
"""
visulize a 2D curve
"""

import matplotlib.pyplot as plt

#fname = '../sineCurve.txt'
fname = '../centerline_yarn.txt'
plt.figure(figsize=(12,3))
with open(fname, 'r') as fin:
    n = int(fin.readline())
    pnt =[]
    for i in range (0,n):
        pnt.append(float(fin.readline()))
        

    plt.plot(pnt, label='centerline')
    plt.legend()
    plt.title("yarn curnterline with ", fontsize=14)
        
    fin.close()
    
print('done!')