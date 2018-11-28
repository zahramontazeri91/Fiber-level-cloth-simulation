# -*- coding: utf-8 -*-
"""
This script file read all the hairs for different frame and render it using 
mitsuba shape: hair
"""

import os

os.chdir("D:/sandbox/fiberSimulation/yarn_generation_project/scene");
dataset = "1224"
#### procedural renderings
for i in range (1,2):
    f = i*10 
    fn = []    
    for y in range (0,8):
        fn.append("../YarnGeneration/output/" + dataset + "/genYarn_" + str(f) + "_" + str(y)+ ".txt")
    os.system('D:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn0="%s" -D fn1="%s" -D fn2="%s" -D fn3="%s" -D fn4="%s" -D fn5="%s" -D fn6="%s" -D fn7="%s"  fibers_grid.xml' % (fn[0], fn[1], fn[2], fn[3], fn[4], fn[5], fn[6], fn[7] ));
    os.rename("fibers_grid.exr", '../results/' + dataset + '/proc_' + str(f) + '.exr')

### procedural renderings without compression
#for i in range (1,124):
#    f = i*10
#    filename = "../YarnGeneration/output/" + dataset + "/genYarn_wo_" + str(f) + ".txt"
#    os.system('D:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers_grid.xml' % (filename));
#    os.rename("fibers_grid.exr", '../results/' + dataset + '/proc_wo_' + str(f) + '.exr')

#### procedural renderings using NN
#for i in range (1,124):
#    f = i*10
#    filename = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + ".txt"
#    os.system('D:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers_grid.xml' % (filename));
#    os.rename("fibers_grid.exr", '../results/' + dataset + '/proc_NN_' + str(f) + '.exr')

#### Simulated renderings
#path = "D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/data/" + dataset 
#for i in range (107,309):
#    j = i * 10
#    filename = path + '/simul_frame_' + str(j) + '.txt'
#    os.system('D:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers_grid.xml' % (filename));
#    os.rename("fibers_grid.exr", "../results/" + dataset + "/simul_" + str(j) + ".exr")

