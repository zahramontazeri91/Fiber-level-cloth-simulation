import os 

os.chdir('F:/sandbox/fiberSimulation/yarn_generation_project/scene')
# In[]:
"""
return cylinder position given simulation parameters
nonlinear move
"""
    ###
    # read cylinder center from mesh
#    fn = "F:/yarn11_yarn_1.6/frame_" + str(f).zfill(7) + "distance_field.obj"
#    x_sum = 0
#    y_sum = 0
#    z_sum = 0
#    n = 1695
#    with open(fn, 'r') as fin:  
#        for i in range (0, n):
#            info = fin.readline().split()
#            v = info[0]
#            x = float(info[1])
#            y = float(info[2])
#            z = float(info[3])
#            x_sum = x_sum + x
#            y_sum = y_sum + y
#            z_sum = z_sum + z
#            
#    x_sum = x_sum/n
#    y_sum = y_sum/n
#    z_sum = z_sum/n    
##    print (x_sum/4, y_sum/4, z_sum/4)
#    pos = y_sum/4
    ###
    
def cubic_ease_function( t, t0, t1, L):
    ta=t0+(t1-t0)/3 
    tb=t1-(t1-t0)/3
    yh = (L * 2.0) / (t1 - t0 + tb - ta)
    if(t < t0 or t > t1):
        return 0.0
    else:
        if(t < ta):
            return (yh * (t0 - t) * (t0 - t) * (t0 - 3.0 * ta + 2.0 * t)) / ((t0 - ta) * (t0 - ta) * (t0 - ta))
        elif (t > tb):
            return (yh * (t1 - t) * (t1 - t) * (t1 - 3.0 * tb + 2.0 * t)) / ((t1 - tb) * (t1 - tb) * (t1 - tb))
        else:
            return yh
        
t0 = 0.0
t1 = 1.0
L = 0.8
firstFrame = 0
lastFrame = 20000
skipFactor = 500
dt = (t1-t0)/( (lastFrame-firstFrame)/skipFactor )
pos = 0
t = t0
for f in range (firstFrame, lastFrame+1, skipFactor):
    t = t + dt
    yh = cubic_ease_function( t, t0, t1, L)
    pos = pos + yh*dt
    print('frame:',f,'t', t, 'position', pos)
    
# In[]:
"""
return cylinder position given simulation parameters
Linear movement (doesn't match)
"""    
init_y = 0.8/4.0 #initial cylinder position
final_y = 0
durationTime = 1.0
dt = 0.00004
output_step = 500
totalFrame = durationTime / dt
outpuFrame = totalFrame / output_step
dy = abs(init_y - final_y) / outpuFrame

for f in range (firstFrame, lastFrame+1, output_step):
    h1 = init_y - ( dy * int(f/output_step) )
    h2 = -h1
    print (h1, h2)
# In[]:
#dataset = "pattern/yarn8/spacing1.0x/00011"
#dataset = 'twist/yarn4/damp2_500'
#dataset = 'pattern/yarn4/spacing0.5x/10/Raymond'
#dataset = 'yarn_stretch'
#dataset = 'single_yarn/yarn8/stretch'
#dataset =  'woven/yarn9/8x8'
#dataset = 'single_yarn/yarn4/teeth/4_1.2'
dataset = 'single_yarn/yarn8/teeth/4_1.6'
#dataset = 'single_yarn/yarn4/teeth/4_1.2_00110'
#dataset = 'woven/arbitrary_pattern/100x100'
#dataset = 'woven/arbitrary_pattern/100x100'

skipFactor = 1000
yarnNum = 1

for i in range (12000/skipFactor, 17000/skipFactor+1):
    f = i*skipFactor
    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
#        fn = "../YarnGeneration/genYarn.txt"
        print(fn)
        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
        os.rename("fibers.exr", '../results/' + dataset + '/NN_' + str(f) + '_' + str(y) + '_flicker_bcsdf.exr')


# In[]:

#dataset = 'single_yarn/yarn11/teeth/4_1.6'
dataset = 'single_yarn/yarn8/teeth/4_1.2'
#dataset = 'single_yarn/yarn4/teeth/4_1.2_00110'
#dataset = 'single_yarn/yarn8/stretch'
#dataset = 'single_yarn/yarn11/stretch'
#dataset = 'single_yarn/yarn4/teeth/1.2_110'
#dataset = 'pattern/yarn100/spacing2.5x/10'
skipFactor = 1000
yarnNum = 1
frame0 = 7000
frame1 = 20000

for i in range (frame0/skipFactor,frame1/skipFactor+1):
    f = i*skipFactor
    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"

        step = 0.03 #always
        h1 = 0.25 - (i)*step
        h2 = -h1
        
        print(h1, h2, fn)
        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D h1=%.2f -D h2=%.2f -D fn="%s" fibers.xml' % (h1, h2, fn) )
        os.rename("fibers.exr", '../results/' + dataset + '/PlyReg2_' + str(f) + '.exr')

#for i in range (frame0/skipFactor,frame1/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        step = 0.03 #always
#        h1 = 0.25 - (i)*step
#        h2 = -h1
#        
#        print(h1, h2, fn)
#        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ "_us.txt"
##        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ "_us.txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn) )
#        os.rename("fibers.exr", '../results/' + dataset + '/new_NN_' + str(f) + '_khodaya_teaser.exr')
#        

#for i in range (frame0/skipFactor,frame1/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
##        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
##        fn = "../YarnGeneration/output/" + dataset + "/genYarn_wo_" + str(f) + "_" + str(y)+ ".txt"
#        fn = "../YarnGeneration/output/" + dataset + "/centerlines/curve_" + str(f) + "_" + str(y)+ "_ds_mitsuba.txt"
#        
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D h1=0 -D h2=0 -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_wo_' + str(f) + '_' + str(y) + '_centerline_blue.exr')
#    
################################################################
#for i in range (frame0/skipFactor,frame1/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
##        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_wo_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D h1=0 -D h2=0 -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_wo_' + str(f) + '_' + str(y) + '_figure2.exr')

#for i in range (frame/skipFactor,frame/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
##        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_fly_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_fly_' + str(f) + '_' + str(y) + '.exr')


#for i in range (frame/skipFactor,frame/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
##        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_wo_fly_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_wo_fly_' + str(f) + '_' + str(y) + '.exr')
#

###############
# In[]:
dataset = 'woven/release/yarn9/8x8'

##### procedural renderings
skipFactor = 1000
yarnNum = 16

for i in range (0/skipFactor, 0/skipFactor+1):
    f = i*skipFactor
    fn_a = []
    for y in range(0,yarnNum):
#        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
        fn_a.append(fn)
        print(fn)
        

    os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba \
              -D fn0="%s" -D fn1="%s" -D fn2="%s" -D fn3="%s" -D fn4="%s" -D fn5="%s" -D fn6="%s" -D fn7="%s" -D fn8="%s" -D fn9="%s" \
              -D fn10="%s" -D fn11="%s" -D fn12="%s" -D fn13="%s" -D fn14="%s" -D fn15="%s" \
              woven_8x8.xml'% (fn_a[0], fn_a[1], fn_a[2], fn_a[3], fn_a[4], fn_a[5], fn_a[6], fn_a[7], fn_a[8], fn_a[9], \
              fn_a[10], fn_a[11], fn_a[12], fn_a[13], fn_a[14], fn_a[15] ))

#    os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn0="%s" woven_8x8.xml'% (fn_a[0] ))        
    os.rename("woven_8x8.exr", '../results/' + dataset + '/proc_NN_' + str(f) + '_x_1_z_1.exr')
        
# In[]:
#dataset = 'single_yarn/yarn11/teeth/4_1.6'
#dataset = 'single_yarn/yarn8/teeth/4_1.2_00110'
#dataset = 'single_yarn/yarn9/teeth/4_1.2'
#dataset = 'single_yarn/yarn100'
dataset = 'single_yarn/yarn4/stretch'
#dataset = 'single_yarn/yarn100/stretch'
f0 = 0
f1 = 0
skipfactor = 1000
#
for f in range (f0,f1+1, skipfactor):
    #fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_0.txt"
    #os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba tiled_cylinder.xml -D fn="%s"' %fn )    
    #os.rename("tiled_cylinder.exr", '../results/' + dataset + '/simul_tiled_' + str(f) + '_missing.exr')
    #
    fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_0.txt"
    os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D h1=0 -D h2=0 -D fn="%s" tiled_cylinder.xml' %fn )    
    os.rename("tiled_cylinder.exr", '../results/' + dataset + '/proc_NN_tiled_' + str(f) + '_khodaya.exr')
    
    #fn = "../YarnGeneration/output/" + dataset + "/genYarn_wo_" + str(f) + "_0.txt"
    #print(fn)
    #os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba tiled_cylinder.xml -D fn="%s"' %fn )    
    #os.rename("tiled_cylinder.exr", '../results/' + dataset + '/proc_wo_tiled_' + str(f) + '_missing.exr')


# In[]:
dataset = 'woven/6x6'

f0 = 6000
f1 = 6000
#for f in range (f0, f1+1, 1000):
#    os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -b 4 woven/6x6_' +str(f) +'.xml')  
##    os.system('F:/sandbox/fiberSimulation/mitsuba_ct2/mitsuba-ct2/dist/mitsuba -b 8 -D testinfo="F:/sandbox/fiberSimulation/yarn_generation_project/scene/woven/test_info_0.txt" -D obj="F:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/output/woven/push/yarn4/100x100/volume_" + str(f) +".obj" fibers_ct2.xml')
#    os.rename("woven/6x6_"+str(f) +".exr", '../results/' + dataset + '/simul_new_' + str(f) + '_relaxed_us.exr')
#

##### procedural renderings
#f0 = 1000000
#f1 = 1000000
#for f in range (f0, f1+1, 1000):
#    os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -b 8 woven/6x6_' +str(f) +'.xml')  
##    os.system('F:/sandbox/fiberSimulation/mitsuba_ct2/mitsuba-ct2/dist/mitsuba -b 8 -D testinfo="F:/sandbox/fiberSimulation/yarn_generation_project/scene/woven/test_info_0.txt" -D obj="F:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/output/woven/push/yarn4/100x100/volume_" + str(f) +".obj" fibers_ct2.xml')
#    os.rename("woven/6x6_"+str(f) +".exr", '../results/' + dataset + '/simul_' + str(f) + '_final_sample10_2.exr')

#f0 = 15000
#f1 = 15000
#for f in range (f0, f1+1, 1000):
#    os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -b 8 woven/6x6_' +str(f) +'.xml')  
##    os.system('F:/sandbox/fiberSimulation/mitsuba_ct2/mitsuba-ct2/dist/mitsuba -b 8 -D testinfo="F:/sandbox/fiberSimulation/yarn_generation_project/scene/woven/test_info_0.txt" -D obj="F:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/output/woven/push/yarn4/100x100/volume_" + str(f) +".obj" fibers_ct2.xml')
#    os.rename("woven/6x6_"+str(f) +".exr", '../results/' + dataset + '/proc_NN_' + str(f) + '_final_sample10_2.exr')
#    
#f0 = 10000
#f1 = 10000
#for f in range (f0, f1+1, 1000):
#    os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -b 8 woven/6x6_' +str(f) +'.xml')  
##    os.system('F:/sandbox/fiberSimulation/mitsuba_ct2/mitsuba-ct2/dist/mitsuba -b 8 -D testinfo="F:/sandbox/fiberSimulation/yarn_generation_project/scene/woven/test_info_0.txt" -D obj="F:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/output/woven/push/yarn4/100x100/volume_" + str(f) +".obj" fibers_ct2.xml')
#    os.rename("woven/6x6_"+str(f) +".exr", '../results/' + dataset + '/simul_' + str(f) + '_final_sample10_2.exr')
########## tiled 6x6:
os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -b 8 woven/simul_frame_6000_full_khodaya.xml')  
os.rename("woven/simul_frame_10000_full.exr", '../results/' + dataset + '/simul_6000_tiled.exr')

os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -b 8 woven/genYarn_NN_wo_6000_full_khodaya.xml')  
os.rename("woven/genYarn_NN_15000_full.exr", '../results/' + dataset + '/NN_wo_6000_tiled.exr')
   
  
# In[]:
dataset = 'woven/push/yarn8/100x100'
#dataset = 'woven/stretch/yarn4/100x100'
#dataset = 'woven/arbitrary_pattern/100x100'
#dataset = 'woven/arbitrary_pattern/512x512'

##### procedural renderings
f0 = 0
f1 = 0
for f in range (f0, f1+1, 200):
    testinfo = "../YarnGeneration/output/"+dataset+"/testinfo_" + str(f) +".txt"
    obj = "../YarnGeneration/output/"+dataset+"/volume_" + str(f) + ".obj"
    print(testinfo)
    print(obj)
    os.system('F:/sandbox/fiberSimulation/mitsuba_ct2/mitsuba-ct2/dist/mitsuba fibers_ct2.xml' )
#    os.system('F:/sandbox/fiberSimulation/mitsuba_ct2/mitsuba-ct2/dist/mitsuba fibers_ct2.xml')
    os.rename("fibers_ct2.exr", '../results/' + dataset + '/proc_NN_fly_' + str(f) + '.exr')
        
###############
# In[]:
#dataset = "woven/yarn4/spacing1.0x/00011"
#dataset = "fall/yarn4/fall"
#dataset = "ball_fall"
#dataset = "shear/shear0403"
#dataset = 'woven/yarn4/spacing1.0x/00011/shear'
dataset = 'woven/release/yarn4/100x100'

##### procedural renderings
skipFactor = 1000
yarnNum = 46

for i in range (1000/skipFactor, 1000/skipFactor+1):
    f = i*skipFactor
    fn_a = []
    for y in range(0,yarnNum):
        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_fly_" + str(f) + "_" + str(y)+ ".txt"
        fn_a.append(fn)
        

    os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba \
              -D fn0="%s" -D fn1="%s" -D fn2="%s" -D fn3="%s" -D fn4="%s" -D fn5="%s" -D fn6="%s" -D fn7="%s" -D fn8="%s" -D fn9="%s" \
              -D fn10="%s" -D fn11="%s" -D fn12="%s" -D fn13="%s" -D fn14="%s" -D fn15="%s" -D fn16="%s" -D fn17="%s" -D fn18="%s" -D fn19="%s" \
              -D fn20="%s" -D fn21="%s" -D fn22="%s" -D fn23="%s" -D fn24="%s" -D fn25="%s" -D fn26="%s" -D fn27="%s" -D fn28="%s" -D fn29="%s" \
              -D fn30="%s" -D fn31="%s" -D fn32="%s" -D fn33="%s" -D fn34="%s" -D fn35="%s" -D fn36="%s" -D fn37="%s" -D fn38="%s" -D fn39="%s" \
              -D fn40="%s" -D fn41="%s" -D fn42="%s" -D fn43="%s" -D fn44="%s" -D fn45="%s"\
              fibers_grid.xml'% (fn_a[0], fn_a[1], fn_a[2], fn_a[3], fn_a[4], fn_a[5], fn_a[6], fn_a[7], fn_a[8], fn_a[9], \
              fn_a[10], fn_a[11], fn_a[12], fn_a[13], fn_a[14], fn_a[15], fn_a[16], fn_a[17], fn_a[18], fn_a[19], \
              fn_a[20], fn_a[21], fn_a[22], fn_a[23], fn_a[24], fn_a[25], fn_a[26], fn_a[27], fn_a[28], fn_a[29], \
              fn_a[30], fn_a[31], fn_a[32], fn_a[33], fn_a[34], fn_a[35], fn_a[36], fn_a[37], fn_a[38], fn_a[39], \
              fn_a[40], fn_a[41], fn_a[42], fn_a[43], fn_a[44], fn_a[45] ))

##    os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn0="%s" fibers_grid.xml'% (fn_a[0] ))        
    os.rename("fibers_grid.exr", '../results/' + dataset + '/proc_NN_' + str(f) + '.exr')
        
        
        
###################
# In[]:
#dataset = "stretch/yarn4/stretch"
#
### procedural renderings
#skipFactor = 1000
#yarnNum = 46
#for i in range (64000/skipFactor,70000/skipFactor+1):
#    fn_a = []
#    f = i*skipFactor
#    for y in range (0, yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_wo_" + str(f) + "_" + str(y)+ ".txt"
#        fn_a.append(fn)
#        
#
#    os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba \
#              -D fn0="%s" -D fn1="%s" -D fn2="%s" -D fn3="%s" -D fn4="%s" -D fn5="%s" -D fn6="%s" -D fn7="%s" -D fn8="%s" -D fn9="%s" \
#              -D fn10="%s" -D fn11="%s" -D fn12="%s" -D fn13="%s" -D fn14="%s" -D fn15="%s" -D fn16="%s" -D fn17="%s" -D fn18="%s" -D fn19="%s" \
#              -D fn20="%s" -D fn21="%s" -D fn22="%s" -D fn23="%s" -D fn24="%s" -D fn25="%s" -D fn26="%s" -D fn27="%s" -D fn28="%s" -D fn29="%s" \
#              -D fn30="%s" -D fn31="%s" -D fn32="%s" -D fn33="%s" -D fn34="%s" -D fn35="%s" -D fn36="%s" -D fn37="%s" -D fn38="%s" -D fn39="%s" \
#              -D fn40="%s" -D fn41="%s" -D fn42="%s" -D fn43="%s" -D fn44="%s" -D fn45="%s"\
#              fibers_grid.xml'% (fn_a[0], fn_a[1], fn_a[2], fn_a[3], fn_a[4], fn_a[5], fn_a[6], fn_a[7], fn_a[8], fn_a[9], \
#              fn_a[10], fn_a[11], fn_a[12], fn_a[13], fn_a[14], fn_a[15], fn_a[16], fn_a[17], fn_a[18], fn_a[19], \
#              fn_a[20], fn_a[21], fn_a[22], fn_a[23], fn_a[24], fn_a[25], fn_a[26], fn_a[27], fn_a[28], fn_a[29], \
#              fn_a[30], fn_a[31], fn_a[32], fn_a[33], fn_a[34], fn_a[35], fn_a[36], fn_a[37], fn_a[38], fn_a[39], \
#              fn_a[40], fn_a[41], fn_a[42], fn_a[43], fn_a[44], fn_a[45] ))
#
##    fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(0)+ ".txt"
##    os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn0="%s" fibers_grid.xml' % (fn));
#    os.rename("fibers_grid.exr", '../results/' + dataset + '/proc_NN_wo_' + str(f) + '.exr')
     
# In[]:
## procedural renderings
#skipFactor = 200
#yarn0 = 0
#yarn1 = 1
#dataset = "speedtest/1.0_both"
#for i in range (0/skipFactor,15000/skipFactor+1):
#    f = i*skipFactor
#    for y in range (yarn0, yarn1):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_' + str(f) + '_' + str(y) + '.exr')
#dataset = "speedtest/0411/0.3"
#for i in range (0/skipFactor,10000/skipFactor+1):
#    f = i*skipFactor
#    for y in range (yarn0, yarn1):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_' + str(f) + '_' + str(y) + '.exr')
#dataset = "speedtest/0411/0.4"
#for i in range (0/skipFactor,10000/skipFactor+1):
#    f = i*skipFactor
#    for y in range (yarn0, yarn1):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_' + str(f) + '_' + str(y) + '.exr')
#dataset = "speedtest/0411/0.5"
#for i in range (0/skipFactor,10000/skipFactor+1):
#    f = i*skipFactor
#    for y in range (yarn0, yarn1):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_' + str(f) + '_' + str(y) + '.exr')
#dataset = "speedtest/0411/0.6"
#for i in range (0/skipFactor,10000/skipFactor+1):
#    f = i*skipFactor
#    for y in range (yarn0, yarn1):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_' + str(f) + '_' + str(y) + '.exr')

# In[]:
#dataset = "stretch/yarn4/stretch"
#
### procedural renderings
#skipFactor = 500
#yarn0 = 10
#yarn1 = 11
#for i in range (10000/skipFactor,74000/skipFactor+1):
#    
#    f = i*skipFactor
#    for y in range (yarn0, yarn1):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_' + str(f) + '_' + str(y) + '.exr')

# In[]:
#skipFactor = 2000
#dataset = "twist/yarn4/damp2_500"
#yarnNum = 1
#for i in range (24000/skipFactor,99500/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_fly_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_' + str(f/100) + '.exr')
# In[]:
##################################
        #### patterns ####
###############
## In[]:
#dataset = "pattern/yarn4/spacing0.5x/10"
#yarnNum = 1
#for i in range (12000/skipFactor,14000/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_' + str(f/100) + '_9_batch32.exr')
#        
#dataset = "pattern/yarn4/spacing0.5x/00011"
#yarnNum = 1
#for i in range (14000/skipFactor,16000/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_fly_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_fly_' + str(f/100) + '.exr')
#        
#dataset = "pattern/yarn4/spacing0.5x/10100"
#yarnNum = 1
#for i in range (13000/skipFactor,15000/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_fly_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_fly_' + str(f/100) + '.exr')
#        
#dataset = "pattern/yarn4/spacing0.5x/11110"
#yarnNum = 1
#for i in range (13000/skipFactor,15000/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_fly_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_fly_' + str(f/100) + '.exr')
#        
#        
###############
# In[]:
#dataset = "pattern/yarn4/spacing1.0x/10"
#yarnNum = 1
#for i in range (12500/skipFactor,14500/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_' + str(f/100) + '_9_batch32.exr')
#        
#dataset = "pattern/yarn4/spacing1.0x/00011"
#yarnNum = 1
#for i in range (15000/skipFactor,17000/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_fly_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_fly_' + str(f/100) + '.exr')
#        
#dataset = "pattern/yarn4/spacing1.0x/10100"
#yarnNum = 1
#for i in range (13500/skipFactor,15500/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_fly_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_fly_' + str(f/100) + '.exr')
#        
#dataset = "pattern/yarn4/spacing1.0x/11110"
#yarnNum = 1
#for i in range (14000/skipFactor,16000/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_fly_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_fly_' + str(f/100) + '.exr')
#        
#        
###############
## In[]:
#dataset = "pattern/yarn4/spacing1.5x/10"
#yarnNum = 1
#for i in range (13000/skipFactor,15000/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_' + str(f/100) + '_9_batch16_pattern10.exr')
#        
#dataset = "pattern/yarn4/spacing1.5x/00011"
#yarnNum = 1
#for i in range (15500/skipFactor,17500/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_fly_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_fly_' + str(f/100) + '.exr')
#        
#dataset = "pattern/yarn4/spacing1.5x/10100"
#yarnNum = 1
#for i in range (14000/skipFactor,16000/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_fly_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_fly_' + str(f/100) + '.exr')
#        
#dataset = "pattern/yarn4/spacing1.5x/11110"
#yarnNum = 1
#for i in range (14500/skipFactor,16500/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_fly_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_fly_' + str(f/100) + '.exr')
#        
#        
#################
# In[]:
#dataset = "pattern/yarn4/spacing2.0x/10"
#yarnNum = 1
#for i in range (14000/skipFactor,16000/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_' + str(f/100) + '_9.exr')
#        
#dataset = "pattern/yarn4/spacing2.0x/00011"
#yarnNum = 1
#for i in range (16000/skipFactor,18000/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_fly_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_fly_' + str(f/100) + '.exr')
#        
#dataset = "pattern/yarn4/spacing2.0x/10100"
#yarnNum = 1
#for i in range (15000/skipFactor,17000/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_fly_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_fly_' + str(f/100) + '.exr')
#        
#dataset = "pattern/yarn4/spacing2.0x/11110"
#yarnNum = 1
#for i in range (15500/skipFactor,17500/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_fly_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_fly_' + str(f/100) + '.exr')
#        
################
# In[]:
#dataset = "pattern/yarn4/spacing2.5x/10"
#yarnNum = 1
#for i in range (14500/skipFactor,16500/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_' + str(f/100) + '_9.exr')
#        
#dataset = "pattern/yarn4/spacing2.5x/00011"
#yarnNum = 1
#for i in range (16500/skipFactor,18500/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_fly_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_fly_' + str(f/100) + '.exr')
#        
#dataset = "pattern/yarn4/spacing2.5x/10100"
#yarnNum = 1
#for i in range (15500/skipFactor,17500/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_fly_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_fly_' + str(f/100) + '.exr')
#        
#dataset = "pattern/yarn4/spacing2.5x/11110"
#yarnNum = 1
#for i in range (16500/skipFactor,18500/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_fly_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_fly_' + str(f/100) + '.exr')
#
################
# In[]:
#dataset = "pattern/yarn4/spacing3.0x/10"
#yarnNum = 1
#for i in range (14500/skipFactor,16500/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_' + str(f/100) + '_9.exr')
#        
#dataset = "pattern/yarn4/spacing3.0x/00011"
#yarnNum = 1
#for i in range (17000/skipFactor,19000/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_fly_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_fly_' + str(f/100) + '.exr')
#        
#dataset = "pattern/yarn4/spacing3.0x/10100"
#yarnNum = 1
#for i in range (18000/skipFactor,20000/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_fly_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_fly_' + str(f/100) + '.exr')
#        
#dataset = "pattern/yarn4/spacing3.0x/11110"
#yarnNum = 1
#for i in range (17000/skipFactor,19000/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_fly_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_fly_' + str(f/100) + '.exr')
#
#
#
#
#
# In[]:
# In[]:
# In[]:
# In[]:
##################################
        #### patterns #### - simulate (yarn4)
###############
# In[]:
##################################
        #### patterns ####
###############
## In[]:
#dataset = "pattern/yarn11/spacing0.5x/10"
#yarnNum = 1
#for i in range (12000/skipFactor,14000/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
#        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/simul_' + str(f/100) + '.exr')
#        
#dataset = "pattern/yarn11/spacing0.5x/00011"
#yarnNum = 1
#for i in range (14000/skipFactor,16000/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_fly_" + str(f) + "_" + str(y)+ ".txt"
#        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/simul_' + str(f/100) + '.exr')
#        
#dataset = "pattern/yarn11/spacing0.5x/10100"
#yarnNum = 1
#for i in range (13000/skipFactor,15000/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_fly_" + str(f) + "_" + str(y)+ ".txt"
#        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/simul_' + str(f/100) + '.exr')
#        
#dataset = "pattern/yarn11/spacing0.5x/11110"
#yarnNum = 1
#for i in range (13000/skipFactor,15000/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_fly_" + str(f) + "_" + str(y)+ ".txt"
#        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/simul_' + str(f/100) + '.exr')
#        
#        
###############
## In[]:
#dataset = "pattern/yarn11/spacing1.0x/10"
#yarnNum = 1
#for i in range (12500/skipFactor,14500/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
#        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/simul_' + str(f/100) + '.exr')
#        
#dataset = "pattern/yarn11/spacing1.0x/00011"
#yarnNum = 1
#for i in range (15000/skipFactor,17000/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_fly_" + str(f) + "_" + str(y)+ ".txt"
#        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/simul_' + str(f/100) + '.exr')
#        
#dataset = "pattern/yarn11/spacing1.0x/10100"
#yarnNum = 1
#for i in range (13500/skipFactor,15500/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_fly_" + str(f) + "_" + str(y)+ ".txt"
#        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/simul_' + str(f/100) + '.exr')
#        
#dataset = "pattern/yarn11/spacing1.0x/11110"
#yarnNum = 1
#for i in range (14000/skipFactor,16000/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_fly_" + str(f) + "_" + str(y)+ ".txt"
#        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/simul_' + str(f/100) + '.exr')
#        
#        
###############
## In[]:
#dataset = "pattern/yarn11/spacing1.5x/10"
#yarnNum = 1
#for i in range (13000/skipFactor,15000/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
#        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/simul_' + str(f/100) + '.exr')
#        
#dataset = "pattern/yarn11/spacing1.5x/00011"
#yarnNum = 1
#for i in range (15500/skipFactor,17500/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_fly_" + str(f) + "_" + str(y)+ ".txt"
#        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/simul_' + str(f/100) + '.exr')
#        
#dataset = "pattern/yarn11/spacing1.5x/10100"
#yarnNum = 1
#for i in range (14000/skipFactor,16000/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_fly_" + str(f) + "_" + str(y)+ ".txt"
#        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/simul_' + str(f/100) + '.exr')
#        
#dataset = "pattern/yarn11/spacing1.5x/11110"
#yarnNum = 1
#for i in range (14500/skipFactor,16500/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_fly_" + str(f) + "_" + str(y)+ ".txt"
#        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/simul_' + str(f/100) + '.exr')
#        
        
################
## In[]:
#dataset = "pattern/yarn11/spacing2.0x/10"
#yarnNum = 1
#for i in range (14000/skipFactor,16000/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
#        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/simul_' + str(f/100) + '_9.exr')
#        
#dataset = "pattern/yarn11/spacing2.0x/00011"
#yarnNum = 1
#for i in range (16000/skipFactor,18000/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_fly_" + str(f) + "_" + str(y)+ ".txt"
#        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/simul_' + str(f/100) + '.exr')
#        
#dataset = "pattern/yarn11/spacing2.0x/10100"
#yarnNum = 1
#for i in range (15000/skipFactor,17000/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_fly_" + str(f) + "_" + str(y)+ ".txt"
#        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/simul_' + str(f/100) + '.exr')
#        
#dataset = "pattern/yarn11/spacing2.0x/11110"
#yarnNum = 1
#for i in range (15500/skipFactor,17500/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_fly_" + str(f) + "_" + str(y)+ ".txt"
#        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/simul_' + str(f/100) + '.exr')
#        
################
## In[]:
#dataset = "pattern/yarn11/spacing2.5x/10"
#yarnNum = 1
#for i in range (14500/skipFactor,16500/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
#        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/simul_' + str(f/100) + '.exr')
#        
#dataset = "pattern/yarn11/spacing2.5x/00011"
#yarnNum = 1
#for i in range (16500/skipFactor,18500/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_fly_" + str(f) + "_" + str(y)+ ".txt"
#        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/simul_' + str(f/100) + '.exr')
#        
#dataset = "pattern/yarn11/spacing2.5x/10100"
#yarnNum = 1
#for i in range (15500/skipFactor,17500/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_fly_" + str(f) + "_" + str(y)+ ".txt"
#        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/simul_' + str(f/100) + '.exr')
#        
#dataset = "pattern/yarn11/spacing2.5x/11110"
#yarnNum = 1
#for i in range (16500/skipFactor,18500/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_fly_" + str(f) + "_" + str(y)+ ".txt"
#        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/simul_' + str(f/100) + '.exr')

###############
# In[]:
#dataset = "pattern/yarn11/spacing3.0x/10"
#yarnNum = 1
#for i in range (14500/skipFactor,16500/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
#        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/simul_' + str(f/100) + '.exr')
#        
#dataset = "pattern/yarn11/spacing3.0x/00011"
#yarnNum = 1
#for i in range (17000/skipFactor,19000/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_fly_" + str(f) + "_" + str(y)+ ".txt"
#        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/simul_' + str(f/100) + '.exr')
#        
#dataset = "pattern/yarn11/spacing3.0x/10100"
#yarnNum = 1
#for i in range (18000/skipFactor,20000/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_fly_" + str(f) + "_" + str(y)+ ".txt"
#        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/simul_' + str(f/100) + '.exr')
#        
#dataset = "pattern/yarn11/spacing3.0x/11110"
#yarnNum = 1
#for i in range (17000/skipFactor,19000/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_fly_" + str(f) + "_" + str(y)+ ".txt"
#        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/simul_' + str(f/100) + '.exr')


# In[]:
# In[]:
# In[]:
# In[]:
##################################
        #### patterns #### - procedural (yarn11)
###############
# In[]:
#dataset = "pattern/yarn11/spacing0.5x/10"
#yarnNum = 1
#for i in range (12000/skipFactor,14000/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
##        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_' + str(f/100) + '.exr')
#        
#dataset = "pattern/yarn11/spacing0.5x/00011"
#yarnNum = 1
#for i in range (14000/skipFactor,16000/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
##        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_' + str(f/100) + '.exr')
#        
#dataset = "pattern/yarn11/spacing0.5x/10100"
#yarnNum = 1
#for i in range (13000/skipFactor,15000/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
##        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_' + str(f/100) + '.exr')
#        
#dataset = "pattern/yarn11/spacing0.5x/11110"
#yarnNum = 1
#for i in range (13000/skipFactor,15000/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
##        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_' + str(f/100) + '.exr')
#        
#        
###############
## In[]:
#dataset = "pattern/yarn11/spacing1.0x/10"
#yarnNum = 1
#for i in range (12500/skipFactor,14500/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
##        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_' + str(f/100) + '.exr')
#        
#dataset = "pattern/yarn11/spacing1.0x/00011"
#yarnNum = 1
#for i in range (15000/skipFactor,17000/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
##        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_' + str(f/100) + '.exr')
#        
#dataset = "pattern/yarn11/spacing1.0x/10100"
#yarnNum = 1
#for i in range (13500/skipFactor,15500/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
##        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_' + str(f/100) + '.exr')
#        
#dataset = "pattern/yarn11/spacing1.0x/11110"
#yarnNum = 1
#for i in range (14000/skipFactor,16000/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
##        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_' + str(f/100) + '.exr')
#        
#        
###############
## In[]:
#dataset = "pattern/yarn11/spacing1.5x/10"
#yarnNum = 1
#for i in range (13000/skipFactor,15000/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
##        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_' + str(f/100) + '.exr')
#        
#dataset = "pattern/yarn11/spacing1.5x/00011"
#yarnNum = 1
#for i in range (15500/skipFactor,17500/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
##        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_' + str(f/100) + '.exr')
#        
#dataset = "pattern/yarn11/spacing1.5x/10100"
#yarnNum = 1
#for i in range (14000/skipFactor,16000/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
##        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_' + str(f/100) + '.exr')
#        
#dataset = "pattern/yarn11/spacing1.5x/11110"
#yarnNum = 1
#for i in range (14500/skipFactor,16500/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
##        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_' + str(f/100) + '.exr')
#        
#        
#################
## In[]:
#dataset = "pattern/yarn11/spacing2.0x/10"
#yarnNum = 1
#for i in range (14000/skipFactor,16000/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
##        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_' + str(f/100) + '_9.exr')
#        
#dataset = "pattern/yarn11/spacing2.0x/00011"
#yarnNum = 1
#for i in range (16000/skipFactor,18000/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
##        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_' + str(f/100) + '.exr')
#        
#dataset = "pattern/yarn11/spacing2.0x/10100"
#yarnNum = 1
#for i in range (15000/skipFactor,17000/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
##        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_' + str(f/100) + '.exr')
#        
#dataset = "pattern/yarn11/spacing2.0x/11110"
#yarnNum = 1
#for i in range (15500/skipFactor,17500/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
##        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_' + str(f/100) + '.exr')
#        
################
## In[]:
#dataset = "pattern/yarn11/spacing2.5x/10"
#yarnNum = 1
#for i in range (14500/skipFactor,16500/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
##        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_' + str(f/100) + '.exr')
#        
#dataset = "pattern/yarn11/spacing2.5x/00011"
#yarnNum = 1
#for i in range (16500/skipFactor,18500/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
##        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_' + str(f/100) + '.exr')
#        
#dataset = "pattern/yarn11/spacing2.5x/10100"
#yarnNum = 1
#for i in range (15500/skipFactor,17500/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
##        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_' + str(f/100) + '.exr')
#        
#dataset = "pattern/yarn11/spacing2.5x/11110"
#yarnNum = 1
#for i in range (16500/skipFactor,18500/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
##        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_' + str(f/100) + '.exr')
#
################
## In[]:
#dataset = "pattern/yarn11/spacing3.0x/10"
#yarnNum = 1
#for i in range (14500/skipFactor,16500/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
##        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_' + str(f/100) + '.exr')
#        
#dataset = "pattern/yarn11/spacing3.0x/00011"
#yarnNum = 1
#for i in range (17000/skipFactor,19000/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
##        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_' + str(f/100) + '.exr')
#        
#dataset = "pattern/yarn11/spacing3.0x/10100"
#yarnNum = 1
#for i in range (18000/skipFactor,20000/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
##        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_' + str(f/100) + '.exr')
#        
#dataset = "pattern/yarn11/spacing3.0x/11110"
#yarnNum = 1
#for i in range (17000/skipFactor,19000/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
##        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_NN_' + str(f/100) + '.exr')

##################################
# In[]:
dataset = 'single_yarn/yarn4/teeth/4_1.6'
#dataset = 'single_yarn/yarn4/teeth/4_1.2'
#dataset = 'single_yarn/yarn4/teeth/4_1.2_00110'
#dataset = 'single_yarn/yarn100/stretch'
#dataset = 'single_yarn/yarn11/stretch'
#dataset = 'single_yarn/yarn4/teeth/1.2_110'
#dataset = 'pattern/yarn100/spacing2.5x/10'
skipFactor = 1000
yarnNum = 1
frame0 = 12000
frame1 = 12000

#for i in range (frame0/skipFactor,frame1/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
##        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
##        step = 0.0297
#        step = 0.03 #1.6
#        h1 = (0.25) - i*step #8 is number of frames
#        h2 = -h1
#        print(h1, h2, fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D h1=0 -D h2=0 -D fn="%s" fibers_crop.xml' % (fn));
#        os.rename("fibers_crop.exr", '../results/' + dataset + '/simul_' + str(f) + '_' + str(y) + '_cropped2.exr')
#
#for i in range (frame0/skipFactor,frame1/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
##        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ "_us.txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D h1=0.005 -D h2=0.005 -D fn="%s" fibers_crop.xml' % (fn));
#        os.rename("fibers_crop.exr", '../results/' + dataset + '/proc_NN_' + str(f) + '_' + str(y) + '_cropped2.exr')
##        
#
#for i in range (frame0/skipFactor,frame1/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
##        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_wo_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D h1=0 -D h2=0 -D fn="%s" fibers_crop.xml' % (fn));
#        os.rename("fibers_crop.exr", '../results/' + dataset + '/proc_wo_' + str(f) + '_' + str(y) + '_cropped2.exr')

dataset = 'single_yarn/yarn8/teeth/4_1.6'
#dataset = 'single_yarn/yarn8/teeth/4_1.2'
#dataset = 'single_yarn/yarn8/teeth/4_1.2_00110'
#dataset = 'single_yarn/yarn100/stretch'
#dataset = 'single_yarn/yarn11/stretch'
#dataset = 'single_yarn/yarn4/teeth/1.2_110'
#dataset = 'pattern/yarn100/spacing2.5x/10'
skipFactor = 1000
yarnNum = 1
frame0 = 12000
frame1 = 12000

#for i in range (frame0/skipFactor,frame1/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
##        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
##        step = 0.0297
#        step = 0.03 #1.6
#        h1 = (0.25) - i*step #8 is number of frames
#        h2 = -h1
#        print(h1, h2, fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D h1=0 -D h2=0 -D fn="%s" fibers_crop.xml' % (fn));
#        os.rename("fibers_crop.exr", '../results/' + dataset + '/simul_' + str(f) + '_' + str(y) + '_cropped3.exr')

for i in range (frame0/skipFactor,frame1/skipFactor+1):
    f = i*skipFactor
    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ "_us.txt"
        print(fn)
        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D h1=0.005 -D h2=0.005 -D fn="%s" fibers_crop.xml' % (fn));
        os.rename("fibers_crop.exr", '../results/' + dataset + '/proc_NN_' + str(f) + '_' + str(y) + '_cropped3.exr')
#        

#for i in range (frame0/skipFactor,frame1/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
##        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_wo_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D h1=0 -D h2=0 -D fn="%s" fibers_crop.xml' % (fn));
#        os.rename("fibers_crop.exr", '../results/' + dataset + '/proc_wo_' + str(f) + '_' + str(y) + '_cropped3.exr')
##
#dataset = 'single_yarn/yarn11/teeth/4_1.6'
##dataset = 'single_yarn/yarn11/teeth/4_1.2'
##dataset = 'single_yarn/yarn11/teeth/4_1.2_00110'
##dataset = 'single_yarn/yarn100/stretch'
##dataset = 'single_yarn/yarn11/stretch'
##dataset = 'single_yarn/yarn4/teeth/1.2_110'
##dataset = 'pattern/yarn100/spacing2.5x/10'
#skipFactor = 1000
#yarnNum = 1
#frame0 = 12000
#frame1 = 12000
#
#for i in range (frame0/skipFactor,frame1/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
##        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
##        step = 0.0297
#        step = 0.03 #1.6
#        h1 = (0.25) - i*step #8 is number of frames
#        h2 = -h1
#        print(h1, h2, fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D h1=0 -D h2=0 -D fn="%s" fibers_crop.xml' % (fn));
#        os.rename("fibers_crop.exr", '../results/' + dataset + '/simul_' + str(f) + '_' + str(y) + '_cropped3.exr')
#
#for i in range (frame0/skipFactor,frame1/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
##        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D h1=0.005 -D h2=0.005 -D fn="%s" fibers_crop.xml' % (fn));
#        os.rename("fibers_crop.exr", '../results/' + dataset + '/proc_NN_' + str(f) + '_' + str(y) + '_cropped3.exr')
##        
#
#for i in range (frame0/skipFactor,frame1/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
##        fn = "../YarnGeneration/data/" + dataset + "/simul_frame_" + str(f) + "_" + str(y)+ ".txt"
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_wo_" + str(f) + "_" + str(y)+ ".txt"
#        print(fn)
#        os.system('F:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D h1=0 -D h2=0 -D fn="%s" fibers_crop.xml' % (fn));
#        os.rename("fibers_crop.exr", '../results/' + dataset + '/proc_wo_' + str(f) + '_' + str(y) + '_cropped3.exr')
#        