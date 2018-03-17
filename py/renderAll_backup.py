import os 

os.chdir('D:/sandbox/fiberSimulation/yarn_generation_project/scene')
#dataset = "pattern/yarn4/spacing1.0x/00011"
dataset = "twist/yarn4/damp"

### procedural renderings
#skipFactor = 500
#yarnNum = 1
#for i in range (49000/skipFactor,49500/skipFactor+1):
#    f = i*skipFactor
#    for y in range (0,yarnNum):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
#        os.system('D:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_' + str(f) + '_' + str(y) + '_dg2.exr')


###############
dataset = "woven/yarn4/spacing1.0x/00011"
##### procedural renderings
skipFactor = 500
yarnNum = 46

for i in range (5000/skipFactor,8000/skipFactor+1):
    f = i*skipFactor
    fn_a = []
    for y in range(0,yarnNum):
        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
        fn_a.append(fn)
        

    os.system('D:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba \
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

#    os.system('D:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn0="%s" fibers_grid.xml'% (fn_a[0] ))        
    os.rename("fibers_grid.exr", '../results/' + dataset + '/proc_NN_' + str(f) + '.exr')
        
        
        
###################

dataset = "woven/yarn4/spacing1.0x/00011"

### procedural renderings
#skipFactor = 500
#yarnNum = 46
#for i in range (0/skipFactor,2500/skipFactor+1):
#    f = i*skipFactor
#    for y in range (23,24):
#        fn = "../YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(f) + "_" + str(y)+ ".txt"
#        os.system('D:/sandbox/fiberSimulation/dist_fiber_mitsuba/dist/mitsuba -D fn="%s" fibers_grid.xml' % (fn));
#        os.rename("fibers.exr", '../results/' + dataset + '/proc_' + str(f) + '_' + str(y) + '_damp_new.exr')
#








