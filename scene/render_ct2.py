'''
Generate the testinfo files needed for volume rendering
Note that this needs python 2.7 and doens't work with 3.5 <<<<
'''

import os
import os.path


    # In[]
def generate_ct2files(dataset, frame0, frame1, yarnNum, skipfactor):
    path = "../output/" + dataset
    ####for generating bounding box in obj
 
    '''
    generate obj mesh
    '''
    for f in range (frame0, frame1+1, skipfactor):
        argin = path + "/volume_" + str(f)  + ".vol"
    
        if not ( os.path.isfile(argin) and os.access(argin, os.R_OK) ) :
            print (argin)
            print (">>>>>>>>> Either file is missing or is not readable <<<<<<<<<<")
        
        argout = path + "/volume_" + str(f) 
        os.system('python createBound.py %s -o %s' %(argin,argout) ) #generate bounding mesh in obj format from AABB volume
        
    
    '''
    generate testinfo
    '''
    for f in range (frame0, frame1+1, skipfactor):
        print(f)
    #    test_info = dataset + "/../test_info_" + str(f) + ".txt"
        test_info = path + "/testinfo_" + str(f) + ".txt"
        with open(test_info, "w") as fout:
            fout.writelines(str(yarnNum) + "\n")
            for y in range (0,yarnNum):
                ### <<< Note that none of the color channels should be 0.0 >>> #### 
                # 150x100
    #            if (y<150):
    #                fout.writelines("0.094 0.176 0.764\n")
    #            else:
    #                fout.writelines("0.931 0.984 0.905\n")
                    
                fout.writelines("0.894 0.894 0.001\n") #knitted
    #            if (y < yarnNum/2):
    #                fout.writelines("0.882 0.762 0.148 \n") #stretch
    #                fout.writelines("0.254 0.623 0.262\n") #push  
    #                fout.writelines("0.927 0.152 0.447\n") #flower100x100
                     
    #            else:
    #                fout.writelines("0.859 0.058 0.191 \n") #stretch
    #                 fout.writelines("0.839 0.764 0.360\n") #push
    #                 fout.writelines("0.99 0.99 0.99\n") #flower100x100
                    
                fn = "F:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/output/" + dataset + "/genYarn_wo_" + str(f) + "_" + str(y) +".txt"
                fout.writelines('%s \n' %fn)
            fout.close()

# In[]
''' render using volume rendering using mitsuba-ct2'''
#
#dataset = 'woven/push/yarn8/100x100'
#dataset = 'woven/stretch/yarn4/100x100'
#dataset = 'woven/arbitrary_pattern/150x100'
#dataset = 'woven/arbitrary_pattern/100x100'
#dataset = 'woven/arbitrary_pattern/512x512'
dataset = 'knitted/slipStitch_stretch'
#dataset = 'knitted/sweater_stretch'
#dataset = 'woven/knitted_push_cube'
f0 = 0
f1 = 0
yarnNum = 1
skipfactor = 2000

generate_ct2files(dataset, f0, f1, yarnNum, skipfactor)

for f in range (f0, f1+1, skipfactor):
    
    testinfo = "../output/"+dataset+"/testinfo_" + str(f) +".txt"
    obj = "../output/"+dataset+"/volume_" + str(f) + ".obj"
    print(testinfo)
    print(obj)
    os.system('F:/sandbox/fiberSimulation/mitsuba_ct2/dist/mitsuba -D testinfo="%s" -D obj="%s" fibers_ct2.xml' % (testinfo, obj) )
    os.rename("fibers_ct2.exr", '../../results/' + dataset + '/50x50_' + str(f) + '_4.exr')
    

#F:/sandbox/fiberSimulation/mitsuba_ct2/dist/mitsuba -D testinfo="../output/knitted/slipStitch_stretch/testinfo_0.txt" -D obj="../output/knitted/slipStitch_stretch/volume_0.obj" fibers_ct2.xml