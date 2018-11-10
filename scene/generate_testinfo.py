
import os
import os.path
#
#dataset = 'woven/push/yarn8/100x100'
#dataset = 'woven/stretch/yarn4/100x100'
#dataset = 'woven/arbitrary_pattern/150x100'
#dataset = 'woven/arbitrary_pattern/100x100'
#dataset = 'woven/arbitrary_pattern/512x512'
dataset = 'woven/knitted'
#dataset = 'woven/knitted_push_cube'


path = "../YarnGeneration/output/" + dataset
frame0 = 0
frame1 = 40000
yarn = 1
skipfactor = 2000

####for generating bounding box in obj

for f in range (frame0, frame1+1, skipfactor):
    argin = path + "/volume_" + str(f)  + ".vol"

    if not ( os.path.isfile(argin) and os.access(argin, os.R_OK) ) :
        print ">>>>>>>>> Either file is missing or is not readable <<<<<<<<<<"
    
    argout = path + "/volume_" + str(f) 
    os.system('python createBound.py %s -o %s' %(argin,argout) )
    

# In[]
    
for f in range (frame0, frame1+1, skipfactor):
    print(f)
#    test_info = dataset + "/../test_info_" + str(f) + ".txt"
    test_info = path + "/testinfo_" + str(f) + ".txt"
    with open(test_info, "w") as fout:
        fout.writelines(str(yarn) + "\n")
        for y in range (0,yarn):
            
            # 150x100
#            if (y<150):
#                fout.writelines("0.094 0.176 0.764\n")
#            else:
#                fout.writelines("0.931 0.984 0.905\n")
                
            fout.writelines("0.431 0.584 0.905\n") #knitted
#            if (y < yarn/2):
#                fout.writelines("0.882 0.762 0.148 \n") #stretch
#                fout.writelines("0.254 0.623 0.262\n") #push  
#                fout.writelines("0.927 0.152 0.447\n") #flower100x100
                 
#            else:
#                fout.writelines("0.859 0.058 0.191 \n") #stretch
#                 fout.writelines("0.839 0.764 0.360\n") #push
#                 fout.writelines("0.99 0.99 0.99\n") #flower100x100
                
#            fn = "D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/output/" + dataset + "/genYarn_NN_fly_wo_" + str(f) + "_" + str(y) +".txt"
#            fn = "D:/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/output/" + dataset + "/centerlines/curve_" + str(f) + "_" + str(y) +"_ds_mitsuba.txt"
                
            fn = "../../output/genYarn_NN_fly_" + str(f) + "_" + str(y) +".txt"
#            fn = "../../twist/genYarn_NN_wo_" + str(f) + "_" + str(y) +".txt"
#            fn = "../../centerlines/curve_" + str(f) + "_" + str(y) +"_ds_mitsuba.txt" 
                
#            fn = "centerlines/curve_" + str(f) + "_" + str(y) +"_ds_mitsuba.txt" 
            fout.writelines('%s \n' %fn)
        fout.close()