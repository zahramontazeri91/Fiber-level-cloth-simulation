'''
Generate the testinfo files needed for volume rendering
Note that this needs python 2.7 and doens't work with 3.5 <<<<
And in case of using server to render, use linux-subsystem
'''

import os
import os.path

    # In[]
def generate_ct2files(dataset, frame0, frame1, y0, y1, skipfactor, data):
    path = "../output/" + dataset
    ####for generating bounding box in obj
 
    '''
    generate obj mesh
    '''
#    if (not useServer):
    for f in range (frame0, frame1+1, skipfactor):
        argin = path + "/volume_" + str(f)  + ".vol"
    
        if not ( os.path.isfile(argin) and os.access(argin, os.R_OK) ) :
            print (argin)
            print (">>>>>>>>> Either file is missing or is not readable <<<<<<<<<<")
        
        argout = path + "/volume_" + str(f) 
        os.system('python2.7 createBound.py %s -o %s' %(argin,argout) ) #generate bounding mesh in obj format from AABB volume
            
    '''
    generate testinfo
    '''
    for f in range (frame0, frame1+1, skipfactor):
        test_info = path + "/testinfo_" + str(f) + ".txt"
        with open(test_info, "w") as fout:
            fout.writelines(str(y1) + "\n")
            for y in range (y0, y1):
                ### <<< Note that none of the color channels should be 0.0 >>> #### 
                # 150x100
    #            if (y<150):
    #                fout.writelines("0.094 0.176 0.764\n")
    #            else:
    #                fout.writelines("0.931 0.984 0.905\n")
                    
                fout.writelines("0.890625 0.726562 0.3828\n") #knitted
    #            if (y < yarnNum/2):
    #                fout.writelines("0.882 0.762 0.148 \n") #stretch
    #                fout.writelines("0.254 0.623 0.262\n") #push  
    #                fout.writelines("0.927 0.152 0.447\n") #flower100x100
                     
    #            else:
    #                fout.writelines("0.859 0.058 0.191 \n") #stretch
    #                 fout.writelines("0.839 0.764 0.360\n") #push
    #                 fout.writelines("0.99 0.99 0.99\n") #flower100x100
    
                if (useServer):
                    #write same for all frames because when using server, I'll only copy for that particular frame. So we can use cached.dat 
                    fout.writelines('../render/output/genYarn_%d.txt\n' %( y) )
                else:
                    fn = "../output/" + dataset + "/genYarn_" + cat +"_" + str(f) + "_" + str(y) +".txt" 
                    fout.writelines('%s \n' %fn)
            fout.close()
        print('successfully generated testinfo!')

# In[]
''' render using volume rendering using mitsuba-ct2'''
useServer = 1
server = 'ampere00'
cat = 'NN'
isFirst = 1

if (useServer):
    os.chdir('/mnt/f/sandbox/fiberSimulation/yarn_generation_project/YarnGeneration/scene')
    os.system('scp fibers_ct2.xml hquser@'+server+'.cs.columbia.edu:~/external/libelasto/render/scene/fibers_ct2.xml')


#info = 'yarn8/slipstitchrib_push 24000 24000 0 1 -k 500 -t 0 -v 85832 -s 1 -vol 1'
info = 'yarn8/slipstitchrib_stretch 0 0 0 1 -k 1000 -t 0 -v 85832 -s 1 -vol 1'

split = info.split()
dataset = split[0]
f0 = int(split[1])
f1 = int(split[2])
y0 = int(split[3])
y1 = int(split[4])
skipfactor = int(split[6])

generate_ct2files(dataset, f0, f1, y0, y1, skipfactor, cat)

for f in range (f0, f1+1, skipfactor):
            
    testinfo = "../output/"+dataset+"/testinfo_" + str(f) +".txt"
    obj = "../output/"+dataset+"/volume_" + str(f) + ".obj"
    print(testinfo)
    print(obj)      
        
    if (useServer):    
        for y in range (y0, y1):
            fn = "../output/" + dataset + "/genYarn_" + cat +"_" + str(f) + "_" + str(y) +".txt"              
            os.system('scp %s hquser@%s.cs.columbia.edu:~/external/libelasto/render/output/genYarn_%d.txt\n' %(fn, server, y) ) 
            
        os.system('scp %s hquser@%s.cs.columbia.edu:~/external/libelasto/render/output/volume.obj ' % (obj, server))
#        if (isFirst):
        if (1):
            os.system('scp %s hquser@%s.cs.columbia.edu:~/external/libelasto/render/output/testinfo.txt ' % (testinfo, server))
            command = "-r 100 -D f=0 -D testinfo='../render/output/testinfo.txt' -D obj='../render/output/volume.obj' ../render/scene/fibers_ct2.xml"
            isFirst = 0
        else:
            command = "-r 100 -D f=1 -D testinfo='cached.dat' -D obj='../render/output/volume.obj' ../render/scene/fibers_ct2.xml"
        os.system('ssh hquser@'+server+'.cs.columbia.edu '\
                  + "'cd external/libelasto/mitsuba-ct2; . setpath.sh; dist/mitsuba %s'" %command  )
        os.system('scp hquser@%s.cs.columbia.edu:~/external/libelasto/render/scene/fibers_ct2.exr /mnt/f/sandbox/fiberSimulation/yarn_generation_project/results/%s/test/%s_color2_%d.exr'  %(server, dataset, cat, f))
    else:
        os.system('F:/sandbox/fiberSimulation/mitsuba_ct2/dist/mitsuba -r 60 -D f=0 -D testinfo="%s" -D obj="%s" fibers_ct2.xml' % (testinfo, obj) )
        os.rename("fibers_ct2.exr", '../../results/' + dataset + '/'+ cat +'_' + str(f) + '.exr')
    







 




























