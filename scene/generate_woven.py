# -*- coding: utf-8 -*-
"""
Created on Mon Apr 30 17:07:09 2018

@author: zahra
"""

isSimul = 0

def wovenXMLgenerator(dataset, frame, sz):
    outname = 'woven/' + str(sz) +'x'+str(sz) + '_' + str(frame) + '_wo.xml'
    print(outname)
    with open(outname, "w") as fout:
        
    #    fout.writelines('<?xml version='1.0' encoding='utf-8'?>')
        fout.writelines('<scene version="0.4.4"> \n')
    
        fout.writelines('\t<integrator type="volpath_simple"> \n')
        fout.writelines('\t\t<integer name="maxDepth" value="25"/> \n')
        fout.writelines('\t\t<integer name="rrDepth" value="30"/> \n')
        fout.writelines('\t</integrator> \n')
        fout.writelines('\n') 
        
        for i in range (0,sz+sz):
            if (i==0 or i==5 or i==6 or i==11): #remove the boundry yarns
                continue
            ### shape hair
            fout.writelines('\t<shape type="hair">\n')
            if (isSimul):
                fn = "F:/YarnGeneration/fibersim/" + dataset + "/simul_frame_" + str(frame) + "_" + str(i)+ ".txt"
            else:
#                fn = "F:/YarnGeneration/output/" + dataset + "/genYarn_NN_" + str(frame) + "_" + str(i)+ "_us.txt"
                fn = "F:/YarnGeneration/output/" + dataset + "/genYarn_wo_" + str(frame) + "_" + str(i)+ ".txt"               
            
            fout.writelines('\t\t<string name="filename" value="%s"/>\n' % (fn) )
            fout.writelines('\t\t<float name="radius" value="0.002"/> 	\n')	
            fout.writelines('\t\t<float name="angleThreshold" value="0.001"/>\n')
            
            # use bsdf
            fout.writelines('\t\t<bsdf type="roughplastic">\n')
            if (i<sz):
                fout.writelines('\t\t\t<spectrum name="diffuseReflectance" value="0.796875, 0.796875, 0.0"/>\n')
            else:
                fout.writelines('\t\t\t<spectrum name="diffuseReflectance" value="0.1, 0.1, 0.8997"/>\n')
            fout.writelines('\t\t\t<float name="alpha" value="0.15"/>\n')
            fout.writelines('\t\t</bsdf> \n')


            #use actualbrdf
            #yarn8
#            fout.writelines('\t\t<subsurface type="fibershader">\n')
#            fout.writelines('\t\t\t<boolean name="useRandomInteractionPoint" value="true"/>\n')
#            fout.writelines('\t\t\t<boolean name="sampleInteractionPointFromCircumference" value="false"/>\n')
#            fout.writelines('\t\t\t<fiberscat type="simpfabric5">\n')
#            fout.writelines('\t\t\t\t<float name="kD" value="0"/>\n')
#            fout.writelines('\t\t\t\t<spectrum name="colorD" value="0.99,0.99,0.99"/>\n')
#            fout.writelines('\t\t\t\t<spectrum name="colorR" value="0.1,0.1,0.05"/>\n')
#            if (i<sz):
#                fout.writelines('\t\t\t\t<spectrum name="colorTT" value="0.882,0.762,0.148"/>\n')
#            else:
#                fout.writelines('\t\t\t\t<spectrum name="colorTT" value="0.859,0.058,0.191"/>\n')
#            fout.writelines('\t\t\t\t<float name="betaR" value="0.2"/>\n')
#            fout.writelines('\t\t\t\t<float name="betaTT" value="27"/>\n')
#            fout.writelines('\t\t\t\t<float name="gammaTT" value="38"/>\n')
#            fout.writelines('\t\t\t\t<float name="alpha" value="5"/>\n')
#            fout.writelines('\t\t\t</fiberscat>\n')
#            fout.writelines('\t\t</subsurface> \n')            
            
            
            fout.writelines('\t</shape>\n')
            fout.writelines('\n')
    
        
        ####  light
        fout.writelines('\t<emitter type="constant">\n')
        fout.writelines('\t\t<spectrum name="radiance" value="0.0"/> \n')
        fout.writelines('\t</emitter>\n')

        fout.writelines('\t<shape type="rectangle"> \n')
        fout.writelines('\t\t<transform name="toWorld"> \n')
        fout.writelines('\t\t\t<rotate x="1" angle="90"/> \n')
        fout.writelines('\t\t\t<scale x="30" y="30" z="30"/> \n')
        fout.writelines('\t\t\t<translate x="0.0" y="200.0" z="0.0"/> \n')
        fout.writelines('\t\t</transform> \n')
        fout.writelines('\t\t<emitter type="area"> \n')
        fout.writelines('\t\t\t<spectrum name="radiance" value="80"/>  \n')
        fout.writelines('\t\t</emitter> \n')
        fout.writelines('\t</shape> \n')
        
#        fout.writelines('\t<shape type="sphere">\n')
#        fout.writelines('\t\t<point name="center" x="-100.0" y="50.0" z ="0.0"/>\n')
#        fout.writelines('\t\t<float name="radius" value="10.0"/>\n')
#        fout.writelines('\t\t<emitter type="area">\n')
#        fout.writelines('\t\t\t<spectrum name="radiance" value="200"/> \n')
#        fout.writelines('\t\t</emitter>\n')
#        fout.writelines('\t</shape>\n')
#        fout.writelines('\n')
        
        ### sensor
        fout.writelines('\t<sensor type="perspective">\n')
        fout.writelines('\t\t<string name="fovAxis" value="smaller"/>\n')
        fout.writelines('\t\t<transform name="toWorld">\n')
        if (isSimul):
            fout.writelines('\t\t\t<lookAt origin="0.25 10 0.0" target="0.25 0 0.0" up="1 0 0"/>  \n')  #for simul xml
        else:
            fout.writelines('\t\t\t<lookAt origin="0.24 10 0.0" target="0.24 0 0.0" up="1 0 0"/>  \n')  #for NN xml
            
        fout.writelines('\t\t</transform>\n')
        
        if(isSimul):
            fout.writelines('\t\t<float name="fov" value="3.0"/>\n') #for NN xml
        else:
            fout.writelines('\t\t<float name="fov" value="3.8"/>\n') #for simul xml
        fout.writelines('\n')
        
        ### sampler
        fout.writelines('\t\t<sampler type="ldsampler">\n')
#        fout.writelines('\t\t\t<integer name="sampleCount" value="1024"/>\n')
        fout.writelines('\t\t\t<integer name="sampleCount" value="16"/>\n')
        fout.writelines('\t\t</sampler>\n')
        fout.writelines('\t\t<film id="film" type="hdrfilm">\n')
        fout.writelines('\t\t\t<integer name="width" value="512"/>\n')
        fout.writelines('\t\t\t<integer name="height" value="512"/>\n')
        fout.writelines('\t\t\t<rfilter type="gaussian"/>\n')
        fout.writelines('\t\t\t<boolean name="banner" value="false" />\n')
        fout.writelines('\t\t</film>\n')
        fout.writelines('\t</sensor>\n')
        fout.writelines('\n')
    
        fout.writelines('</scene>')
    
    fout.close()


# In[]:
#dataset = 'woven/push/yarn4/100x100'
#dataset = 'woven/stretch/yarn4/100x100'
#dataset = 'woven/6x6'
dataset = 'knitted/sweater_stretch'
    


# for generating xml file
sz = 17
if (isSimul):
#    f = 10000
#    f = 1000000 #for frame-2
#    f = 6000
    f = 1000
else:
#    f = 6000
    f = 0
#    f = 15000

for i in range (f, f+1, 500):
    print (i)
    frame = i
    wovenXMLgenerator(dataset, frame, sz)
