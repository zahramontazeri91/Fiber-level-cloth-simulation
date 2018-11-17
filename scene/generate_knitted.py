# -*- coding: utf-8 -*-
"""
Generate mitsuba scene file for knitted example
Note that this uses python 2.7 and doens't work with 3.5 because 
createBoun.py works with 2.7

This uses fiber-mitsuba to render knitted sample
@author: zahra
"""

def knittedXMLgenerator(dataset, frame, yarnNum, outname):
    print(outname)
    with open(outname, "w") as fout:
        
    #    fout.writelines('<?xml version='1.0' encoding='utf-8'?>')
        fout.writelines('<scene version="0.4.4"> \n')
    
        fout.writelines('\t<integrator type="volpath_simple"> \n')
        fout.writelines('\t\t<integer name="maxDepth" value="25"/> \n')
        fout.writelines('\t\t<integer name="rrDepth" value="30"/> \n')
        fout.writelines('\t</integrator> \n')
        fout.writelines('\n') 
        
        for i in range (0,yarnNum):
            ### shape hair
            fout.writelines('\t<shape type="hair">\n')
#            fn = "F:/YarnGeneration/output/" + dataset + "/genYarn_wo_" + str(frame) + "_" + str(i)+ ".txt"               
            
            fout.writelines('\t\t<string name="filename" value="$fn%d"/>\n' % (i) )
            fout.writelines('\t\t<float name="radius" value="0.002"/> 	\n')	
            fout.writelines('\t\t<float name="angleThreshold" value="0.001"/>\n')
            
            # use bsdf
            fout.writelines('\t\t<bsdf type="roughplastic">\n')
            fout.writelines('\t\t\t<spectrum name="diffuseReflectance" value="0.856875, 0.796875, 0.0"/>\n')
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
#            fout.writelines('\t\t\t\t<spectrum name="colorTT" value="0.859,0.058,0.191"/>\n')
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
        fout.writelines('\t\t<spectrum name="radiance" value="0.15"/> \n')
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
        fout.writelines('\t\t\t<lookAt origin="1 0 -10" target="1 0 -1.5" up="0 1 0"/>   \n')  #for sweater
            
        fout.writelines('\t\t</transform>\n')
        fout.writelines('\t\t<float name="fov" value="50"/>\n') #3.8 for simul xml
        fout.writelines('\n')
        
        ### sampler
        fout.writelines('\t\t<sampler type="ldsampler">\n')
#        fout.writelines('\t\t\t<integer name="sampleCount" value="1024"/>\n')
        fout.writelines('\t\t\t<integer name="sampleCount" value="32"/>\n')
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
dataset = 'knitted/slipStitch_stretch'
outname = 'slipStitch_stretch.xml'
f0 = 0
f1 = 0
yarnNum = 1

# for generating xml file
for i in range (f0, f1+1, 500):
    print (i)
    frame = i
    knittedXMLgenerator(dataset, frame, yarnNum, outname)
