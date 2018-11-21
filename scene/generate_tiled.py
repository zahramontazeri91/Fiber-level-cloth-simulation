"""
Created on Thu Oct 18 10:49:14 2018
generate xml file for mitsuba for tiling the single-yarn
@author: zahra
"""

def genScene (xmlfile, spp, yarnType=""):
    print('writing xml file ...\n' + xmlfile)
    isYarn4=0
    isYarn8=0
    isYarn11=0
    isYarn100=0 
    if (yarnType=='yarn4'):
        isYarn4=1
    elif (yarnType=='yarn8'):
        isYarn8=1
    elif (yarnType=='yarn11'):
        isYarn11=1
    elif (yarnType=='yarn100'):
        isYarn100=1

    with open(xmlfile, "w") as fout:

        sz = 33 #cylinder (for stretch case)
#        sz = 22 #flat tiling (for teeth cases)

        fout.writelines('<scene version="0.4.4"> \n')    
        fout.writelines('\t<integrator type="path"> \n')
        fout.writelines('\t\t<integer name="maxDepth" value="25"/> \n')
        fout.writelines('\t\t<integer name="rrDepth" value="30"/> \n')
        fout.writelines('\t</integrator> \n')
        fout.writelines('\n') 

        if (isYarn8):
            #### cylinders 1.2 00110
            fout.writelines('\t<shape type="cylinder">\n')
            fout.writelines('\t\t<float name="radius" value="0.16"/>\n')
            fout.writelines('\t\t<point name="p0" x="-1.8" y="$h2" z="0"/>\n')
            fout.writelines('\t\t<point name="p1" x="1.8" y="$h2" z="0"/>\n')
            fout.writelines('\t\t<bsdf type="diffuse">\n')
            fout.writelines('\t\t\t<spectrum name="reflectance" value="0.3"/>\n')
            fout.writelines('\t\t</bsdf>	\n')
            fout.writelines('\t</shape>\n')
            fout.writelines('\t<shape type="cylinder">\n')
            fout.writelines('\t\t<float name="radius" value="0.16"/>\n')
            fout.writelines('\t\t<point name="p0" x="-1.8" y="$h2" z="0.48"/>\n')
            fout.writelines('\t\t<point name="p1" x="1.8" y="$h2" z="0.48"/>\n')
            fout.writelines('\t\t<bsdf type="diffuse">\n')
            fout.writelines('\t\t\t<spectrum name="reflectance" value="0.3"/>\n')
            fout.writelines('\t\t</bsdf>	\n')
            fout.writelines('\t</shape>\n')
            fout.writelines('\t<shape type="cylinder">\n')
            fout.writelines('\t\t<float name="radius" value="0.16"/>\n')
            fout.writelines('\t\t<point name="p0" x="-1.8" y="$h2" z="-0.48"/>\n')
            fout.writelines('\t\t<point name="p1" x="1.8" y="$h2" z="-0.48"/>\n')
            fout.writelines('\t\t<bsdf type="diffuse">\n')
            fout.writelines('\t\t\t<spectrum name="reflectance" value="0.3"/>\n')
            fout.writelines('\t\t</bsdf>	\n')
            fout.writelines('\t</shape>\n')
            fout.writelines('\t<shape type="cylinder">\n')
            fout.writelines('\t\t<float name="radius" value="0.16"/>\n')
            fout.writelines('\t\t<point name="p0" x="-1.8" y="$h1" z="0.96"/>\n')
            fout.writelines('\t\t<point name="p1" x="1.8" y="$h1" z="0.96"/>\n')
            fout.writelines('\t\t<bsdf type="diffuse">\n')
            fout.writelines('\t\t\t<spectrum name="reflectance" value="0.3"/>\n')
            fout.writelines('\t\t</bsdf>	\n')
            fout.writelines('\t</shape>\n')
            fout.writelines('\t<shape type="cylinder">\n')
            fout.writelines('\t\t<float name="radius" value="0.16"/>\n')
            fout.writelines('\t\t<point name="p0" x="-1.8" y="$h1" z="-0.96"/>\n')
            fout.writelines('\t\t<point name="p1" x="1.8" y="$h1" z="-0.96"/>\n')
            fout.writelines('\t\t<bsdf type="diffuse">\n')
            fout.writelines('\t\t\t<spectrum name="reflectance" value="0.3"/>\n')
            fout.writelines('\t\t</bsdf>	\n')
            fout.writelines('\t</shape>\n')
            
        if (isYarn11):  
            ######## cylinder 1.6
            fout.writelines('\t<shape type="cylinder">\n')
            fout.writelines('\t\t<float name="radius" value="0.16"/>\n')
            fout.writelines('\t\t<point name="p0" x="-1.8" y="$h1" z="0"/>\n')
            fout.writelines('\t\t<point name="p1" x="1.8" y="$h1" z="0"/>\n')
            fout.writelines('\t\t<bsdf type="diffuse">\n')
            fout.writelines('\t\t\t<spectrum name="reflectance" value="0.3"/>\n')
            fout.writelines('\t\t</bsdf>	\n')
            fout.writelines('\t</shape>\n')
        
            fout.writelines('\t<shape type="cylinder">\n')
            fout.writelines('\t\t<float name="radius" value="0.16"/>\n')
            fout.writelines('\t\t<point name="p0" x="-1.8" y="$h2" z="0.64"/>\n')
            fout.writelines('\t\t<point name="p1" x="1.8" y="$h2" z="0.64"/>\n')
            fout.writelines('\t\t<bsdf type="diffuse">\n')
            fout.writelines('\t\t\t<spectrum name="reflectance" value="0.3"/>\n')
            fout.writelines('\t\t</bsdf>	\n')
            fout.writelines('\t</shape>\n')
        
            fout.writelines('\t<shape type="cylinder">\n')
            fout.writelines('\t\t<float name="radius" value="0.16"/>\n')
            fout.writelines('\t\t<point name="p0" x="-1.8" y="$h2" z="-0.64"/>\n')
            fout.writelines('\t\t<point name="p1" x="1.8" y="$h2" z="-0.64"/>\n')
            fout.writelines('\t\t<bsdf type="diffuse">\n')
            fout.writelines('\t\t\t<spectrum name="reflectance" value="0.3"/>\n')
            fout.writelines('\t\t</bsdf>	\n')
            fout.writelines('\t</shape>\n')

        ### shape hair
        fout.writelines('\t<shape type="shapegroup" id="myShapeGroup">\n')
        fout.writelines('\t\t<shape type="hair">\n')
        fout.writelines('\t\t\t<string name="filename" value="$fn"/>\n'  )
        fout.writelines('\t\t\t<float name="radius" value="0.002"/> 	\n')	
        fout.writelines('\t\t\t<float name="angleThreshold" value="0.001"/>\n')
        
        if (isYarn4):
            #yarn4
            fout.writelines('\t\t<subsurface type="fibershader">\n')
            fout.writelines('\t\t\t<boolean name="useRandomInteractionPoint" value="true"/>\n')
            fout.writelines('\t\t\t<boolean name="sampleInteractionPointFromCircumference" value="false"/>\n')
            fout.writelines('\t\t\t<fiberscat type="simpfabric5">\n')
            fout.writelines('\t\t\t\t<float name="kD" value="0"/>\n')
            fout.writelines('\t\t\t\t<spectrum name="colorD" value="0.99,0.99,0.99"/>\n')
            fout.writelines('\t\t\t\t<spectrum name="colorR" value="0.25,0.25,0.3"/>\n')
            fout.writelines('\t\t\t\t<float name="betaR" value="1.23799781194"/>\n')
            fout.writelines('\t\t\t\t<float name="betaTT" value="9.99999999999"/>\n')
            fout.writelines('\t\t\t\t<float name="gammaTT" value="25.9890848941"/>\n')
            fout.writelines('\t\t\t\t<spectrum name="colorTT" value="0.45,0.86,0.95"/>\n')
            fout.writelines('\t\t\t\t<float name="alpha" value="5"/>\n')
            fout.writelines('\t\t\t</fiberscat>\n')
            fout.writelines('\t\t</subsurface> \n')
        elif (isYarn8):
            #yarn8
            fout.writelines('\t\t<subsurface type="fibershader">\n')
            fout.writelines('\t\t\t<boolean name="useRandomInteractionPoint" value="true"/>\n')
            fout.writelines('\t\t\t<boolean name="sampleInteractionPointFromCircumference" value="false"/>\n')
            fout.writelines('\t\t\t<fiberscat type="simpfabric5">\n')
            fout.writelines('\t\t\t\t<float name="kD" value="0"/>\n')
            fout.writelines('\t\t\t\t<spectrum name="colorD" value="0.99,0.99,0.99"/>\n')
            fout.writelines('\t\t\t\t<spectrum name="colorR" value="0.1,0.1,0.05"/>\n')
            fout.writelines('\t\t\t\t<float name="betaR" value="0.2"/>\n')
            fout.writelines('\t\t\t\t<float name="betaTT" value="27"/>\n')
            fout.writelines('\t\t\t\t<float name="gammaTT" value="38"/>\n')
            fout.writelines('\t\t\t\t<spectrum name="colorTT" value="0.93,0.53,0.01"/>\n')
            fout.writelines('\t\t\t\t<float name="alpha" value="5"/>\n')
            fout.writelines('\t\t\t</fiberscat>\n')
            fout.writelines('\t\t</subsurface> \n')
        elif (isYarn11):
            #yarn11
            fout.writelines('\t\t<subsurface type="fibershader">\n')
            fout.writelines('\t\t\t<boolean name="useRandomInteractionPoint" value="true"/>\n')
            fout.writelines('\t\t\t<boolean name="sampleInteractionPointFromCircumference" value="false"/>\n')
            fout.writelines('\t\t\t<fiberscat type="simpfabric5">\n')
            fout.writelines('\t\t\t\t<float name="kD" value="0"/>\n')
            fout.writelines('\t\t\t\t<spectrum name="colorD" value="0.99,0.99,0.99"/>\n')
            fout.writelines('\t\t\t\t<spectrum name="colorR" value="0.1,0.1,0.05"/>\n')
            fout.writelines('\t\t\t\t<float name="betaR" value="4.0"/>\n')
            fout.writelines('\t\t\t\t<float name="betaTT" value="10"/>\n')
            fout.writelines('\t\t\t\t<float name="gammaTT" value="20"/>\n')
            fout.writelines('\t\t\t\t<spectrum name="colorTT" value="0.88,0.83,0.01"/>\n')
            fout.writelines('\t\t\t\t<float name="alpha" value="5"/>\n')
            fout.writelines('\t\t\t</fiberscat>\n')
            fout.writelines('\t\t</subsurface> \n')
        elif (isYarn100):
            #yarn100
            fout.writelines('\t\t<subsurface type="fibershader">\n')
            fout.writelines('\t\t\t<boolean name="useRandomInteractionPoint" value="true"/>\n')
            fout.writelines('\t\t\t<boolean name="sampleInteractionPointFromCircumference" value="false"/>\n')
            fout.writelines('\t\t\t<fiberscat type="simpfabric5">\n')
            fout.writelines('\t\t\t\t<float name="kD" value="0"/>\n')
            fout.writelines('\t\t\t\t<spectrum name="colorD" value="0.99,0.99,0.99"/>\n')
            fout.writelines('\t\t\t\t<spectrum name="colorR" value="0.1,0.1,0.1"/>\n')
            fout.writelines('\t\t\t\t<float name="betaR" value="1.23799781194"/>\n')
            fout.writelines('\t\t\t\t<float name="betaTT" value="9.99999999999"/>\n')
            fout.writelines('\t\t\t\t<float name="gammaTT" value="25.9890848941"/>\n')
            fout.writelines('\t\t\t\t<spectrum name="colorTT" value="0.9,0.24,0.48"/>\n')
            fout.writelines('\t\t\t\t<float name="alpha" value="5"/>\n')
            fout.writelines('\t\t\t</fiberscat>\n')
            fout.writelines('\t\t</subsurface> \n')   
        else:
            # use bsdf
            fout.writelines('\t\t<bsdf type="roughplastic">\n')
            fout.writelines('\t\t\t<spectrum name="diffuseReflectance" value="0.8, 0.1, 0.1"/>\n')
            fout.writelines('\t\t\t<float name="alpha" value="0.05"/>\n')
            fout.writelines('\t\t</bsdf> \n')
    
    
        
        fout.writelines('\t\t</shape>\n')
        fout.writelines('\t</shape>\n')
        fout.writelines('\n')
        
        #use cylinder tiling        
        for i in range (0,sz):
            ### shape hair
            fout.writelines('\t<shape type="instance">\n')
            fout.writelines('\t\t<ref id="myShapeGroup"/>\n' )     
            fout.writelines('\t\t<transform name="toWorld">\n')	
            
#            x_t = 0.11 * float(i)  #cylinder tiling for synthesized result
                        
            x_t = 0.06 * float(i)
            y_t = 1.0 - x_t * x_t/1.8
            fout.writelines('\t\t\t<translate x="%.6f" y="%.6f" />\n' % (x_t, y_t) )#cylinder tiling
#            fout.writelines('\t\t\t<translate x="%.6f" />\n' % (1.2*x_t))#flat tiling
            fout.writelines('\t\t</transform> \n')
            #use actualbrdf
            fout.writelines('\t</shape>\n')
            fout.writelines('\n')
    
        for i in range (1,sz):
            ### shape hair
            fout.writelines('\t<shape type="instance">\n')
            fout.writelines('\t\t<ref id="myShapeGroup"/>\n' )     
            fout.writelines('\t\t<transform name="toWorld">\n')	
            
#            x_t = -1.0 * 0.11 * float(i) #cylinder tiling for synthesized result
            
            x_t = -1.0 * 0.06 * float(i)
            y_t = 1.0 - x_t * x_t/1.8
            fout.writelines('\t\t\t<translate x="%.6f" y="%.6f" />\n' % (x_t, y_t) )#cylinder tiling
#            fout.writelines('\t\t\t<translate x="%.6f" />\n' % (1.2*x_t))#flat tiling
            fout.writelines('\t\t</transform> \n')
            fout.writelines('\t</shape>\n')
            fout.writelines('\n')
        
        ####  rectangle
        fout.writelines('\t<shape type="rectangle">\n')
        fout.writelines('\t\t<transform name="toWorld">\n')
        fout.writelines('\t\t\t<rotate x="1" angle="-90"/>\n')
        fout.writelines('\t\t\t<translate y="-0.1"/>\n')
        fout.writelines('\t\t\t<scale x="6" y="6" z="6"/>\n')
        fout.writelines('\t\t</transform>\n')
        fout.writelines('\t\t<bsdf type="diffuse">\n')
        fout.writelines('\t\t\t<texture type="bitmap" name="reflectance">\n')
        fout.writelines('\t\t\t\t<string name="filename" value="concrete.jpg"/>\n')
        fout.writelines('\t\t\t</texture>\n')
        fout.writelines('\t\t</bsdf>\n')
        fout.writelines('\t</shape>\n')
        
        ####  light
        fout.writelines('\t<emitter type="constant">\n')
        fout.writelines('\t\t<spectrum name="radiance" value="0.3"/> \n')
        fout.writelines('\t</emitter>\n')
        
        ####  light
    #    fout.writelines('\t<shape type="sphere">\n')
    #    fout.writelines('\t\t<point name="center" x="-5.0" y="80.0" z="0.0"/>\n')
    #    fout.writelines('\t\t<float name="radius" value="10.0"/>\n')
    #    fout.writelines('\t\t<emitter type="area">\n')
    #    fout.writelines('\t\t\t<spectrum name="radiance" value="40"/> \n')
    #    fout.writelines('\t\t</emitter>\n')
    #    fout.writelines('\t</shape>\n')
    #    fout.writelines('\n')
        
        ####  light
        fout.writelines('\t<shape type="rectangle"> \n')
        fout.writelines('\t\t<transform name="toWorld"> \n')
        fout.writelines('\t\t\t<rotate x="1" angle="90"/> \n')
        fout.writelines('\t\t\t<scale x="1" y="40" z="20"/> \n')
        fout.writelines('\t\t\t<translate x="0.0" y="100.0" z="0.0"/> \n')
        fout.writelines('\t\t</transform> \n')
        fout.writelines('\t\t<emitter type="area"> \n')
        fout.writelines('\t\t\t<spectrum name="radiance" value="50"/>  \n')
        fout.writelines('\t\t</emitter> \n')
        fout.writelines('\t</shape> \n')
        
        
    #    ### sensor
        fout.writelines('\t<sensor type="perspective">\n')
        fout.writelines('\t\t<string name="fovAxis" value="smaller"/>\n')
        
        fout.writelines('\t\t<transform name="toWorld">\n')
    #    fout.writelines('\t\t\t<lookAt origin="-17 40 15" target="0 0.9 1.1" up="0 1 0"/>  \n') #9:16 ratio for video
    #    fout.writelines('\t\t\t<lookAt origin="-17 40 15" target="0 -0.3 -0.0" up="0 1 0"/>  \n') #cylinder tiling for snthesized result
        fout.writelines('\t\t\t<lookAt origin="-17 40 15" target="0 -0.9 -0.4" up="0 1 0"/>  \n') #cylinder tiling
#        fout.writelines('\t\t\t<lookAt origin="-17 30 5" target="0 0 0" up="0 1 0"/>  \n') #flat tiling for video
        fout.writelines('\t\t</transform>\n')
        
    #    fout.writelines('\t\t<float name="fov" value="7.1"/>\n')#cylinder tiling for snthesized result
        fout.writelines('\t\t<float name="fov" value="6.3"/>\n')
        fout.writelines('\t\t<sampler type="ldsampler">\n')
        fout.writelines('\t\t\t<integer name="sampleCount" value="%d"/>\n' %spp) #256
        fout.writelines('\t\t</sampler>\n')
        fout.writelines('\t\t<film id="film" type="hdrfilm">\n')
        fout.writelines('\t\t\t<integer name="width" value="1024"/>\n') # for video
        fout.writelines('\t\t\t<integer name="height" value="576"/>\n')
        fout.writelines('\t\t\t<rfilter type="gaussian"/>\n')
        fout.writelines('\t\t\t<boolean name="banner" value="false" />\n')
        fout.writelines('\t\t</film>\n')
        fout.writelines('\t</sensor>\n')
        fout.writelines('\n')
    
        fout.writelines('</scene>')
    
    fout.close()



