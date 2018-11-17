import numpy as np
import copy

def isInside(p, lb, ub):
    for i in range(0, 3):
        if p[i] < lb[i] or p[i] > ub[i]:
            return False
    return True

def clip(p, q, lb, ub):
    assert not isInside(p, lb, ub) and isInside(q, lb, ub)
    L = p
    R = q
    while np.linalg.norm(L - R) > 1e-8:
        mid = 0.5*(L + R)
        if isInside(mid, lb, ub):
            R = mid
        else:
            L = mid
    return R

center = np.array([0.24, 0.0, 0.0])
d = 10
ntile = 3

#path_genYarn = "../../YarnGeneration/data/woven/6x6/simul_frame_6000_"
#basename_in = path_genYarn + "%d.txt"
#basename_out = path_genYarn + "%d_tiled.txt"
#basename_xml = "genYarn_NN_wo_6000.xml"
## fov = 3.8
#fov = 3.4
#nsmooth = 25


#path_genYarn = "../../YarnGeneration/output/woven/6x6/genYarn_NN_wo_6000_"
#basename_in = path_genYarn + "%d.txt"
#basename_out = path_genYarn + "%d_tiled.txt"
#basename_xml = "genYarn_NN_wo_6000.xml"
## fov = 3.8
#fov = 3.4
#nsmooth = 25
############################### uncomment for bottomline
#path_genYarn = "../../YarnGeneration/output/woven/6x6/genYarn_NN_15000_"
#basename_in = path_genYarn + "%d.txt"
#basename_out = path_genYarn + "%d_tiled.txt"
#basename_xml = "genYarn_NN_15000.xml"
## fov = 3.8
#fov = 3.4
#nsmooth = 25

############################### uncomment for NN
#path_genYarn = "../../YarnGeneration/output/woven/6x6/genYarn_NN_15000_"
#basename_in = path_genYarn + "%d.txt"
#basename_out = path_genYarn + "%d_tiled.txt"
#basename_xml = "genYarn_NN_15000.xml"
## fov = 3.8
#fov = 3.4
#nsmooth = 25


############################### uncomment for simul
#path_simul = "../../YarnGeneration/data/woven/6x6/simul_frame_10000_"
#basename_in = path_simul + "%d.txt"
#basename_out = path_simul + "%d_tiled.txt"
#basename_xml = "simul_frame_10000.xml"
#fov = 2.6
##fov = 3.3
#nsmooth = 30


path_simul = "../../YarnGeneration/data/woven/6x6/simul_frame_6000_"
basename_in = path_simul + "%d.txt"
basename_out = path_simul + "%d_tiled.txt"
basename_xml = "simul_frame_6000.xml"
fov = 2.6
#fov = 3.3
nsmooth = 30


u = d*np.tan(0.5*np.pi*fov/180.0)
print(u)
delta = np.array([u, 1e5, u])

lb = center - delta
ub = center + delta

indices = [1, 2, 3, 4, 7, 8, 9, 10]
if True:
    for index in indices:
        fname_in = basename_in % index

        yarns = []
        with open(fname_in, "r") as fin:
            n = int(fin.readline().strip())
            for i in range(0, n):
                vertices = []
                m = int(fin.readline().strip())
                for j in range(0, m):
                    vtx = np.array([float(x) for x in fin.readline().strip().split()])
                    vertices.append(vtx)

                j = 0
                while j < m and not isInside(vertices[j], lb, ub):
                    j = j + 1
                if j < m:
                    k = j + 1
                    while k < m and isInside(vertices[k], lb, ub):
                        k = k + 1

                    assert j > 0
                    vertices[j - 1] = clip(vertices[j - 1], vertices[j], lb, ub)
                    j = j - 1

                    assert k < m
                    vertices[k] = clip(vertices[k], vertices[k - 1], lb, ub)
                    k = k + 1

                    nvertices = vertices[j : k - 1]
                    m1 = len(nvertices)
                    for j in range(0, nsmooth):
                        nvertices1 = copy.deepcopy(nvertices)
                        for k in range(0, m1):
                            if index <= 5:
                                nvertices[k][0 : 2] = 0.25*nvertices1[(k + m1 - 1) % m1][0 : 2] + 0.5*nvertices1[k][0 : 2] + 0.25*nvertices1[(k + 1) % m1][0 : 2]
                            else:
                                nvertices[k][1 : 3] = 0.25*nvertices1[(k + m1 - 1) % m1][1 : 3] + 0.5*nvertices1[k][1 : 3] + 0.25*nvertices1[(k + 1) % m1][1 : 3]

                    nvertices2 = []
                    for j in range(-ntile, ntile + 1):
                        nvertices1 = copy.deepcopy(nvertices)
                        for k in range(0, m1):
                            if index <= 5:
                                nvertices1[k][2] = nvertices1[k][2] + 2.0*u*j
                            else:
                                nvertices1[k][0] = nvertices1[k][0] + 2.0*u*j
                        nvertices2.extend(nvertices1)
                    yarns.append(nvertices2)

        if len(yarns) > 0:
            fname_out = basename_out % index
            with open(fname_out, "w") as fout:
                print >> fout, len(yarns)
                for yarn in yarns:
                    print >> fout, len(yarn)
                    for vtx in yarn:
                        print >> fout, "%.6lf %.6lf %.6lf" % (vtx[0], vtx[1], vtx[2])

with open(basename_xml, "w") as fout:
    print >> fout, "<scene version=\"0.4.4\">"
    
    for index in indices:
        if index <= 5:
            c = np.array([0.882,0.762,0.148])
        else:
            c = np.array([0.859,0.058,0.191])
        print >> fout, """<shape type="shapegroup" id="yarn_%d">
    <shape type="hair">
        <string name="filename" value="%s"/>
        <float name="radius" value="0.002"/>
        <float name="angleThreshold" value="0.001"/>
        <subsurface type="fibershader">
    			<boolean name="useRandomInteractionPoint" value="true"/>
    			<boolean name="sampleInteractionPointFromCircumference" value="false"/>
    			<fiberscat type="simpfabric5">
    				<float name="kD" value="0"/>
    				<spectrum name="colorD" value="0.99,0.99,0.99"/>
    				<spectrum name="colorR" value="0.1,0.1,0.05"/>
    				<spectrum name="colorTT" value="%.6f, %.6f, %.6f"/>
    				<float name="betaR" value="0.2"/>
    				<float name="betaTT" value="27"/>
    				<float name="gammaTT" value="38"/>
    				<float name="alpha" value="5"/>
    			</fiberscat>
    		</subsurface> 
            </shape>
            </shape>""" % (index, basename_out % index, c[0], c[1], c[2])


#        <bsdf type="roughplastic">
#            <spectrum name="diffuseReflectance" value="%.6f, %.6f, %.6f"/>
#            <float name="alpha" value="0.15"/>
#        </bsdf>




    for i in range(-ntile, ntile + 1):
        for index in indices:
            if index <= 5:
                t = np.array([2.0*u*i, 0.0, 0.0])
            else:
                t = np.array([0.0, 0.0, 2.0*u*i])

            print >> fout, """<shape type="instance">
    <ref id="yarn_%d"/>
    <transform name="toWorld">
        <translate x="%.6f" y="%.6f" z="%.6f"/>
    </transform>
</shape>""" % (index, t[0], t[1], t[2])

    print >> fout, "</scene>"
