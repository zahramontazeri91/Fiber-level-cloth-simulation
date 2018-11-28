import numpy as np
import argparse
import copy


def parse_file1(filename):
    hairs=[]
    current_hair=[]
    with open(filename) as f:
        content = f.readlines()
        # content = content.split()
        # content = map(float(content))
        for data in content:
            data = data.split()
            data = list(map(float, data))
            if not len(data)==3:
                if len(current_hair) > 0:
                    hairs.append(current_hair)
                    current_hair=[]
            else:
                current_hair.append(data)
        hairs.append(current_hair)
    return hairs

def rescale_hairs(hairs):
    hairs = hairs * 4
    return hairs

def downsample_hairs(hairs, downscale):
    return hairs[:, ::downscale, :]


#stretchingMultiplier = 0.0217

def write_head(f, record_filename):
    f.write('<scene>\n')
    f.write('<description text="Nothing to say"/>\n')
    f.write('<scenetag tag="yarn_simulation"/>\n')
    f.write('<duration time="2.0"/>\n')
    f.write('<integrator type="linearized-implicit-euler" dt="0.00004" apic="1" criterion="1e-3"/>\n') #old dt="0.00004"
    f.write('<collision type="continuous-time"/>\n')
    if (0):
        f.write('<bucketinfo levelsetdistance="0.2" size="0.64" numcells="16" record_filename="{}" output_step="500"/>\n'.format(record_filename)) #yarn100
    else:
        f.write('<bucketinfo size="1.28" numcells="16" record_filename="{}" output_step="200"/>\n'.format(record_filename)) 
    f.write('<simplegravity fx="0.0" fy="0.0" fz="0.0" />\n')
    f.write('<StrandParameters>\n')   
    f.write('<poissonRatio value="0.35" />\n')
    f.write('<attachMultiplier value="1e-2" />\n')
    f.write('<density value="1.32" /> \n')
    if (0):
        f.write('<radius value="0.02" />\n') #yarn100
        f.write('<youngsModulus value="12.4e5"/>\n')  #yarn100
        f.write('<collisionMultiplier value="8.2e-5"/>\n') #yarn100
        f.write('<stretchingMultiplier value="0.0231" />\n') #yarn100
        f.write('<bendingMultiplier value="0.1314" />\n') #yarn100 
    else:
        f.write('<radius value="0.12" />\n') 
#        f.write('<youngsModulus value="9.6e5"/>\n')  #yarn4
        f.write('<youngsModulus value="8.5e5"/>\n')  #yarn8
#        f.write('<youngsModulus value="7.9e5"/>\n')  #yarn9
#        f.write('<youngsModulus value="8.8e5"/>\n')  #yarn11
#        f.write('<collisionMultiplier value="10.2e-5"/>\n')  
        f.write('<collisionMultiplier value="0.001"/>\n')  
        f.write('<stretchingMultiplier value="0.0181" />\n')
        f.write('<bendingMultiplier value="0.0278" />\n')         
    f.write('<twistingMultiplier value="1" />\n')
    f.write('<viscosity value="1e7" />\n')
    f.write('<baseRotation value="0.0"/>\n')
    f.write('<accumulateWithViscous value="0"/>\n')
    f.write('<accumulateViscousOnlyForBendingModes value="0"/>\n')
    f.write('<friction_angle value="60.0"/>\n')
    f.write('</StrandParameters>\n')
    f.write('<savetag geometry="y" twist="y" fe="y"/>\n')
#    f.write('<savetag geometry="y" twist="y" fe="y" distance_field="y"/>\n')

def write_script(f):
    
    ###stretch all fixedpoints away from the center
#    f.write('<script type="scale" x="2.5" y="2.5" z="1.0" w="0.0" start="0.0" end="2.0" group="1"/> \n')
#    f.write('<script type="scale" x="3.0" y="3.0" z="1.0" w="0.0" start="0.0" end="2.0" group="1"/> \n') #for knitted 5x10
    
    ###push on a sphere
    f.write('<script type="translate" x="0.0" y="0.0" z="-30.0" w="0.0" start="0.0" end="2.0" group="1"/> \n')   
    
    ###stretch 4 edges
#    f.write('<script type="translate" x="-10.0" y="0.0" z="0.0" w="0.0" start="0.0" end="2.0" group="1"/> \n')
#    f.write('<script type="translate" x="10.0" y="0.0" z="0.0" w="0.0" start="0.0" end="2.0" group="2"/> \n')
#    f.write('<script type="translate" x="0.0" y="-10.0" z="0.0" w="0.0" start="0.0" end="2.0" group="3"/> \n')
#    f.write('<script type="translate" x="0.0" y="10.0" z="0.0" w="0.0" start="0.0" end="2.0" group="4"/> \n')


def write_sphere(f):
	print('<distancefield type="capsule" cx="0.0" cy="0.0" cz="-9.0" rx="0.0" ry="0.0" rz="0.0" rw="1.5707963268" radius="8.0" halflength="0.0" group="5" sampled="0"/>' , file=f)
    # capsule has two hemisphere at both ends with radius and a cylinder in the middle with halflength

   
def write_cylinder(f, spacing):
	cylinder_loc=[]
	current_z = 0.0
	isup = True
	while 1:
		if current_z == 0:
			cylinder_loc.append((current_z, isup))
		else:
			cylinder_loc.append((current_z, isup))
			cylinder_loc.append((-current_z, isup))

		isup = not isup
		current_z += spacing
		if current_z > 4.5:
			break
	
	for loc in cylinder_loc:
		if loc[1]:
			# print('<distancefield type="capsule" cx="0.0" cy="0.3" cz="{}" rx="0.0" ry="0.0" rz="0.0" rw="1.5707963268" radius="0.156" halflength="0.3" group="1"/>'.format(loc[0]), file=f)
			print('<distancefield type="capsule" cx="0.0" cy="0.3" cz="{}" rx="0.0" ry="0.0" rz="0.0" rw="1.5707963268" radius="0.624" halflength="0.3" group="1"/>'.format(loc[0]), file=f)
		else:
			# print('<distancefield type="capsule" cx="0.0" cy="-0.3" cz="{}" rx="0.0" ry="0.0" rz="0.0" rw="1.5707963268" radius="0.156" halflength="0.3" group="2"/>'.format(loc[0]), file=f)
			print('<distancefield type="capsule" cx="0.0" cy="-0.3" cz="{}" rx="0.0" ry="0.0" rz="0.0" rw="1.5707963268" radius="0.624" halflength="0.3" group="2"/>'.format(loc[0]), file=f)


def write_cylinder_with_pattern(f, spacing, pattern):
    cylinder_loc=[]
    current_z = 0.0
    pattern_len = len(pattern)
    i_count = 0
    while 1:
        if current_z == 0:
            cylinder_loc.append((current_z, pattern[i_count]))
            # cylinder_loc.append((current_z, 1-pattern[i_count]))

        else:
            cylinder_loc.append((current_z, pattern[i_count]))
            cylinder_loc.append((-current_z, pattern[-i_count]))
            # cylinder_loc.append((current_z, 1-pattern[i_count]))
            # cylinder_loc.append((-current_z, 1-pattern[i_count]))
        i_count += 1
        if(i_count>=pattern_len):
            i_count = 0
        current_z += spacing
        if current_z > 4.5:
            break

  

    for loc in cylinder_loc:
        if loc[1]:
            print('<distancefield type="capsule" cx="0.0" cy="0.8" cz="{}" rx="0.0" ry="0.0" rz="0.0" rw="1.5707963268" radius="0.64" halflength="0.3" group="1"/>'.format(loc[0]), file=f)
        else:
            print('<distancefield type="capsule" cx="0.0" cy="-0.8" cz="{}" rx="0.0" ry="0.0" rz="0.0" rw="1.5707963268" radius="0.64" halflength="0.3" group="2"/>'.format(loc[0]), file=f)


def write_hair_collections(f, hair_counter):
	start_id = 0
	for i in range(len(hair_counter)):
		count = hair_counter[i]
		f.write('<hair params="0" start="'+str(start_id)+'" count="'+str(count)+'"/>\n')
		start_id += count


def write_end(f):
    print('writing scene is done!')
    f.write('</scene>\n')


def write_fibers(f, hairs, hair_counter, fid, is_tip_fixed):
	for hair in hairs:
		hair_counter = write_single_hair(f, hair, hair_counter, fid, is_tip_fixed)
	return hair_counter


def write_fibers_given_fixes(f, hairs, hair_counter, fid, fixed_idx, fixed_group):
	for hair in hairs:
		hair_counter = write_single_hair_given_fixes(f, hair, hair_counter, fid, fixed_idx, fixed_group)
	return hair_counter

def write_single_hair_given_fixes(f, hair, hair_counter, fid, fixed_idx, fixed_group):
    for i in range(hair.shape[0]):
        if (i in fixed_idx):
            group = fixed_group[ fixed_idx.tolist().index(i) ]
            f.write('<particle x="'+str(hair[i, 0])+' '+str(hair[i, 1])+' '+str(hair[i, 2])+'" v="0.0 0.0 0.0" fixed="1" group="'+str(group)+'" fid="'+str(fid)+'"/>\n')
        else:
            f.write('<particle x="'+str(hair[i, 0])+' '+str(hair[i, 1])+' '+str(hair[i, 2])+'" v="0.0 0.0 0.0" fixed="0" fid="'+str(fid)+'"/>\n')

    f.write('\n')

    hair_counter.append(hair.shape[0])
    return hair_counter

def write_single_hair(f, hair, hair_counter, fid, is_tip_fixed):
    for i in range(hair.shape[0]):
        if (is_tip_fixed and (i==0)):
            f.write('<particle x="'+str(hair[i, 0])+' '+str(hair[i, 1])+' '+str(hair[i, 2])+'" v="0.0 0.0 0.0" fixed="1" group="3" fid="'+str(fid)+'"/>\n')
        elif (is_tip_fixed and (i==hair.shape[0]-1)):
            f.write('<particle x="'+str(hair[i, 0])+' '+str(hair[i, 1])+' '+str(hair[i, 2])+'" v="0.0 0.0 0.0" fixed="1" group="4" fid="'+str(fid)+'"/>\n')
        else:
            f.write('<particle x="'+str(hair[i, 0])+' '+str(hair[i, 1])+' '+str(hair[i, 2])+'" v="0.0 0.0 0.0" fixed="0" fid="'+str(fid)+'"/>\n')

    f.write('\n')

    hair_counter.append(hair.shape[0])
    return hair_counter

def centralizehair(hair):
	mean_hair = np.mean(hair, axis=0)
	center = np.mean(mean_hair, axis=0)
	center_hair = copy.deepcopy(hair)
	center_hair[:, :, 0] -= center[0]
	center_hair[:, :, 1] -= center[1]
	center_hair[:, :, 2] -= center[2]
	return center_hair


def rotatehairs(hairs, theta_x, theta_y, theta_z):
	Rz = np.array([[np.cos(theta_z), -np.sin(theta_z), 0], [np.sin(theta_z), np.cos(theta_z), 0], [0, 0, 1]])
	Ry = np.array([[np.cos(theta_y), 0, np.sin(theta_y)], [0, 1, 0], [-np.sin(theta_y), 0, np.cos(theta_y)]])
	Rx = np.array([[1, 0, 0], [0, np.cos(theta_x), -np.sin(theta_x)], [0, np.sin(theta_x), np.cos(theta_x)]])

	R = Rz.dot(Ry).dot(Rx)
	mhairs = copy.deepcopy(hairs)

	for hair in mhairs:
		for p in hair:
			p0 = p.reshape(3,1)
			p0 = R.dot(p0)
			p[0] = p0[0]
			p[1] = p0[1]
			p[2] = p0[2]
	return mhairs

def translatehairs(hairs, dx, dy, dz):
	mhairs = copy.deepcopy(hairs)
	for hair in mhairs:
		for p in hair:
			p[0] += dx
			p[1] += dy
			p[2] += dz
	return mhairs

def generate_hair1():
    numv = 100
    length = 11.96
    # length = 12.12
    hairs = np.zeros((numv, 3))
    zs = np.linspace(-length/2, length/2, numv)
    hairs[:, 2] = zs
    return hairs

def parse_obj(filename):
    vertices = []
    edges = []
    fibers = []
    with open(filename, 'r') as fin:
        for line in fin:
            data = line.split()
            if len(data)==0:
                break;
            if data[0] == 'v':
                data = list(map(float, data[1:]))
                assert(len(data) == 3 or len(data) == 4)
                vertices.append(data)
            elif data[0] == 'l':
                data = list(map(int, data[1:]))
                data[0] -= 1
                data[1] -= 1
                assert(len(data) == 2)
                edges.append(data)
                pass
            else:
                break

    # now create individual fiber
    for e in edges:
        ll = len(fibers)
        for ii in range(ll-1,-1,-1):
            if fibers[ii][-1] == e[0]:
                fibers[ii].append( e[1] )
                break
        else:
            fibers.append( [e[0], e[1]] )

    #check consistency
    for ii in range(1, len(fibers)):
        assert len( fibers[ii] ) == len(fibers[0])

    num_fibers = len(fibers)
    num_v_per_fiber = len(fibers[0])

    print(num_fibers, num_v_per_fiber)

    hairs = np.zeros((num_fibers, num_v_per_fiber, 3))
    # print(fibers[0][99])
    for i in range(num_fibers):
        for j in range(num_v_per_fiber):
            # print(i,j, vertices[fibers[i][j]])
            hairs[i, j, :] = vertices[fibers[i][j]][:3]
    print(num_fibers, 'strand loaded')
    return hairs


def main_batch(args):
    patterns = [[0,0,0,1,1], [1,1,1,1,0], [1,0,1,0,0], [1,0]]
    names = ['00011', '11110', '10100', '10']
    spacings = [0.3, 0.4, 0.5, 0.6, 0.7, 0.8]
    spacing_names = ['0.5x','1.0x','1.5x','2.0x','2.5x','3.0x']

    outbase = args['output']

    for (i, pattern) in enumerate(patterns):
        for (j, spacing) in enumerate(spacings):
            outname = outbase+names[i]+'_'+spacing_names[j]+'.xml'
            print(outname)
            
            hair1 = generate_hair1()
            hair1 = hair1.reshape((1,-1,3))
            if args['input']:
                hair1 = parse_obj(args['input'])
            hair_counter=[]
            f = open(outname, 'w+')
            record_name = 'results/yarn11/{}/{}/'.format(names[i], spacing_names[j])
            # print(record_name)
            write_head(f, record_filename=record_name)
            # write_head(f, record_filename='results/speed/1.5_10/1.0_both/')
            write_cylinder_with_pattern(f, spacing=spacing, pattern=pattern)
            write_script(f)
            hair_counter = write_fibers(f=f, hairs=hair1, hair_counter=hair_counter, fid=0, is_tip_fixed=True)
            write_hair_collections(f, hair_counter)
            write_end(f)


# In[]
#simulConfig = 'slipstitchrib_stretch'   
simulConfig = 'slipstitchrib_push' 
output_file = simulConfig + '.xml'
hair_file = 'tileKnit/slipstitchrib.obj'
fixed_file = 'tileKnit/fixedPnts.txt'

hair1 = generate_hair1()
hair1 = hair1.reshape((1,-1,3))
hair1 = parse_obj(hair_file)
data = np.loadtxt(fixed_file, dtype='int')
fixed_idx = data[:,0]
fixed_group = data[:,1]

hair_counter=[]
f = open(output_file, 'w+')
write_head(f, record_filename='results/'+ simulConfig + '/')
write_script(f)
write_sphere(f)
hair_counter = write_fibers_given_fixes(f=f, hairs=hair1, hair_counter=hair_counter, fid=0, fixed_idx=fixed_idx, fixed_group=fixed_group)
write_hair_collections(f, hair_counter)
write_end(f)
