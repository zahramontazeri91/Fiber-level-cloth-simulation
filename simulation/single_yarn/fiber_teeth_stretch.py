import numpy as np
import argparse
import copy


def parse_fiber_bundle(filename):
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

def rescale_hairs(hairs, scale=4):
    hairs = hairs * scale
    return hairs

def downsample_hairs(hairs, downscale):
    return hairs[:, ::downscale, :]



def write_head(f, record_filename):
	f.write('<scene>\n')
	f.write('<description text="Nothing to say"/>\n')
	f.write('<scenetag tag="recording"/>\n')
	f.write('<duration time="2.0"/>\n')
	f.write('<integrator type="linearized-implicit-euler" dt="0.00004" apic="1" criterion="1e-3"/>\n')
	f.write('<collision type="continuous-time"/>\n')               
	# f.write('<bucketinfo size="0.64" numcells="16" record_filename="results/spacing/train/spacing3.5x/fiber/"/>\n') #spacing change
	f.write('<bucketinfo size="0.64" numcells="16" record_filename="{}" output_step="100"/>\n'.format(record_filename))  #pattern change
	f.write('<simplegravity fx="0.0" fy="0.0" fz="0.0" />\n')
	f.write('<StrandParameters>\n')
	f.write('<radius value="0.002" />\n')
	f.write('<youngsModulus value="9.6e5"/>\n')
	f.write('<poissonRatio value="0.35" />\n')
	f.write('<collisionMultiplier value="1e-3"/>\n')
	f.write('<attachMultiplier value="1e-3" />\n') 
	f.write('<density value="1.32" /> \n')
	f.write('<viscosity value="1e7" />\n')
	f.write('<baseRotation value="0.0"/>\n') 
	f.write('<accumulateWithViscous value="0"/>\n')
	f.write('<accumulateViscousOnlyForBendingModes value="0"/>\n')
	f.write('<friction_angle value="60.0"/>\n')
	f.write('</StrandParameters>\n')
	f.write('<savetag geometry="y"/>\n')

def write_cylinder2(f, spacing):
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
			print('<distancefield type="capsule" cx="0.0" cy="0.3" cz="{}" rx="0.0" ry="0.0" rz="0.0" rw="1.5707963268" radius="0.1" halflength="0.3" group="1"/>'.format(loc[0]), file=f)
		else:
			print('<distancefield type="capsule" cx="0.0" cy="-0.3" cz="{}" rx="0.0" ry="0.0" rz="0.0" rw="1.5707963268" radius="0.1" halflength="0.3" group="2"/>'.format(loc[0]), file=f)

def write_cylinder(f):
	f.write('<distancefield type="capsule" cx="0.0" cy="0.3" cz="0.0" rx="0.0" ry="0.0" rz="0.0" rw="1.5707963268" radius="0.1" halflength="0.3" group="1"/>\n')
	f.write('<distancefield type="capsule" cx="0.0" cy="0.3" cz="-3.0" rx="0.0" ry="0.0" rz="0.0" rw="1.5707963268" radius="0.1" halflength="0.3" group="1"/>\n')
	f.write('<distancefield type="capsule" cx="0.0" cy="0.3" cz="3.0" rx="0.0" ry="0.0" rz="0.0" rw="1.5707963268" radius="0.1" halflength="0.3" group="1"/>\n')
	f.write('<distancefield type="capsule" cx="0.0" cy="-0.3" cz="1.5" rx="0.0" ry="0.0" rz="0.0" rw="1.5707963268" radius="0.1" halflength="0.3" group="2"/>\n')
	f.write('<distancefield type="capsule" cx="0.0" cy="-0.3" cz="-1.5" rx="0.0" ry="0.0" rz="0.0" rw="1.5707963268" radius="0.1" halflength="0.3" group="2"/>\n')
	f.write('<distancefield type="capsule" cx="0.0" cy="-0.3" cz="4.5" rx="0.0" ry="0.0" rz="0.0" rw="1.5707963268" radius="0.1" halflength="0.3" group="2"/>\n')
	f.write('<distancefield type="capsule" cx="0.0" cy="-0.3" cz="-4.5" rx="0.0" ry="0.0" rz="0.0" rw="1.5707963268" radius="0.1" halflength="0.3" group="2"/>\n')



	# space = 1.1
	# for i in range(1,2):
		# f.write('<distancefield type="capsule" cx="0.0" cy="0.35" cz="'+str(i*space)+'" rx="0.0" ry="0.0" rz="0.0" rw="1.5707963268" radius="0.1" halflength="0.3" group="1"/>\n')
		# f.write('<distancefield type="capsule" cx="0.0" cy="0.35" cz="'+str(-i*space)+'" rx="0.0" ry="0.0" rz="0.0" rw="1.5707963268" radius="0.1" halflength="0.3" group="1"/>\n')
	# 	# f.write('<distancefield type="capsule" cx="0.0" cy="-0.35" cz="'+str(i*space+space/2)+'" rx="0.0" ry="0.0" rz="0.0" rw="1.5707963268" radius="0.1" halflength="0.3" group="2"/>\n')
	# 	# f.write('<distancefield type="capsule" cx="0.0" cy="-0.35" cz="'+str(-i*space-space/2)+'" rx="0.0" ry="0.0" rz="0.0" rw="1.5707963268" radius="0.1" halflength="0.3" group="2"/>\n')
	return

def write_cylinder_with_pattern(f, spacing, pattern):
    cylinder_loc=[]
    current_z = 0.0
    pattern_len = len(pattern)
    i_count = 0
    while 1:
        if current_z == 0:
            cylinder_loc.append((current_z, pattern[i_count]))
        else:
            cylinder_loc.append((current_z, pattern[i_count]))
            cylinder_loc.append((-current_z, pattern[-i_count]))
        i_count += 1
        if(i_count>=pattern_len):
            i_count = 0
        current_z += spacing
        if current_z > 4.5:
            break

    for loc in cylinder_loc:
        if loc[1]:
            # print('<distancefield type="capsule" cx="0.0" cy="0.3" cz="{}" rx="0.0" ry="0.0" rz="0.0" rw="1.5707963268" radius="0.156" halflength="0.3" group="1"/>'.format(loc[0]), file=f)
            print('<distancefield type="capsule" cx="0.0" cy="0.8" cz="{}" rx="0.0" ry="0.0" rz="0.0" rw="1.5707963268" radius="0.64" halflength="0.3" group="1"/>'.format(loc[0]), file=f)
        else:
            # print('<distancefield type="capsule" cx="0.0" cy="-0.3" cz="{}" rx="0.0" ry="0.0" rz="0.0" rw="1.5707963268" radius="0.156" halflength="0.3" group="2"/>'.format(loc[0]), file=f)
            print('<distancefield type="capsule" cx="0.0" cy="-0.8" cz="{}" rx="0.0" ry="0.0" rz="0.0" rw="1.5707963268" radius="0.64" halflength="0.3" group="2"/>'.format(loc[0]), file=f)


def write_script(f):
    f.write('<script type="translate" x="0.0" y="-0.8" z="0.0" w="0.0" start="0.0" end="1.0" group="1"/>\n')
    f.write('<script type="translate" x="0.0" y="0.8" z="0.0" w="0.0" start="0.0" end="1.0" group="2"/>\n')
    f.write('<script type="translate" x="0.0" y="0.0" z="-1.0" w="0.0" start="1.1" end="2.0" group="3"/>\n')
    f.write('<script type="translate" x="0.0" y="0.0" z="1.0" w="0.0" start="1.1" end="2.0" group="4"/>\n')   

    return	

def write_hair_collections(f, hair_counter):
	start_id = 0
	for i in range(len(hair_counter)):
		count = hair_counter[i]
		f.write('<hair params="0" start="'+str(start_id)+'" count="'+str(count)+'"/>\n')
		start_id += count
	

def write_end(f):
    f.write('</scene>\n')


def write_fibers(f, hairs, hair_counter, fid, is_tip_fixed):
	for hair in hairs:
		hair_counter = write_single_hair(f, hair, hair_counter, fid, is_tip_fixed)
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



def write_strand(f, s, fixed_s):
	counts = 0
	start = 0
	for strand in s:
		for p in strand:
			f.write('<particle x="'+str(p[0])+' '+str(p[1])+' '+str(p[2])+'" v="0.0 0.0 0.0" fixed="0"/>\n')

	for strand in fixed_s:
		for i, p in enumerate(strand):
			if i==0 or i==strand.shape[0]-1:
				f.write('<particle x="'+str(p[0])+' '+str(p[1])+' '+str(p[2])+'" v="0.0 0.0 0.0" fixed="1"/>\n')
			else:
				f.write('<particle x="'+str(p[0])+' '+str(p[1])+' '+str(p[2])+'" v="0.0 0.0 0.0" fixed="0"/>\n')

	for strand in s: 
		counts = strand.shape[0]
		f.write('<hair params="0" start="'+str(start)+'" count="'+str(counts)+'"/>\n')
		start += counts

	for strand in fixed_s: 
		counts = strand.shape[0]
		f.write('<hair params="0" start="'+str(start)+'" count="'+str(counts)+'"/>\n')
		start += counts

def centralizehair(hair):
	mean_hair = np.mean(hair, axis=0)
	center = np.mean(mean_hair, axis=0)
	center_hair = copy.deepcopy(hair)
	center_hair[:, :, 0] -= center[0]
	center_hair[:, :, 1] -= center[1]
	center_hair[:, :, 2] -= center[2]
	return center_hair, center


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

if __name__ == '__main__':
	parser = argparse.ArgumentParser()
	parser.add_argument("-i", "--input", required=True, help="path to input fiber pattern")
	parser.add_argument("-o", "--output", required=True, help="output file name")
	args = vars(parser.parse_args())
	input_file = args['input']
	output_file = args['output']
	hairs =  parse_fiber_bundle(input_file)
	print(len(hairs))
	nphair = np.array(hairs)
	nphair = rescale_hairs(nphair)
	hair1, meanval = centralizehair(nphair)
	hair_counter=[]
	f = open(output_file, 'w+')
	write_head(f, record_filename='results/single_yarn/yarn11_fiber_1.6/')
	write_cylinder_with_pattern(f, spacing=2.56, pattern=[1,0])
	#write_cylinder_with_pattern(f, spacing=1.92, pattern=[0,0,1,1,0])
	write_script(f)
	hair_counter = write_fibers(f=f, hairs=hair1, hair_counter=hair_counter, fid=0, is_tip_fixed=True)
	write_hair_collections(f, hair_counter)
	write_end(f)










# spacing0.5x:0.3
# spacing1.0x:0.4  
# spacing1.5x:0.5 
# spacing2.0x:0.6 
# spacing2.5x:0.7 
# spacing3.0x:0.8 

