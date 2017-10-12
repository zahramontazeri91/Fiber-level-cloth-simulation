"""
This read a the deforemed yarn and generates the average fiber

"""

from itertools import islice
hairs_path = "../data/hairs/"

with open(hairs_path + 'frame00001_hairs.txt') as f:
    num_vertices = int(f.readline()) #read number of vertices 

for h in range (1,30): #iterate over all hair frames
    if h<10: 
        frame_num = 'frame0000' + str(h)
    else:
        frame_num = 'frame000' + str(h)

    list_fibers = []        
    with open(hairs_path + frame_num + '_hairs.txt') as f:
        while True:
            next_n_lines = list(islice(f, int(num_vertices + 1)) )
            list_fibers.append(next_n_lines)
            if not next_n_lines:
                break
    
    yarn_file = open(hairs_path + frame_num + '_yarn.txt','w') 
    
    yarn_file.write(str(1480) ) 
    for v in range (1, num_vertices): #iterate over vertices 
        yarn_ij_x = 0.0
        yarn_ij_y = 0.0
        yarn_ij_z = 0.0      
        for f in range (0,len(list_fibers)-1): #iterate over fibers
               yarn_ij_x = (yarn_ij_x + float(list_fibers[f][v].split()[0]) )/len(list_fibers)
               yarn_ij_y = (yarn_ij_y + float(list_fibers[f][v].split()[1]) )/len(list_fibers)
               yarn_ij_z = (yarn_ij_z + float(list_fibers[f][v].split()[2]) )/len(list_fibers)
               
        yarn_file.write('\n')
        yarn_file.writelines([str(yarn_ij_x), ' ', str(yarn_ij_y), ' ', str(yarn_ij_z)])
    yarn_file.close()
    print ('Center-fiber for frame ' + str(h) + ' is generated.')


      
