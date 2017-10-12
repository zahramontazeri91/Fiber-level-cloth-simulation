"""
This read a the deforemed yarn and generates the average fiber

"""
from itertools import islice
hairs_path = "../../../hairs"

vrtx_num = 0

for h in range (29,30): #iterate over all hair frames
    if h<10: 
        frame_num = 'frame0000' + str(h)
    else:
        frame_num = 'frame000' + str(h)

    list_fibers = []        
    with open(hairs_path + '/hairs/' + frame_num + '_hairs.txt') as f:
        while True:
            cnt = list(islice(f, 1) )
            if not cnt:
                break
            vrtx_num = int( cnt[0] )
            next_n_lines = list(islice(f, vrtx_num) )
            list_fibers.append(next_n_lines)
    
    
    yarn_file = open(hairs_path + '/avgHairs/' + frame_num + '_avg.txt','w') 
    
    yarn_file.write(str(vrtx_num) ) 
    for v in range (0, vrtx_num): #iterate over vertices 
        yarn_ij_x = 0.0
        yarn_ij_y = 0.0
        yarn_ij_z = 0.0      
        for f in range (0,len(list_fibers)): #iterate over fibers
            yarn_ij_x = yarn_ij_x + float(list_fibers[f][v].split()[0])*0.25 
            yarn_ij_y = yarn_ij_y + float(list_fibers[f][v].split()[1])*0.25
            yarn_ij_z = yarn_ij_z + float(list_fibers[f][v].split()[2])*0.25

        yarn_file.write('\n')
        yarn_file.writelines([str(yarn_ij_x/len(list_fibers) ), ' ', str(yarn_ij_y/len(list_fibers) ), ' ', str(yarn_ij_z/len(list_fibers) )])
    yarn_file.close()
    print ('Center-fiber for frame ' + str(h) + ' is generated.')


      
