"""
This read a the deforemed yarn and generates the average fiber

"""
from itertools import islice
hairs_path = "../../../hairs"

vrtx_num = 0

for h in range (1,30): #iterate over all hair frames
    if h<10: 
        frame_num = 'frame0000' + str(h)
    else:
        frame_num = 'frame000' + str(h)

    list_fibers = []        
    with open(hairs_path + '/scaledHairs/' + frame_num + '_scaled.txt') as f:
        fiber_cnt = list(islice(f, 1) )
        fiber_cnt = int (fiber_cnt[0])

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
            yarn_ij_x = yarn_ij_x + float(list_fibers[f][v].split()[0])
            yarn_ij_y = yarn_ij_y + float(list_fibers[f][v].split()[1])
            yarn_ij_z = yarn_ij_z + float(list_fibers[f][v].split()[2])

        yarn_file.write('\n')
        avg_x = yarn_ij_x/len(list_fibers)
        avg_y = yarn_ij_y/len(list_fibers)
        avg_z = yarn_ij_z/len(list_fibers)
        yarn_file.writelines('%.6f %.6f %.6f' % (avg_x, avg_y, avg_z) )
    yarn_file.close()
    print ('Center-fiber for frame ' + str(h) + ' is generated.')


      
