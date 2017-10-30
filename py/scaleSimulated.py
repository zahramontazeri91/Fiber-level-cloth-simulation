"""
This read a the simulated data and devide the values by 4

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
    
    fiber_num = len(list_fibers)
    yarn_file = open(hairs_path + '/scaledHairs/' + frame_num + '_scaled_1ply.txt','w') 
    
    yarn_file.write( '%s \n' % str(fiber_num/2) )     
    for f in range (0,fiber_num/2):
        yarn_file.write( '%s \n' % str(vrtx_num) )  
        for v in range (0,vrtx_num):
            scale_x = float(list_fibers[f][v].split()[0])*0.25 
            scale_y = float(list_fibers[f][v].split()[1])*0.25 
            scale_z = float(list_fibers[f][v].split()[2])*0.25 
            yarn_file.writelines('%.6f %.6f %.6f \n' % (scale_x, scale_y, scale_z) )  
         
    yarn_file.close()
    print ('Center-fiber for frame ' + str(h) + ' is generated.')


      
