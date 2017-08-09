MADEBY = "This 3D code was produced by Vivid"





def obj_exporter(surf, filename):
    f = open(filename, 'w')
    f.write('# '+MADEBY +'\n\n\n')
    f.write('o ' + surf['label'] +'\n')
    surf_vertex = surf['vertres_filtered']
    for i in range(surf_vertex.shape[0]):
        f.write('v ' + '%.7g %.7g %.7g \n' % (surf_vertex[i][0], surf_vertex[i][1], surf_vertex[i][2]))

    for i in range(len(surf['faces_filtered'])):
        f.write('f ' + ' '.join(map(lambda a: '%d' % a, surf['faces_filtered'][i]))+'\n')

    f.close
    
def obj_and_mtl_exporter(surf, filename, vmin, vmax):
    import matplotlib.pyplot as plt
    import matplotlib.colors as colors

    if vmin == None:
        vmin = min(surf['quen']['list_quen'])
    if vmax == None:
        vmax = max(surf['quen']['list_quen'])

    f = open(filename+'.obj', 'w')
    m = open(filename+'.mtl', 'w')
    f.write('# '+MADEBY +'\n\n\n')
    f.write('mtllib '+filename+'.mtl\n')
    f.write('o ' + surf['label'] +'\n')
    surf_vertex = surf['vertres_filtered']
    for i in range(surf_vertex.shape[0]):
        f.write('v ' + '%.7g %.7g %.7g \n' % (surf_vertex[i][0], surf_vertex[i][1], surf_vertex[i][2]))
    
    cmap = plt.cm.get_cmap('cool') #TODO find a way to get 'Plasma' color map to work
    #norm = matplotlib.colors.Normalize(vmin = vmin, vmax = vmax) #TODO make it possible so the user can choose
    norm =colors.LogNorm(vmin = vmin, vmax = vmax)


    color_d = {'list_color': []}
    for i in surf['quen']['list_quen']:
        color = cmap(norm(i))
        color = '%0.03g %0.03g %0.03g' %(color[0], color[1], color[1])
        if color in color_d:
            color_d[color].append(i)
        else:
            color_d['list_color'].append(color)
            color_d[color] = [i, ]

    counter = 0
    for i in color_d['list_color']: #TODO change i to i_colors

        m.write('newmtl '+ str(counter)+'\n' +\
                'Ns 96.078\n'+\
                'Ka 1.000 1.000 1.000 \n' +\
                'Kd '+ i +'\n' +\
                'Ks 0.000 0.000 0.000\n'+\
                'Ni 1.000000\n'+\
                'd 1.0\n'+\
                'illum 0\n'+\
                'em 0.000000\n\n\n')
        f.write('usemtl '+ str(counter)+'\n')
        for k in color_d[i]:
            for j in surf['quen'][str(k)]:
                f.write('f ' + ' '.join(map(lambda a: '%d' % a, surf['faces_filtered'][j]))+'\n')
        counter += 1
    f.close
    m.close
    
def vornoy_emporter_elad(surcename):
    import numpy as np
    f = open(surcename, 'r')
    points_len = int(f.readline())
    array_points = np.zeros((points_len,3), dtype =np.float64, order='C')
    for i in range(points_len):
        line = f.readline().split()
        array_points[i]=np.asfarray(line)
    faces_len = int(f.readline())
    array_faces = np.empty((faces_len), dtype = np.dtype(object), order='C')
    array_faces_indx = np.empty((faces_len), dtype= np.dtype(object), order='C')
    for i in range(faces_len):
        line = f.readline().split()
        if line == []:
            print 'line empty'
        array_faces[i] = tuple(map(lambda p: int(p), line[1:1 + int(line[0])]))
        lst = list(array_faces[i])
        for j in range(len(lst)):
            lst[j]+= 1
        array_faces[i]=tuple(lst)
        array_faces_indx[i] = (int(line[-1]), int(line[-2]))
        if int(line[-1]) == 1017610 or int(line[-2]) == 1017610:
            print 'ERROW', i, line
    f.close()
    output = {'vertres':array_points, 'faces':array_faces, 'faces_indx':array_faces_indx}
    return output


def clear_faces_mask(surf, mask):
    import numpy as np
    list_faces = []
    mask_len = mask.size
    for i in range(surf['faces'].size):
        # TODO consulte with tomer about upgrading
        f_i = surf['faces_indx'][i][0]
        f_j = surf['faces_indx'][i][1]
        if (mask_len <= f_i and mask_len <= f_j):
            pass
        elif (mask_len <= f_i or mask_len <= f_j):
            if (mask_len <= f_i and mask[f_j] == True):
                 list_faces.append(surf['faces'][i])
            elif (mask_len <= f_j and mask[f_i] == True):
                 list_faces.append(surf['faces'][i])
        elif mask[f_i] != mask[f_j]:
            list_faces.append(surf['faces'][i])


    array_faces = np.array(list_faces)
    surf['faces_filtered'] = array_faces
    return surf

def clear_faces_color(surf, mask, quen):
    import numpy as np
    list_faces = []
    quen_d = {}
    quen_d['list_quen'] = []
    mask_len = mask.size
    counter = 0
    for i in range(surf['faces'].size):
        # TODO consulte with tomer about upgrading
        f_i = surf['faces_indx'][i][0]
        f_j = surf['faces_indx'][i][1]
        #if (mask_len <= f_i and mask_len <= f_j):
        #    pass
        if (mask_len <= f_i or mask_len <= f_j):
            if (mask_len <= f_i and mask[f_j] == True):
                list_faces.append(surf['faces'][i])
                color = np.float64(quen[f_j])
                if str(color) in quen_d:
                    quen_d[str(color)].append(counter)
                else:
                    quen_d['list_quen'].append(color)
                    quen_d[str(color)] = [counter, ]
                    counter += 1
            elif (mask_len <= f_j and mask[f_i] == True):
                list_faces.append(surf['faces'][i])
                color = np.float64(quen[f_i])
                if str(color) in quen_d:
                    quen_d[str(color)].append(counter)
                else:
                    quen_d['list_quen'].append(color)
                    quen_d[str(color)] = [counter, ]
                    counter += 1
        elif mask[f_i] != mask[f_j]:
            list_faces.append(surf['faces'][i])
            color = (np.float64(quen[f_i]) + np.float64(quen[f_j])) / 2
            if str(color) in quen_d:
                quen_d[str(color)].append(counter)
            else:
                quen_d['list_quen'].append(color)
                quen_d[str(color)] = [counter, ]
                counter += 1
    quen_d['list_quen']= np.array(quen_d['list_quen'])
    surf['quen'] = quen_d
    array_faces = np.array(list_faces)
    surf['faces_filtered'] = array_faces
    return surf

def clear_points(surf):
    import numpy as np
    used_points = {}
    for i in range(len(surf['faces_filtered'])):
        for j in range(len(surf['faces_filtered'][i])):
            if surf['faces_filtered'][i][j] not in used_points:
                used_points[str(surf['faces_filtered'][i][j])] = True
    list_points = []
    old_new = {}
    counter = 1
    for i in range(surf['vertres'].size):
        if str(i+1) in used_points:
            list_points.append(surf['vertres'][i])
            old_new[str(i+1)] = counter
            counter += 1
    print surf['faces_filtered'].size, len(surf['faces_filtered'])
    print len(surf['faces_filtered'][0])
    print surf['faces_filtered']
    for i in range(len(surf['faces_filtered'])):
        face = []
        for j in range(len(surf['faces_filtered'][i])):
            face.append(old_new[str(surf['faces_filtered'][i][j])])
        surf['faces_filtered'][i] = tuple(face)
    array_points = np.zeros((len(list_points), 3), dtype=np.float64, order='C')
    for i in range(len(list_points)):
        array_points[i] = list_points[i]
    surf['vertres_filtered'] = array_points
    return surf

def write_to_vorn(mesh_par, box_R, vornoi_input):
    import pickle

    f = open(vornoi_input, 'w')
    f.write('-%.7g -%.7g -%.7g \n' % (box_R, box_R, box_R))
    f.write('%.7g %.7g %.7g \n' % (box_R, box_R, box_R))

    for i in range(len(mesh_par)):
        f.write('%.7g %.7g %.7g \n' % (mesh_par[i][0], mesh_par[i][1], mesh_par[i][2]))
    
    f.close

def run_vorn(vornoi_input, vornoi_output):
    import os
    os.system(r'D:\alpa\OBJ\Rich3D.exe '+vornoi_input+' '+vornoi_output)

def full_func_no_color(mesh_par, mask, box_R, output_file, label='vivid 3d obj'):

    vornoi_input = r'code_used_files\_vornoi_input'
    vornoi_output = r'code_used_files\_vornoi_output'
    
    write_to_vorn(mesh_par, box_R, vornoi_input)
    run_vorn(vornoi_input, vornoi_output)
    surf = vornoy_emporter_elad(vornoi_output)
    surf['label']= label
    surf = clear_faces_mask(surf, mask) #clear_faces_color(surf, mask, quen)
    surf = clear_points(surf)
    obj_exporter(surf, output_file)

def full_func_color(mesh_par, mask, box_R, quen, output_file,label='vivid 3d obj', vmin = None, vmax = None):   #, quen , vmin = None, vmax = None

    vornoi_input = r'code_used_files\_vornoi_input'
    vornoi_output = r'code_used_files\_vornoi_output'

    write_to_vorn(mesh_par, box_R, vornoi_input)
    run_vorn(vornoi_input, vornoi_output)
    surf = vornoy_emporter_elad(vornoi_output)
    surf['label']= label
    #surf = clear_faces_mask(surf, mask)
    surf = clear_faces_color(surf, mask, quen)
    surf = clear_points(surf)
    #obj_exporter(surf, output_file)
    obj_and_mtl_exporter(surf, output_file, vmin, vmax)


    
#def second_vorn():
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
