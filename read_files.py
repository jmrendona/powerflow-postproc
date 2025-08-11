from scipy.io import netcdf
from copy import deepcopy
import numpy as np
import sys
from collections import OrderedDict
import pdb
import matplotlib.pyplot as plt 
import os

def read_mass_flux(directory,probe_file):
    directory = directory
    probe_file = probe_file

    f1 = netcdf.netcdf_file('{0:s}/{1:s}'.format(directory,probe_file), 'r')
    measurements = f1.variables['measurements'][()]
    st1 = f1.variables['start_time'][()]
    et1 = f1.variables['end_time'][()]

    lx_scales = f1.variables['lx_scales']
    lx_offsets = f1.variables['lx_offsets']
    face_area = f1.variables['face_area'][()]
    face_ids = f1.variables['face'][()]


    if sys.version_info[0] < 3:
        vnames = f1.variables['variable_short_names']
        variable_names = deepcopy(vnames[0: -1]) # copy without last 0
        variable_names[variable_names == ''] = ' ' # replace 0 by space
        var_list = ''.join(variable_names).split(' ') # join and split
    else:
        x = f1.variables['variable_short_names'][()]
        var_list = x.tostring().decode('utf-8').split('\x00')[:-1]

    print(var_list)

    if sys.version_info[0] < 3:
        fnames = f1.variables['face_names']
        face_names = deepcopy(fnames[0: -1]) # copy without last 0
        face_names[face_names == ''] = ' ' # replace 0 by space
        face_list = ''.join(face_names).split(' ') # join and split
    else:
        x = f1.variables['face_names'][()]
        face_list = x.tostring().decode('utf-8').split('\x00')[:-1]

    print(face_list)

    print("Scaling factors")
    coeff_dx = lx_scales[0]
    CFL_number = lx_scales[-1]
    coeff_dt = lx_scales[0] / lx_scales[3]
    dt = CFL_number * coeff_dt
    coeff_press = lx_scales[1] * lx_scales[3]**2
    coeff_vel = lx_scales[3]
    coeff_density = lx_scales[1]
    offset_pressure = lx_offsets[4]

    print("weight factor")
    t_lat = 1.0/3.0 # lattive temperature
    r_lat = 1.0 # lattice r_gas constant
    rTlat = r_lat*t_lat # weight for density to pressure
    irTlat = 1.0/rTlat # weight for pressure to density

    # Parse var name to find index
    ipress=-1
    irho=-1
    ivx=-1
    for iv,var in enumerate(var_list):
        if var=='x_velocity':
            ivx=iv
        elif var=='y_velocity':
            ivy=iv
        elif var=='z_velocity':
            ivz=iv
        elif var=='static_pressure':
            ipress=iv
        elif var=='density':
            irho=iv
        elif var=='mass_flux':
            imflux=iv
    var_list=[]
    var_idx=[]
    if ipress==-1:
        var_list.append('density')
        var_idx.append(irho)
    else:
        var_list.append('static_pressure')
        var_idx.append(ipress)
    if ivx>=0:
        var_list.append('x_velocity')
        var_list.append('y_velocity')
        var_list.append('z_velocity')
        var_idx.append(ivx)
        var_idx.append(ivy)
        var_idx.append(ivz)
    if imflux>=0:
        var_list.append('mass_flux')
        var_idx.append(imflux)



    for num,f_id in enumerate(face_ids):
        fname = face_list[f_id]
        #print(fname)
        export_name = '{0:s}.txt'.format(probe_file.replace('.csnc',''),fname)

        data = OrderedDict()

        data['time'] = 0.5*(st1+et1)*dt

        if ipress==-1:
            data[var_list[0]]=measurements[:,var_idx[0],num] * coeff_density
            data['static_pressure']=(measurements[:,var_idx[0],num] * rTlat + offset_pressure) * coeff_press
        else:
            data[var_list[0]]=(measurements[:,var_idx[0],num]+offset_pressure)*coeff_press
            data['density']=measurements[:,var_idx[0],num] * irTlat * coeff_density

        if ivx>=0:
            data[var_list[1]]=measurements[:,var_idx[1],num] * coeff_vel
            data[var_list[2]]=measurements[:,var_idx[2],num] * coeff_vel
            data[var_list[3]]=measurements[:,var_idx[3],num] * coeff_vel

        if imflux>=0:
            data[var_list[-1]]=measurements[:,var_idx[-1],num] * coeff_density * coeff_vel * coeff_dx**2*face_area[num]


        header = ' '.join(data.keys())
        nf = len(data.keys())
        nt = data['time'].size
        X = np.zeros((nt,nf))
        for iff,fn in enumerate(data.keys()):
            X[:,iff] = data[fn]
        
        #np.savetxt(export_name,X,header=header, comments='#')
        np.savetxt(os.path.join(directory, export_name),X,header=header, comments='#')
        
    f1.close()
    return(directory + export_name)

def read_probe_file(directory,probe_file):
     
    directory = directory
    probe_file = probe_file
    
    f1 = netcdf.netcdf_file('{0:s}/{1:s}'.format(directory,probe_file), 'r')
    meas1 = f1.variables['measurements']
    st1 = f1.variables['start_time'][()]
    et1 = f1.variables['end_time'][()]
    
    lx_scales = f1.variables['lx_scales']
    lx_offsets = f1.variables['lx_offsets']
    fluid_volumes = f1.variables['fluid_volumes']
    
    if sys.version_info[0] < 3:
        vnames = f1.variables['variable_short_names']
        variable_names = deepcopy(vnames[0: -1]) # copy without last 0
        variable_names[variable_names == ''] = ' ' # replace 0 by space
        var_list = ''.join(variable_names).split(' ') # join and split
    else:
        x = f1.variables['variable_short_names'][()]
        var_list = x.tostring().decode('utf-8').split('\x00')[:-1]
    
    print(var_list)
    
    # Parse var name to find index
    ipress=-1
    irho=-1
    ivx=-1
    for iv,var in enumerate(var_list):
        if var=='x_velocity':
            ivx=iv
        elif var=='y_velocity':
            ivy=iv
        elif var=='z_velocity':
            ivz=iv
        elif var=='static_pressure':
            ipress=iv
        elif var=='density':
            irho=iv
    var_list=[]
    var_idx=[]
    if ipress==-1:
        var_list.append('density')
        var_idx.append(irho)
    else:
        var_list.append('static_pressure')
        var_idx.append(ipress)
    if ivx>=0:
        var_list.append('x_velocity')
        var_list.append('y_velocity')
        var_list.append('z_velocity')
        var_idx.append(ivx)
        var_idx.append(ivy)
        var_idx.append(ivz)
    
    print("Scaling factors")
    coeff_dx = lx_scales[0]
    CFL_number = lx_scales[-1]
    coeff_dt = lx_scales[0] / lx_scales[3]
    dt = CFL_number * coeff_dt
    coeff_press = lx_scales[1] * lx_scales[3]**2
    coeff_vel = lx_scales[3]
    coeff_density = lx_scales[1]
    offset_pressure = lx_offsets[4]
    
    print("weight factor")
    t_lat = 1.0/3.0 # lattive temperature
    r_lat = 1.0 # lattice r_gas constant
    rTlat = r_lat*t_lat # weight for density to pressure
    irTlat = 1.0/rTlat # weight for pressure to density
    
    # Average point
    ivol=1.0/float(np.sum(fluid_volumes[:]))
    
    mean_meas=np.sum(meas1[()]*fluid_volumes[:],axis=2)*ivol
    
    data = OrderedDict()
    
    data['time'] = 0.5*(st1+et1)*dt
    
    if ipress==-1:
        data[var_list[0]]=mean_meas[:,var_idx[0]] * coeff_density
        data['static_pressure']=(mean_meas[:,var_idx[0]] * rTlat + offset_pressure) * coeff_press
    else:
        data[var_list[0]]=(mean_meas[:,var_idx[0]]+offset_pressure)*coeff_press
        data['density']=mean_meas[:,var_idx[0]] * irTlat * coeff_density
    
    if ivx>=0:
        data[var_list[1]]=mean_meas[:,var_idx[1]] * coeff_vel
        data[var_list[2]]=mean_meas[:,var_idx[2]] * coeff_vel
        data[var_list[3]]=mean_meas[:,var_idx[3]] * coeff_vel
    
    f1.close()
    
    header = ' '.join(data.keys())
    nf = len(data.keys())
    nt = data['time'].size
    X = np.zeros((nt,nf))
    for iff,fn in enumerate(data.keys()):
        X[:,iff] = data[fn]
    
    
    np.savetxt(os.path.join(directory, '{0:s}.txt'.format(probe_file.replace('.pfnc',''))),X,header=header, comments='#')
    
    return(os.path.join(directory, '{0:s}.txt'.format(probe_file.replace('.pfnc',''))))
    
    
#read_probe_file('Medium_Case/Case_26-02-2023/Post_processing/Data/','mic7.pfnc')
