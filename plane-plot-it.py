import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from scipy.interpolate import griddata
from scipy.spatial import cKDTree
import pdb
from math_operations import *
import os
from vtk_export import VTKExporter

variable = 'vel_mag'
i_theta = 0
n_theta = 180
delta_theta = 10
norm_value = 58.33
max_value = 6000
y_label = '$\omega/RPS~[-]$'
c_map = 'turbo'
degrees_f = np.arange(i_theta+delta_theta, n_theta + 1, delta_theta)
levels = np.linspace(0, 1000, 101)#max_value/norm_value, 101)
#levels = 100

for degree_i in degrees_f:
    run = 't'+str(int(degree_i/delta_theta))
    degrees = np.arange(degree_i, n_theta + 1, delta_theta) 
    degrees = np.append(degrees,np.arange(0,degree_i-delta_theta+1,delta_theta))
    #degrees = np.append(degrees,degrees)
    timesteps = np.arange(1,39)
    os.makedirs(f'images/planes/{variable}/{run}', exist_ok=True)
    os.makedirs(f'images/planes/{variable}/mean', exist_ok=True)

    data=np.genfromtxt(f'1t-radial-cut-0deg-{variable}.txt',skip_header=15)
    df=pd.DataFrame({'x':data[:,0],'y':data[:,1],'z':data[:,2],'var':data[:,3]})
    df2=df.groupby(['x','z'],as_index=False).mean()
    x_g=np.linspace(df2['x'].min(),df2['x'].max(),5000)
    z_g=np.linspace(df2['z'].min(),df2['z'].max(),5000)
    Xg,Yg=np.meshgrid(x_g,z_g)

    Z_all = []

    #with h5py.File(f'meridional-planes_{variable}.h5','w') as f:
        
        #f.create_dataset('X_coord', data=Xg)
        #f.create_dataset('Y_coord', data=Yg)
    for t,degree in zip(timesteps,degrees):

        print(10*'-',f'Computing plane at {degree} deg',10*'-')
        data=np.genfromtxt(f'{t}t-radial-cut-{degree}deg-{variable}.txt',skip_header=15)
        df=pd.DataFrame({'x':data[:,0],'y':data[:,1],'z':data[:,2],'var':data[:,3]})

        x_rot, y_rot, z_rot = rotate_plane(df['x'], df['y'], df['z'],
                                    axis='z',
                                    angle=-degree)
        print(10*'-','Rotation computed',10*'-')

        df = pd.DataFrame({'x':x_rot.round(6),'y':y_rot,'z':z_rot.round(6),'var':data[:,3]})
        df2=df.groupby(['x','z'],as_index=False).mean()

        print(10*'-','Interpolating Data',10*'-')
        Zlinear = griddata((df2['x'], df2['z']),df2['var'],(Xg,Yg), method='linear')
        Znearest = griddata((df2['x'], df2['z']),df2['var'],(Xg,Yg), method='nearest')
        Zcombined = np.where(np.isnan(Zlinear),Znearest,Zlinear)

        print(10*'-','Applying Mask',10*'-')
        pts = np.vstack([df['x'].values,df['z'].values]).T
        tree = cKDTree(pts)
        dist, idxs = tree.query(pts,k=2)
        nn_dist = dist[:,1]
        med = np.median(nn_dist)
        #threshold = 3*med
        threshold = np.percentile(nn_dist, 99.99)
        grid_pts = np.vstack([Xg.ravel(), Yg.ravel()]).T
        dist_grid, idx = tree.query(grid_pts, k=1)
        dist_grid = dist_grid.reshape(Xg.shape)

        mask = dist_grid <= 3*threshold

        Znearest_m = np.where(mask,Znearest,np.nan)
        Zcombined_m = np.where(mask,Zcombined,np.nan)
        #f.create_dataset(f'{variable}_{degree}deg', data=Zcombined_m)
        Z_all.append(Zcombined_m)
        
        print(10*'-', 'Generating Figures',10*'-')
        #plt.contourf(Xg,Yg,Znearest_m/norm_value,levels=levels,cmap=c_map)
        #cbar = plt.colorbar()
        #cbar.set_label(y_label, fontsize=14)
        #cbar.ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
        #plt.xlim([-0.4,0.4])
        #plt.xticks(fontsize=12)
        #plt.ylim([-0.2,0.05])
        #plt.yticks(fontsize=12)
        #print('Saving plot')
        #plt.tight_layout()
        #plt.savefig(os.path.join(f'images/planes/{variable}',f'Znearest-{degree}deg_{t}t_{variable}.png'),dpi=600)
        #plt.close()

        plt.figure(figsize=(10, 2.5))
        plt.contourf(Xg,Yg,Zcombined_m/norm_value,levels=levels,cmap=c_map)
        cbar = plt.colorbar()
        cbar.set_label(y_label, fontsize=14)
        cbar.ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
        plt.xlim([-0.4,0.4])
        plt.xticks(fontsize=12)
        plt.ylim([-0.2,0.05])
        plt.yticks(fontsize=12)
        print('Saving plot')
        plt.tight_layout()
        plt.savefig(os.path.join(f'images/planes/{variable}/{run}',f'Zcombined-{degree}deg_{t}{run}_{variable}.png'),dpi=600)
        plt.close()

        #plt.contourf(Xg,Yg,Znearest_m-Zcombined_m,levels=100,cmap=c_map)
        #cbar = plt.colorbar()
        #cbar.set_label('$\Delta U~[-]$', fontsize=14)
        #cbar.ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
        #plt.xlim([-0.4,0.4])
        #plt.xticks(fontsize=12)
        #plt.ylim([-0.2,0.05])
        #plt.yticks(fontsize=12)
        #print('Saving plot')
        #plt.tight_layout()
        #plt.savefig(os.path.join(f'images/planes/{variable}',f'Z-delta-{degree}deg_{t}t_{variable}.png'),dpi=600)
        #plt.close()

        #plt.contourf(Xg,Yg,np.abs(Znearest_m-Zcombined_m)*100/Znearest_m,levels=100,cmap=c_map)
        #cbar = plt.colorbar()
        #cbar.set_label('$\epsilon U~[\%]$', fontsize=14)
        #cbar.ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
        #plt.xlim([-0.4,0.4])
        #plt.xticks(fontsize=12)
        #plt.ylim([-0.2,0.05])
        #plt.yticks(fontsize=12)
        #print('Saving plot')
        #plt.tight_layout()
        #plt.savefig(os.path.join(f'images/planes/{variable}',f'Z-epsilon-{degree}deg_{t}t_{variable}.png'),dpi=600)
        #plt.close()
        print('\n')

    Z_stack = np.stack(Z_all, axis=0)
    Z_mean = np.nanmean(Z_stack, axis=0)
    #f.create_dataset(f'{variable}_mean', data=Z_mean)

    plt.figure(figsize=(10, 2.5))
    plt.contourf(Xg,Yg,Z_mean/norm_value,levels=levels,cmap=c_map)
    cbar = plt.colorbar()
    cbar.set_label(y_label, fontsize=14)
    cbar.ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
    plt.xlim([-0.4,0.4])
    plt.xticks(fontsize=12)
    plt.ylim([-0.2,0.05])
    plt.yticks(fontsize=12)
    print('Saving plot')
    plt.tight_layout()
    plt.savefig(os.path.join(f'images/planes/{variable}/mean',f'Zmean-{degree_i}deg_{variable}.png'),dpi=600)
    plt.close()

    Z_filled = np.nan_to_num(Z_mean, nan=-1e6)
    exporter = VTKExporter(Xg, Yg, Z_filled, variable)
    exporter.SaveVTI(f'{variable}_{degree_i}deg_mean-plane')
