import numpy as np 
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pandas as pd
import h5py
import os
import re
import pdb
from scipy.signal import savgol_filter

def pf_fnc_visualization(speed:float,diameter:float,lx,ly,xlim:list,ylim:list,path:str,*files:str):
    
    u_tip = 1#(np.pi*(diameter/2)*speed)/60
    levels = np.linspace(0, 40, 101)
    
    for file in files:
        print(f'Reading file: {file}')
        data=np.genfromtxt(os.path.join(path,file),skip_header=15,filling_values=0)
        ny, nx = data.shape
        lx = 5.03808
        ly = 5.34529
        x = np.linspace(-lx/2,lx/2,nx)
        y = np.linspace(-ly/2,ly/2,ny)
        X, Y = np.meshgrid(x, y)
        Z = data
        print('Plotting')
        cp=plt.contourf(X,Y,Z/u_tip,levels=levels,cmap='turbo',vmin=0,vmax=30)
        cbar = plt.colorbar(cp)
        cbar.set_label('$|\omega|/RPS~[-]$', fontsize=14)
        cbar.ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.2f'))
        plt.xlim(xlim)
        plt.xticks(fontsize=12)
        plt.ylim(ylim)
        plt.yticks(fontsize=12)
        plt.tight_layout()
        print('Saving plot')
        plt.savefig(os.path.join(os.path.dirname(path),'images',file.split('.')[0] +'.png'), dpi=600)
        plt.show()
        

def pf_delta_fnc_visualization(speed:float,diameter:float,lx,ly,xlim:list,ylim:list,path:str,ref_file:str,*comp_files:list):
    
    u_tip = 1#(np.pi*(diameter/2)*speed)/60
    levels = np.linspace(0, 0.02, 101)
    
    for file in comp_files:
        print(f'Reading file: {ref_file}')
        ref_data=np.genfromtxt(os.path.join(path,ref_file),skip_header=15,filling_values=0)
        print(f'Reading file: {file}')
        comp_data=np.genfromtxt(os.path.join(path,file),skip_header=15,filling_values=0)
        ny, nx = ref_data.shape
        lx = 5.03808
        ly = 5.34529
        x = np.linspace(-lx/2,lx/2,nx)
        y = np.linspace(-ly/2,ly/2,ny)
        X, Y = np.meshgrid(x, y)
        print('Computing difference')
        Z = np.abs(ref_data/u_tip - comp_data/u_tip)
        print('Plotting')
        cp=plt.contourf(X,Y,Z,levels=levels,cmap='turbo_r')#,vmin=-5,vmax=5)
        cbar = plt.colorbar(cp)
        cbar.set_label('$U/U_{tip}~[-]$', fontsize=14)
        cbar.ax.yaxis.set_major_formatter(ticker.FormatStrFormatter('%.3f'))
        plt.xlim(xlim)
        plt.xticks(fontsize=12)
        plt.ylim(ylim)
        plt.yticks(fontsize=12)
        print('Saving plot')
        plt.tight_layout()
        plt.savefig(os.path.join(os.path.dirname(path),'images','delta_'+ file.split('.')[0] +'.png'), dpi=600)
        plt.show() 
    
def pf_snc_cp(rho:float,char_p:float,rpm:int,lower_r:float,upper_r:float,r_number:int,path:str):
    
    radius = np.append(np.linspace(upper_r,lower_r,r_number),np.linspace(-lower_r,-upper_r,r_number))
    radius = [0.367,0.348,0.309,0.194,0.116]
    os.makedirs(os.path.join(os.path.dirname(path),f'images/cp'), exist_ok=True)

    for r in radius:
        
        name = str(r).split('.')[0] + 'p' + str(np.round(r,decimals=3)).split('.')[1]
        
        with open(os.path.join(path,'cp_r'+ name +'.txt')) as f:
            lines = f.readlines()
            
        start_idx = next(i for i, line in enumerate(lines) if re.match(r'^[\s\-0-9\.]', line))
        content = "".join(lines[start_idx:])
        
        suction_s = content.split('\ne\n\n')[0]
        pressure_s = content.split('\ne\n\n')[1]
        
        lines_ss = [line for line in suction_s.splitlines() if re.match(r'^[\s\-0-9\.]', line)]
        lines_ps = [line for line in pressure_s.splitlines() if re.match(r'^[\s\-0-9\.]', line)]
        
        suction_s_p = np.loadtxt(lines_ss)[:,1]
        suction_s_x = np.loadtxt(lines_ss)[:,0]
        pressure_s_p = np.loadtxt(lines_ps)[:,1]
        pressure_s_x = np.loadtxt(lines_ps)[:,0]
        
        u_inf = (2*np.pi*rpm*np.abs(r))/60
        cp_suction = (suction_s_p - char_p)/(0.5 * rho * u_inf ** 2)
        cp_pressure = (pressure_s_p - char_p)/(0.5 * rho * u_inf ** 2)
        
        cp_s_smooth = savgol_filter(cp_suction, window_length=51, polyorder=3)
        cp_p_smooth = savgol_filter(cp_pressure, window_length=51, polyorder=3)

        if r >= 0:
            plt.plot(-suction_s_x/(0.25121/2),-cp_s_smooth,'k')
            plt.plot(-pressure_s_x/(0.25121/2),-cp_p_smooth,'k')
            
        else:
            plt.plot(suction_s_x/(0.25121/2),-cp_s_smooth,'k')
            plt.plot(pressure_s_x/(0.25121/2),-cp_p_smooth,'k')
        
        plt.grid(True, which='both', ls='--')
        plt.xlabel('$x/C$', fontsize=16)
        plt.ylabel('$-C_p$', fontsize=16)
        # plt.ylim([-1.1,1.5])
        # plt.xlim([-0.12,0.12])
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        print(f'Saving Cp for radius {r}\n')
        print(40*'-')
        plt.tight_layout()
        plt.savefig(os.path.join(os.path.dirname(path),f'images/cp','cp_r'+ name +'.png'))
        # plt.show()
        plt.close()
        
        
        
def pf_snc_cp_comparison(rho:float,char_p:float,rpm:int,lower_r:float,upper_r:float,r_number:int,path:str,iteration:str,*files:str):
    
    radius = np.append(np.linspace(upper_r,lower_r,r_number),np.linspace(-lower_r,-upper_r,r_number))
    radius = [0.122,0.099,0.076,0.053,0.03]
    labels = [95,80,60,40,25]
    os.makedirs(os.path.join(os.path.dirname(path),f'images/cp'), exist_ok=True)

    if iteration == 'radius':
        first_it = files
        second_it = radius
    elif iteration == 'case':
        first_it = radius
        second_it = files

    
    for i in first_it:
        
        k = 0
        cmap = plt.get_cmap('turbo')  
        n_colors = len(second_it)
        colors = [cmap(i / (n_colors - 1)) for i in range(n_colors)]
        
        for r in second_it:
            
            name = str(i).split('r')[0] + 'r0p' + str(np.round(r,decimals=3)).split('.')[1]

            with open(os.path.join(path, name +'.txt')) as f:
                lines = f.readlines()
                
            start_idx = next(i for i, line in enumerate(lines) if re.match(r'^[\s\-0-9\.]', line))
            content = "".join(lines[start_idx:])
            
            suction_s = content.split('\ne\n\n')[0]
            pressure_s = content.split('\ne\n\n')[1]
            
            lines_ss = [line for line in suction_s.splitlines() if re.match(r'^[\s\-0-9\.]', line)]
            lines_ps = [line for line in pressure_s.splitlines() if re.match(r'^[\s\-0-9\.]', line)]
            
            suction_s_p = np.loadtxt(lines_ss)[:,1]
            suction_s_x = np.loadtxt(lines_ss)[:,0]
            pressure_s_p = np.loadtxt(lines_ps)[:,1]
            pressure_s_x = np.loadtxt(lines_ps)[:,0]
            
            u_inf = (2*np.pi*rpm*np.abs(r))/60
            cp_suction = (suction_s_p - char_p)/(0.5 * rho * u_inf ** 2)
            cp_pressure = (pressure_s_p - char_p)/(0.5 * rho * u_inf ** 2)
            
            cp_s_smooth = savgol_filter(cp_suction, window_length=51, polyorder=3)
            cp_p_smooth = savgol_filter(cp_pressure, window_length=51, polyorder=3)
            
            color = colors[k]
            
            if r >= 0:
                plt.plot(-suction_s_x/(0.025),-cp_s_smooth, color=color, label=f'$r={labels[k]}~\%$')
                plt.plot(-pressure_s_x/(0.025),-cp_p_smooth, color=color)
            
            else:
                plt.plot(-suction_s_x/(0.025),-cp_s_smooth, color=color, label=f'$r={labels[k]}~\%$')
                plt.plot(-pressure_s_x/(0.025),-cp_p_smooth, color=color)
                
            k += 1
            
        plt.grid(True, which='both', ls='--')
        plt.xlabel('$x/C$', fontsize=16)
        plt.ylabel('$-C_p$', fontsize=16)
        #plt.ylim([-1.1,1.5])
        #plt.xlim([-0.12,0.12])
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        plt.legend(fontsize=10)
        print(f'Saving Cp for radius {r}\n')
        print(40*'-')
        plt.tight_layout()
        plt.savefig(os.path.join(os.path.dirname(path),f'images/cp','comp'+ name +'.png'))
        # plt.show()
        plt.close()