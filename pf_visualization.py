import numpy as np 
import matplotlib.pyplot as plt
import pandas as pd
import h5py
import os
import re
import pdb

def pf_fnc_visualization():
    
    data=np.genfromtxt('avg-vel-mag-zyplane.txt',skip_header=15,filling_values=0)
    ny, nx = data.shape
    lx = 5.03808
    ly = 5.34529
    x = np.linspace(0,lx,nx)
    y = np.linspace(0,ly,ny)
    X, Y = np.meshgrid(x, y)
    Z = data
    cp=plt.contourf(X,Y,Z)
    plt.colorbar(cp)
    plt.xlim([2,3])
    plt.ylim([1,3])
    plt.show()
    
def pf_snc_cp(rho:float,char_p:float,rpm:int,lower_r:float,upper_r:float,r_number:int,path:str,variable:str):
    
    radius = np.append(np.linspace(upper_r,lower_r,r_number),np.linspace(-lower_r,-upper_r,r_number))
    os.makedirs(os.path.join(os.path.dirname(path),f'images/{variable}'), exist_ok=True)

    for r in radius:
        
        name = str(r).split('.')[0] + 'p' + str(np.round(r,decimals=3)).split('.')[1]
        
        with open(os.path.join(path,'cp_t_r'+ name +'.txt')) as f:
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

        if r >= 0:
            plt.plot(-suction_s_x/(np.max(pressure_s_x)-np.min(pressure_s_x)),-cp_suction,'k')
            plt.plot(-pressure_s_x/(np.max(pressure_s_x)-np.min(pressure_s_x)),-cp_pressure,'k')
            
        else:
            plt.plot(suction_s_x/(np.max(suction_s_x)-np.min(suction_s_x)),-cp_suction,'k')
            plt.plot(pressure_s_x/(np.max(pressure_s_x)-np.min(pressure_s_x)),-cp_pressure,'k')
        
        plt.grid(True, which='both', ls='--')
        plt.xlabel('$x/C$', fontsize=16)
        plt.ylabel('$-C_p$', fontsize=16)
        plt.ylim([-1.2,1.2])
        plt.xlim([-0.6,0.6])
        plt.xticks(fontsize=12)
        plt.yticks(fontsize=12)
        print(f'Saving Cp for radius {r}\n')
        print(40*'-')
        plt.tight_layout()
        plt.savefig(os.path.join(os.path.dirname(path),f'images/{variable}','cp_t_r'+ name +'.png'))
        # plt.show()
        
        # plt.show()
        plt.close()