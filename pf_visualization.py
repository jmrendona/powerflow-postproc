import numpy as np 
import matplotlib.pyplot as plt
import h5py
import os

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