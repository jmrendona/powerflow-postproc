import pdb
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from read_files import *
from post_processing import *
from pf_visualization import *

# ------------------------- Inputs to be changed -------------------------- #

images_path = '/Users/jmrendona/Library/CloudStorage/OneDrive-USherbrooke/PhD/EDAT/Simulations/Iteration1/baseline/images'
exp_data_path = '/Users/jmrendona/Library/CloudStorage/OneDrive-USherbrooke/PhD/EDAT/Experiments/Bell'
sim_data_path = '/home/renj3003/PhD/rotor-alone/UdeS_Case/6e-5_6000rpm/data'
sim_data_path_t = '/Users/jmrendona/Library/CloudStorage/OneDrive-USherbrooke/PhD/rotor-alone/6e-5-6000rpm-transition/data'
sim_data_comp = '/Users/jmrendona/Library/CloudStorage/OneDrive-USherbrooke/PhD/rotor-alone/comparison/data'
#read = int(input('If the .pfnc files want to be converted enter 1, otherwise enter a 0\n'))
mics_sim = [0,60,90,120]
#mics_sim = np.arange(0,180,10)
mics_exp = [16,9]
angles = np.arange(0,360,10)
# angles = [0]
mics=['000001','000013','000029']
geometry = 'casing'
pitch_key1 = 103
pitch_key2 = 109
speed = 6000
n_blades = 2
probes = [4,6,7,8]
duct_pr = np.arange(40,96,5)
# duct_pr = [0]
single_pr = [36]
vanes_pr = np.arange(40,101,5)
single_vanes = [40,45]
duct_z = [36,38,41,44,46,48,51,54,56,57,61]

# pdb.set_trace()

# experimental_check('/run/media/renj3003/BckSmoreau6/Jovan-McGill','95phi-4.h5','Ds01-Time','Ds02-Signal 13','95phi-4-300k')
# experimental_check('/run/media/renj3003/BckSmoreau6/Jovan-McGill','0phi-11ms-3.h5','Ds001-Time','Ds084-Signal 13','0phi-11ms-3-900k')
# experimental_check('/run/media/renj3003/BckSmoreau6/Jovan-McGill','0phi-3.h5','Ds01-Time','Ds14-Signal 13','0phi-3-400k')

# ------------------------- Conversion of PF mics ------------------------- #

# pfnc_convert(1,mics_sim,angles,sim_data_path)

# ------------------------- Timetrace of PF mics -------------------------- #

# time_trace_pfnc(speed,mics_sim,angles,sim_data_path,'Plot')

# ------------------------ Timetrace of PF probes ------------------------- #

# time_trace_psnc(speed,duct_pr,duct_z,'velocity',sim_data_path,images_path,'Plot')
# time_trace_psnc(speed,duct_pr,duct_z,'vorticity',sim_data_path,images_path,'Plot')
# time_trace_psnc(speed,duct_pr,duct_z,'density',sim_data_path,images_path,'Plot')
# time_trace_psnc(speed,duct_pr,duct_z,'pressure',sim_data_path,images_path,'Plot')

# ------------------------- PSD of single PF mics ------------------------- #

# single_welch(n_blades,speed,angles,mics_sim,sim_data_path,images_path,'Plot',0.04,'density')

# ---------------- PSD comparison of multiple sources mics ---------------- #

# psd_comparison(n_blades,speed,pitch_key1,angles,mics_sim,mics_exp,sim_data_path,exp_data_path,images_path,'Plot',0.1,'density')
# welch_all(n_blades,speed,angles,mics_sim,pitch_key1,mics_exp,exp_data_path,sim_data_path,0.1,'density','rotor-full','stator-V2')
# simcomparisonpsd(n_blades,speed,angles,mics_sim,mics_sim,sim1_data_path,sim2_data_path,images_path1,'Plot',0.1,'density')

# ------------------------ PSD of PF probes (psnc) ------------------------ #

# psd_psnc(4,speed,'velocity',duct_pr,duct_z,sim_data_path,images_path,'Plot')
# psd_psnc(4,speed,'vorticity',duct_pr,duct_z,sim_data_path,images_path,'Plot')
# psd_psnc(4,speed,'density',duct_pr,duct_z,sim_data_path,images_path,'Plot')
# psd_psnc(4,speed,'pressure',duct_pr,duct_z,sim_data_path,images_path,'Plot')
# wpf(n_blades,speed,probes,sim2_data_path,images_path2,'Plot',0.1,'density')

# -------------------- Spectrogram of PF mics & probes -------------------- #

# spectrogram(n_blades,speed,mics_sim,angles,sim_data_path,images_path,'fluid',1)


# coherence(n_blades,speed,images_path1,'mic1-pr3',sim1_dataq_path,'mic1_0.txt',sim1_data_path,'duct-pr-4.txt')
# sherfwhwelch(n_blades,speed,mics,sim1_data_path,images_path1,geometry,'Plot','density')

# -------------------- Forces calculation from *.csnc  --------------------- #

# forces_calculation(speed,'M','y','integrated',sim_data_path,'CMF_rotor.txt','Plot',transient_rev=10,n_rev=1)
# forces_calculation(speed,'F','y','integrated',sim_data_path,'CMF_rotor.txt','Plot',transient_rev=10,n_rev=1)

# ------------------ Polar plot calculation from *.csnc  ------------------- #

polar_forces_plot(speed,sim_data_path,'CMF_rotor.txt','Plot',transient_rev=9,n_rev=1)


# induced_velocity_fluidline(0.0229,0.089,20,'FluidPoint-LineGraph-12mm.csv',sim_data_path,images_path)
# psnc_surface_avg(3500,duct_pr,angles,'velocity',sim_data_path,images_path1,'Plot')
# psnc_surface_avg(3500,duct_pr,angles,'pressure',sim_data_path,images_path1,'Plot')
# psnc_surface_avg(3500,duct_pr,angles,'velocity',sim_data_path,images_path1,'Plot',n_rev=1)
# psnc_surface_avg(3500,duct_pr,angles,'pressure',sim_data_path,images_path1,'Plot',n_rev=1)

# -------------------- Mass flow verification of *.csnc ------------------- #

# read_mass_flux(sim_data_path,'CSFM_mass_rotor_inlet.csnc')
# mass_flux(speed,sim_data_path,'CSFM_mass_rotor_inlet.txt','Plot')
# read_mass_flux(sim_data_path,'CSFM_mass_rotor_outlet.csnc')
# mass_flux(speed,sim_data_path,'CSFM_mass_rotor_outlet.txt','Plot')



# pdb.set_trace()


# --------------------------------------- Visualization --------------------------------------#

# pf_fnc_visualization(6000,0.25121,5.03808,5.34529,[-0.4,0.4],[-1.1,-0.5],sim_data_path,'avg-vor-mag-zyplane.txt','avg-vor-mag-zyplane_t.txt')
# pf_delta_fnc_visualization(6000,0.25121,5.03808,5.34529,[-0.4,0.4],[-1.1,-0.5],sim_data_path,'avg-vel-mag-zyplane.txt','avg-vel-mag-zyplane_t.txt')
# pf_snc_cp(1.204,101325,3500,0.018,0.122,10,sim_data_comp)
# pf_snc_cp_comparison(1.204,101325,6000,0.018,0.122,10,sim_data_comp,'radius','cp_r0p122.txt','cp_t_r0p122.txt')

print('Time to do magic')

# pdb.set_trace()
