import pdb
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from read_files import *
from post_processing import *
from pf_visualization import *

images_path = '/Users/jmrendona/Library/CloudStorage/OneDrive-USherbrooke/PhD/EDAT/Simulations/Iteration1/baseline/images'
exp_data_path = '/Users/jmrendona/Library/CloudStorage/OneDrive-USherbrooke/PhD/EDAT/Experiments/Bell'
sim_data_path = '/Users/jmrendona/Library/CloudStorage/OneDrive-USherbrooke/PhD/EDAT/Simulations/Iteration1/baseline/data'
sim_data_path_t = '/Users/jmrendona/Library/CloudStorage/OneDrive-USherbrooke/PhD/rotor-alone/6e-5-6000rpm-transition/data'
sim_data_comp = '/Users/jmrendona/Library/CloudStorage/OneDrive-USherbrooke/PhD/rotor-alone/comparison/data'
#read = int(input('If the .pfnc files want to be converted enter 1, otherwise enter a 0\n'))
# mics_sim = [0,90]
mics_sim = np.arange(105,180,5)
mics_exp = [16,9]
# angles = np.arange(0,180,5)
angles = [0]
mics=['000001','000013','000029']
geometry = 'casing'
pitch_key1 = 103
pitch_key2 = 109
speed = 3500
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
# pfncconvert(1,mics_sim,angles,sim_data_path)
# time_trace_pfnc(speed,angles,mics_sim,sim_data_path,images_path,'Plot')
# time_trace_psnc(speed,duct_pr,duct_z,'velocity',sim_data_path,images_path,'Plot')
# time_trace_psnc(speed,duct_pr,duct_z,'vorticity',sim_data_path,images_path,'Plot')
# time_trace_psnc(speed,duct_pr,duct_z,'density',sim_data_path,images_path,'Plot')
# time_trace_psnc(speed,duct_pr,duct_z,'pressure',sim_data_path,images_path,'Plot')
# single_welch(n_blades,speed,angles,mics_sim,sim_data_path,images_path,'Plot',0.08,'density')
# psd_comparison(n_blades,speed,pitch_key1,angles,mics_sim,mics_exp,sim_data_path,exp_data_path,images_path,'Plot',0.1,'density')
# welch_all(n_blades,speed,angles,mics_sim,pitch_key1,mics_exp,exp_data_path,sim_data_path,0.1,'density','rotor-full','stator-V2')
# simcomparisonpsd(n_blades,speed,angles,mics_sim,mics_sim,sim1_data_path,sim2_data_path,images_path1,'Plot',0.1,'density')
# psd_psnc(4,speed,'velocity',duct_pr,duct_z,sim_data_path,images_path,'Plot')
# psd_psnc(4,speed,'vorticity',duct_pr,duct_z,sim_data_path,images_path,'Plot')
# psd_psnc(4,speed,'density',duct_pr,duct_z,sim_data_path,images_path,'Plot')
# psd_psnc(4,speed,'pressure',duct_pr,duct_z,sim_data_path,images_path,'Plot')
# wpf(n_blades,speed,probes,sim2_data_path,images_path2,'Plot',0.1,'density')
# coherence(n_blades,speed,images_path1,'mic1-pr3',sim1_dataq_path,'mic1_0.txt',sim1_data_path,'duct-pr-4.txt')
# sherfwhwelch(n_blades,speed,mics,sim1_data_path,images_path1,geometry,'Plot','density')
# time,fy_data,filter_time_m,filter_fy_m =
# forces_calculation(speed,'M','z','integrated',images_path,sim_data_path,'CMF-forces-all.txt','Plot',transient_rev=10,n_rev=1)
# time,fx_data,filter_time_m,filter_fx_m =
# forces_calculation(speed,'M','z','integrated',images_path,sim_data_path,'CMF-forces-casing.txt','Plot',transient_rev=10,n_rev=1)
# time,fz_data,filter_time_m,filter_fz_m = 
# time,force_data,filter_time,filter_force = forces_calculation(speed,'F','y','integrated',sim_data_path,'CMF_rotor.txt','Plot',transient_rev=9,n_rev=1)
# time_t,force_data_t,filter_time_t,filter_force_t = forces_calculation(speed,'F','y','integrated',sim_data_path_t,'CMF_rotor.txt','Plot',transient_rev=9,n_rev=1)
# forces_calculation(speed,'M','y','integrated',sim_data_path,'CMF_rotor.txt','Plot',transient_rev=9,n_rev=1)
# forces_calculation(speed,'F','z','integrated',sim_data_path,'CMF-forces-duct.txt','Plot',transient_rev=10,n_rev=1)
# forces_calculation(speed,'M','z','integrated',sim_data_path,'CMF-forces-duct.txt','Plot',transient_rev=10,n_rev=1)
# f,t,filter_time_m,filter_fz_m=forces_calculation(speed,'M','z','integrated',sim_data_path,'CMF-forces-single-stator.txt','Plot',transient_rev=10,n_rev=1)
# forces_calculation(speed,'F','z','integrated',sim_data_path,'CMF-forces-stator.txt','Plot',transient_rev=10,n_rev=1)
# time,force_data,filter_time_f,filter_force_f = forces_calculation(6420,'F','y',images_path,sim_data_path,'Forces.txt','Plot',transient_rev=10,n_rev=4)
# induced_velocity_fluidline(0.0229,0.089,20,'FluidPoint-LineGraph-12mm.csv',sim_data_path,images_path)
# psnc_surface_avg(3500,duct_pr,angles,'velocity',sim_data_path,images_path1,'Plot')
# psnc_surface_avg(3500,duct_pr,angles,'pressure',sim_data_path,images_path1,'Plot')
# psnc_surface_avg(3500,duct_pr,angles,'velocity',sim_data_path,images_path1,'Plot',n_rev=1)
# psnc_surface_avg(3500,duct_pr,angles,'pressure',sim_data_path,images_path1,'Plot',n_rev=1)
# read_mass_flux(sim_data_path,'MassFlowRate_rear.csnc')
mass_flux(3500,images_path,sim_data_path,'MassFlowRate_rear.txt','Plot')
# read_mass_flux(sim_data_path,'MassFlowRate_front.csnc')
mass_flux(3500,images_path,sim_data_path,'MassFlowRate_front.txt','Plot')
# spectrogram(2,6420,mics_sim,angles,sim_data_path,images_path,'fluid',8)
# spectrogram(4,3500,duct_pr,duct_z,sim_data_path,images_path,5)
# polar_forces_plot(speed,sim_data_path,'CMF_rotor.txt','Plot',transient_rev=9,n_rev=1)
# polar_forces_plot(speed,sim_data_path,'CMF-forces-casing.txt','Plot',transient_rev=10,n_rev=2)
# polar_forces_plot(speed,sim_data_path,'CMF-forces-rotor.txt','Plot',transient_rev=10,n_rev=2)
# polar_forces_plot(speed,sim_data_path,'CMF-forces-single-stator.txt','Plot',transient_rev=10,n_rev=4)
# polar_forces_plot(speed,sim_data_path,'CMF-forces-blades.txt','Plot',transient_rev=10,n_rev=2)
# polar_forces_plot(speed,sim_data_path,'CMF-forces-duct.txt','Plot',transient_rev=10,n_rev=2)
# polar_forces_plot(speed,sim_data_path,'CMF-forces-single-blade.txt','Plot',transient_rev=10,n_rev=4)
# polar_forces_plot(speed,sim_data_path,'CMF-forces-stator.txt','Plot',transient_rev=10,n_rev=2)

# pdb.set_trace()


# --------------------------------------- Visualization --------------------------------------#

# pf_fnc_visualization(6000,0.25121,5.03808,5.34529,[-0.4,0.4],[-1.1,-0.5],sim_data_path,'avg-vor-mag-zyplane.txt','avg-vor-mag-zyplane_t.txt')
# pf_delta_fnc_visualization(6000,0.25121,5.03808,5.34529,[-0.4,0.4],[-1.1,-0.5],sim_data_path,'avg-vel-mag-zyplane.txt','avg-vel-mag-zyplane_t.txt')
# pf_snc_cp(1.204,101325,3500,0.018,0.122,10,sim_data_comp)
# pf_snc_cp_comparison(1.204,101325,6000,0.018,0.122,10,sim_data_comp,'radius','cp_r0p122.txt','cp_t_r0p122.txt')

print('Time to do magic')

# pdb.set_trace()

# plt.plot(filter_time_m,filter_fz_m,'k',label='$F_y$')
# # plt.plot(filter_time_m,filter_fz_m,'r',label='$F_z$')
# # plt.plot(filter_time_m,filter_fx_m,'b',label='$F_x$')
# for i,j in zip(range(6),range(6)):
#     plt.axvline(0.17321 + (i * 0.00286),color='red',linestyle='--')
#     # plt.axvline(0.173 + 0.00286/2 + (j*0.00286),color='blue',linestyle='--')
# plt.grid()
# plt.xlabel('Time $[s]$',fontsize=14, style='italic')
# plt.ylabel('Force $[N]$',fontsize=14, style='italic')
# # plt.title('Filtered Thrust of the system')
# plt.xticks(fontsize=14)
# plt.yticks(fontsize=14)
# ax=plt.gca()
# ax.tick_params(axis='y')
# ax.set_yticklabels([])
# ax.set_ylabel('Thrust $[N]$',fontsize=14, style='italic')
# ax2 = ax.twinx()  # Create a second axes sharing the same y-axis
# ax2.plot(filter_time_m,np.abs(filter_fz_m),'b',label='$F_z$')
# ax2.set_ylabel('$F_z/F_x [N]$',fontsize=14, style='italic')
# ax2.tick_params(axis='y')
# # ax3 = ax.twinx()
# ax2.plot(filter_time_m,np.abs(filter_fx_m),'r',label='$F_x$')
# # ax3.set_ylabel('$F_x [N]$', color='red',fontsize=14, style='italic')
# # ax3.tick_params(axis='y', labelcolor='red')
# # ax3.set_yticklabels([])
# ax2.yaxis.tick_right()
# ax3.yaxis.tick_right()
# ax3.spines['right'].set_position(('outward', 50))
# plt.legend()
# plt.tight_layout()
# plt.savefig(os.path.join(images_path, 'forces/filtered', 'All-forces-filtered-4rev-2cale.png'), dpi=600)
# plt.show()
# plt.close()
