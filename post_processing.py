import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as sg
from matplotlib import ticker
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import MaxNLocator
import pdb
import os
import pandas as pd
import itertools
import warnings
from read_files import *

# Here the functions for the different calculations will be done

plt.rc({'text.usetex':True})

def read_mic_file(full_filename,mic):
    """ Function to read microphone file (exported from B&K software)
    Inputs:
        full_filename: str
            filename to be read including relative or absolute path
        fmt : str
            hdf5 or uff string key
    Outputs:
        mat_data_NN : numpy array
            microphone_pressure of size (nmics,nlen)
        dt: float
            sampling time
    """

    mat_data_NN, dt, time = _h5_read_file(full_filename,mic)
    return mat_data_NN, dt, time

def data_filtering(data1,data2,t_int,t_end):
    
    start_time_i = next(x for x, val in enumerate(data1) if val > t_int)
    end_time_i = next(y for y, val in enumerate(data1) if val > t_end)
    data1 = data1[start_time_i:end_time_i]
    data2 = data2[start_time_i:end_time_i]

    return data1, data2, [start_time_i,end_time_i]

def variable_selection(variable,data_full):

    if variable == 'pressure':
        data = data_full[1] - np.mean(data_full[1])
        y_label = 'Pressure $[Pa]$'

    elif variable == 'velocity':
        data = data_full[2]
        y_label = 'Velocity Magnitude $[m/s]$'

    elif variable == 'vorticity':
        data = data_full[3]
        y_label = 'Vorticity Magnitude $[rad/s]$'

    elif variable == 'density':
        data = data_full[4]
        y_label = 'Density $[kg/m^3]$'

    else:
        print('Not a valid variable or not yet implemented')

    return(data,y_label)

def _h5_read_file(full_filename,mic):
    import h5py

    # Fake initialization
    dt = 0.
    mat_data_NN = []

    # Open file in read only mode
    fid = h5py.File(full_filename, 'r')

    # Get measured signals at microphone locations
    # Assumption they are in collector order
    signals = []
    fields = list(fid['Table1'].keys())

    #pdb.set_trace()

    for field in fields:
        if 'Time' in field:
            # Sampling period
            dt = fid['Table1'][field][1]
            # Signal length
            nlen = fid['Table1'][field].size
            time = np.empty((1,nlen))
            time = fid['Table1'][field][()]
        elif 'Signal' in field:
            signals.append(field)

    # # Number of microphone
    # nmics = len(signals)

    # # Collect timetrace signals in a matrix
    # mat_data_NN = np.empty((nmics,nlen))
    # for imic,field in enumerate(signals):
    #     pdb.set_trace()
    #     mat_data_NN[imic,:] = fid['Table1'][field][()]

    # For only one mic
    nmics = 1

    # Collect timetrace signals in a matrix
    mat_data_NN = np.empty((nmics,nlen))
    mat_data_NN[0,:] = fid['Table1'][fields[mic]][()]


    # Close file
    fid.close()



    return mat_data_NN,dt,time

def pfncconvert(read,mics,angles,data_path):

    """
    This function helps to convert in a loop various .pfnc files from PowerFlow\n
    to .txt files for a further analysis.\n
    \n
    The needed inputs are the followings:\n
    - read = given input from the user. 1 to perform conversion, 0 otherwise\n
    - mics = an array with the number of mics to convert (file name)\n
    - angles = an array with the number of angles to convert (file name)\n
    - data_path = the path where the .pfnc files are located \n

    """

    if read == 1:

        print('The conversion of the .pfnc files into .txt files will start\n')

        for i in angles:
            for j in mics:
                print('--------------------------------\n')
                print(f'Converting mic {j} at angle {i}\n')
                globals()[f'mic_{i}_{j}'] = read_probe_file(data_path ,f'mic_{i}_{j}.pfnc')

    else:

        print('No conversion will be done. Check the given input')

def time_trace_pfnc(rpm:int,first_it:list,second_it:list,data_path:str,image_path:str,case:str,transient_rev=10,n_rev=2):

    '''
    This function implements the ploting and analysis of the time trace of fluid probes\n
    converted with the in-house tool.\n
    A complete plot of the signal is done togehter with a close up to a user defined number\n
    of revolutions.\n
    To use it the base filename needs to be adjusted to match with the index you use for\n
    your specific application (2 index name supported)\n
    \n
    - rpm = revolutions per minute of the rotatin machine used for the close up image\n
    - first_it = first variable to iterate over shown in the first promt of the script\n
    - second_it = second variable to iterate over shown in the first promt of the script\n
    - data_path = global path where the data is located\n
    - image_path = global path where the images will be stored\n
    - case = either Plot or Data depending if the graph want to be saved or not\n
    - transient_rev = revolution where the close up will start. Prior knowledge of the signal is required\n
    - n_rev = revolutions to consider for the close up. Prior knowledge of the signal is required\n
    '''

    base_filename = 'mic_0_0.txt' # To compare with input line 121
    print('The base file name is as:', base_filename)
    rev_time = 60/rpm

    for i in first_it:
        for j in second_it:

            file_name = f'mic_{i}_{j}.txt' # Input to change
            image_name = file_name.split('.')[0]

            data = pd.read_csv(os.path.join(data_path, file_name),delimiter='\s+')
            figure_extension = '.png'
            time = data['#time']
            pressure = data['static_pressure']

            time_filter, pressure_filter, index = data_filtering(time,pressure,rev_time*transient_rev,rev_time * (transient_rev + n_rev))
            pressure_filter = pressure_filter - np.mean(pressure_filter)

            if case == 'Plot':

                plt.plot(time,pressure,'k')
                # plt.axvline(initialization_t, color = 'r')
                plt.grid()
                plt.xlabel('Time $[s]$',fontsize=14, style='italic')
                plt.ylabel('Pressure $[Pa]$',fontsize=14, style='italic')
                plt.xticks(fontsize=14)
                plt.yticks(fontsize=14)
                #plt.title('Time trace of ' + title)
                plt.tight_layout()
                print('Saving non-filtered figure\n')
                plt.savefig(os.path.join(image_path,'timetrace/mics/full', image_name + figure_extension), dpi=600)
                #plt.show()
                plt.close()


                plt.plot(time_filter,pressure_filter,'k')
                # plt.axvline(initialization_t, color = 'r')
                plt.grid()
                plt.xlabel('Time $[s]$',fontsize=14, style='italic')
                plt.ylabel('Pressure $[Pa]$',fontsize=14, style='italic')
                plt.xticks(fontsize=14)
                plt.yticks(fontsize=14)
                #plt.title('Time trace of ' + title)
                plt.tight_layout()
                print(f'Saving filtered figure for angle {i}-{j}\n')
                plt.savefig(os.path.join(image_path,'timetrace/mics/filtered', image_name + '-filtered' + figure_extension), dpi=600)
                #plt.show()
                plt.close()

            elif case == 'Data':
                pass
            else:
                print('Not a valid case (Data or Plot)')


    return(time, pressure)

def time_trace_psnc(rpm:int,first_it:list,second_it:list,variable:str,data_path:str,image_path:str,case:str,transient_rev=10,n_rev=2,skiprows=12):

    '''

    This function implements the ploting and analysis of the time trace of surface probes\n
    from PowerFlow simulations converted using the PowerAcoustics tool.\n
    A complete plot of the signal is done together with a close up to a user defined number\n
    of revolutions.\n
    \n
    - rpm = revolutions per minute of the rotatin machine used for the close up image\n
    - first_it = first variable to iterate over shown in the first promt of the script\n
    - second_it = second variable to iterate over shown in the first promt of the script\n
    - variable = physical variable to analyze the time trace, defaults are velocity, vorticity, pressure and density\n
    - data_path = global path where the data is located\n
    - image_path = global path where the images will be stored\n
    - case = either Plot or Data depending if the graph want to be saved or not\n
    - transient_rev = revolution where the close up will start. Prior knowledge of the signal is required\n
    - n_rev = revolutions to consider for the close up. Prior knowledge of the signal is required\n
    - skiprows = number of rows to skip in the data file. Defined due to the conversion form\n

    '''

    base_filename = 'probe_in_vanes_{-}pr_{-}z_impigning_side_export.txt' # To compare with input line 121
    print('The base file name is as:', base_filename)
    rev_time = 60/rpm

    for i in first_it:
        for j in second_it:

            file_name = f'probe_in_vanes_{i}pr_{j}z_impigning_side_export.txt' # Input to change
            image_name = file_name.split('.')[0] + '_' + variable

            data_full = pd.read_csv(os.path.join(data_path,file_name),skiprows=skiprows,delimiter='\s+',header=None)
            figure_extension = '.png'
            time = data_full[0]

            data,y_label = variable_selection(variable,data_full)

            start_time_i = next(x for x, val in enumerate(time) if val > rev_time * transient_rev)
            end_time_i = next(y for y, val in enumerate(time) if val > rev_time * (transient_rev + n_rev))
            data_filter = data[start_time_i:end_time_i]
            time_filter = time[start_time_i:end_time_i]
            data_filter = data_filter - np.mean(data_filter)

            if case == 'Plot':

                print('--- Generating image for', image_name, '---')

                plt.plot(time,data,'k')
                # plt.axvline(initialization_t, color = 'r')
                plt.grid()
                plt.xlabel('Time $[s]$',fontsize=14, style='italic')
                plt.ylabel(y_label,fontsize=14, style='italic')
                plt.xticks(fontsize=14)
                plt.yticks(fontsize=14)
                #plt.title('Time trace of ' + title)
                plt.tight_layout()
                plt.savefig(os.path.join(image_path, 'timetrace/probes/duct/full',image_name + figure_extension), dpi=600)
                #plt.show()
                plt.close()


                plt.plot(time_filter,data_filter,'k')
                # plt.axvline(initialization_t, color = 'r')
                plt.grid()
                plt.xlabel('Time $[s]$',fontsize=14, style='italic')
                plt.ylabel(y_label,fontsize=14, style='italic')
                plt.xticks(fontsize=14)
                plt.yticks(fontsize=14)
                #plt.title('Time trace of ' + title)
                plt.tight_layout()
                plt.savefig(os.path.join(image_path, 'timetrace/probes/duct/filtered',image_name + f'-{n_rev}rev-filtered' + figure_extension), dpi=600)
                #plt.show()
                plt.close()

            elif case == 'Data':
                pass
            else:
                print('Not a valid case (Data or Plot)')


    return(time, data)

def mass_flux(rpm,image_path,path,file_name,case,transient_rev=10,n_rev=1):

    image_name = file_name.split('.')[0]
    data = pd.read_csv(os.path.join(path,file_name),delimiter='\s+')
    figure_extension = '.png'
    time = data['#time']
    mass_flux = np.abs(data['mass_flux'])
    rev_time = 60/rpm

    start_time_i = next(x for x, val in enumerate(time) if val > rev_time * transient_rev)
    end_time_i = next(y for y, val in enumerate(time) if val > rev_time * (transient_rev + n_rev))
    filter_time = time[start_time_i:end_time_i]
    filter_mass = mass_flux[start_time_i:end_time_i]
    mass_flux_mean = np.mean(mass_flux[start_time_i:])

    if case == 'Plot':
        plt.plot(time, mass_flux,'k')
        plt.grid()
        plt.xlabel('Time $[s]$',fontsize=14, style='italic')
        plt.ylabel('Mass Flow $[Kg/s]$',fontsize=14, style='italic')
        #plt.title('Mass flow rate at the ' + title + ' of the system')
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        # ax=plt.gca()
        # ax.set_yticklabels([])
        # ax.yaxis.tick_right()
        plt.tight_layout()
        plt.savefig(os.path.join(image_path,'mass-flux/full',image_name+figure_extension), dpi=600)
        #plt.show()
        plt.close()

        plt.plot(filter_time,filter_mass,'k')
        plt.grid()
        plt.xlabel('Time $[s]$',fontsize=14, style='italic')
        plt.ylabel('Mass Flow $[Kg/s]$',fontsize=14, style='italic')
        #plt.title('Filtered Mass flow rate at the ' + title + ' of the system')
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        # ax=plt.gca()
        # ax.set_yticklabels([])
        # ax.yaxis.tick_right()
        plt.tight_layout()
        plt.savefig(os.path.join(image_path,'mass-flux/filtered',image_name + f'_{n_rev}rev_filtered' + figure_extension), dpi=600)
        #plt.show()
        plt.close()


    elif case == 'Data':
        pass
    else:
        print('Not a valid case (Data or Plot)')


    print('The mean mass flux after the transient data is: ', mass_flux_mean, 'kg/s')
    return(time, mass_flux)

def next_greater_power_of_2(x):
    return round(2**(x-1).bit_length())

def experimental_check(path:str,filename:str,key_name_time:str,key_name_pressure:str,image_name:str):
    
    '''
    This function could be used to quickly evaluate the signals obtained from the B&K system\n
    after they are extracted to HDF5 format.\n
    Previous knowledge of the key names is needed in order to do the  check. The first key name\n
    is related with the time data and is as follows: 'DsX1-Time', where X is the proper zero-padding\n
    to match the number of signals used. If you have less then 100 it will be Ds01-Time, otherwise\n
    Ds-001 is used. The key name pressure refers to the name B&K uses for each signal. Usually it \n
    is defined as 'DsZZ-Signal YY' where XX is related with the signal number saved for the particular\n
    case with its associated zero-padding, if a proper saving is done XX=02 if less than 100 signals are\n
    recorded. YY on the other side refers to the used port on the B&K system, i.e. if the master card is used\n
    the first signal will have YY=1. 
    
    - path = string with the full path where the data is stored. Note that the image will be stored there\n
    - filename = string with the name of the file inclduing the .h5 extention\n
    - key_name_time = string with the key name for the time data\n
    - key_name_pressure = string with the key name for the mic to analyze\n
    - image_name = string with the name for the image that wil be saved\n
     
    '''
    
    import h5py
    plt.rcParams['agg.path.chunksize'] = 10000
    
    data = h5py.File(os.path.join(path,filename),'r')
    time_data = data['Table1'][key_name_time][:]
    pressure_data = data['Table1'][key_name_pressure][:]
    
    dt=time_data[1]-time_data[0]
    fs=1.0/dt
    
    fmin = 10
    n_chunk = 20
    lensg_exp = pressure_data.size
    nperseg_exp = lensg_exp/n_chunk
    nfft_exp = next_greater_power_of_2(int(nperseg_exp))
    noverlap_exp = nperseg_exp/2

    if nperseg_exp > lensg_exp:
        raise RuntimeError('Wrong value for $f_{min}$')

    [f,Pxx]=sg.welch(pressure_data,fs=fs,window='hann',nperseg=nperseg_exp,nfft=nfft_exp,scaling='density')
    
    # Time trace plot
    plt.figure()
    plt.plot(time_data,pressure_data)
    plt.grid(True, which='both', ls='--')
    plt.xlabel('Time $[s]$',fontsize=12, style='italic')
    plt.ylabel('Pressure $[Pa]$',fontsize=12, style='italic')
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig(os.path.join(path, image_name + '-timetrace.png') , dpi=600)
    plt.close()

    # PSD plot
    plt.plot(f,10*np.log10(Pxx/4.0e-10),'k--')
    plt.grid(True, which='both', ls='--')
    plt.xlabel('Frequency $\\left[Hz\\right]$',fontsize=12, style='italic')
    plt.ylabel('PSD $\\left[\\frac{dB}{Hz}\\right]$',fontsize=12, style='italic')
    plt.xlim([100,np.max(f)])
    ax=plt.gca()
    ax.set_xscale('log')
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    print('Saving figure with levels \n')
    print('---------------------------------------------------------\n')
    plt.savefig(os.path.join(path, image_name + '-PSD.png'), dpi=600)
    plt.close()

def single_welch(n_blades:int,speed:int,angles:list,mics:list,data_path:str,image_path:str,case='Plot',filtered=0.1,mode='density'): #Change filtered to transient rev

    '''
    This function calculates the single PSD of far-field probes coming from PowerFlow\n
    simulations and previously converted using the in-house tool.\n
    
    
    '''

    print('The PSD will be calculated for a total of ', len(angles) + len(mics) ,' cases')

    BPF = n_blades * speed / 60

    for i in angles:
        for j in mics:

            print('---------------------------------------------------------\n')
            print(f'The PSD for mic {j} in the angle {i} is being calculated\n')
            print('---------------------------------------------------------\n')

            file_name = f'mic_{j}_{i}.txt'
            image_name = f'mic_{j}_{i}'

            data = pd.read_csv(os.path.join(data_path,file_name),delimiter='\s+')
            figure_extension = '.png'
            time = data['#time']
            pressure = data['static_pressure']

            # time = data['time']
            # pressure = data['pressure']

            n_chunk=4
            lensg=pressure.size
            nperseg=lensg/n_chunk
            nfft=next_greater_power_of_2(int(nperseg))

            dt=time[1]-time[0]
            fs=1.0/dt

            keep = time > filtered
            pressure = pressure[keep]

            [f,Pxx]=sg.welch(pressure,fs=fs,window='hann',nperseg=nperseg,nfft=nfft,scaling=mode)

            figure_extension = '.png'

            if case == 'Plot':

                if mode == 'density':

                    plt.plot(f/BPF,10*np.log10(Pxx/4.0e-10),'k')
                    plt.grid(True, which='both', ls='--')
                    plt.xlabel('Frequency $\\left[Hz/BPF\\right]$',fontsize=12, style='italic')
                    plt.ylabel('PSD $\\left[dB/Hz\\right]$',fontsize=12, style='italic')
                    #plt.ticklabel_format(axis='x', style='sci',scilimits=(0,0))
                    plt.xlim([100/BPF,np.max(f/BPF)])
                    #plt.title('PSD located ' + title + ' the system')
                    ax=plt.gca()
                    ax.set_yticklabels([])
                    ax.yaxis.tick_right()
                    ax.set_xscale('log')
                    ax.xaxis.set_major_formatter(FormatStrFormatter('%1.0f'))
                    plt.xticks(fontsize=12)
                    plt.yticks(fontsize=12)
                    print('Saving confidential figure\n')
                    print('---------------------------------------------------------\n')
                    plt.savefig(os.path.join(image_path, 'PSD/mics/singles', image_name + figure_extension) , dpi=600)
                    #plt.show()
                    plt.close()

                    plt.plot(f,10*np.log10(Pxx/4.0e-10),'k--')
                    plt.grid(True, which='both', ls='--')
                    plt.xlabel('Frequency $\\left[Hz\\right]$',fontsize=12, style='italic')
                    plt.ylabel('PSD $\\left[\\frac{dB}{Hz}\\right]$',fontsize=12, style='italic')
                    #plt.ticklabel_format(axis='x', style='sci',scilimits=(0,0))
                    plt.xlim([100,np.max(f)])
                    #plt.title('PSD located ' + title + ' the system')
                    ax=plt.gca()
                    #ax.set_yticklabels([])
                    #ax.yaxis.tick_right()
                    ax.set_xscale('log')
                    #ax.xaxis.set_major_formatter(FormatStrFormatter('%1.0f'))
                    plt.xticks(fontsize=12)
                    plt.yticks(fontsize=12)
                    print('Saving figure with levels \n')
                    print('---------------------------------------------------------\n')
                    plt.savefig(os.path.join(image_path, 'PSD/mics/singles', image_name + '_levels' + figure_extension), dpi=600)
                    #plt.show()
                    plt.close()

                elif mode == 'spectrum':

                    plt.plot(f,10*np.log10(Pxx/4.0e-10),'k--')
                    plt.grid(True, which='both', ls='--')
                    plt.xlabel('Frequency $\\left[Hz\\right]$',fontsize=12, style='italic')
                    plt.ylabel('SPL $\\left[dB\\right]$',fontsize=12, style='italic')
                    #plt.ticklabel_format(axis='x', style='sci',scilimits=(0,0))
                    plt.xlim([100,np.max(f)])
                    #plt.title('SPL located ' + title + ' the system')
                    ax=plt.gca()
                    ax.set_xscale('log')
                    #ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
                    plt.xticks(fontsize=12)
                    plt.yticks(fontsize=12)
                    print('Saving figure')
                    print('---------------------------------------------------------\n')
                    plt.savefig(image_path + image_name + figure_extension , dpi=1200)
                    #plt.show()
                    plt.close()

                else:
                    print('Not a valid mode (density or spectrum)')

            elif case == 'Data':
                pass

            else:
                print('Not a valid case (Data or Plot)')

            print('Finished process\n')

    return(f,Pxx)

def welch_all(n_blades,speed,first_it,second_it,pk,exp_mics,exp_data_path,sim_data_path,filtered,mode,*geometries):

    os.makedirs(os.path.join(sim_data_path, 'images/psd'), exist_ok=True)

    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']

    BPF = n_blades * speed / 60

    sim_check = int(input('If you want to plot the direct acoustic data enter 1, otherwise enter 0 \n'))
    exp_check = int(input('If you want to plot the experimental acoustic data enter 1, otherwise enter 0 \n'))
    pa_check = int(input('If you want to plot the hybrid acoustic data obtained from PowerAcoustics enter 1, otherwise enter 0 \n'))

    for i in first_it:
        for j,z in zip(second_it,exp_mics):

            #fig, ax1 = plt.subplots(figsize=(15, 4))
            image_name = f'frec_ffn_mic_{i}_{j}_2'

            print('---------------------------------------------------------\n')
            print(f'The PSD for mic {j} in the angle {i} is being calculated\n')
            print('---------------------------------------------------------\n')

            

            file_name = f'mic_{i}_{j}.txt'

            data = pd.read_csv(os.path.join(sim_data_path, file_name),delimiter='\s+')
            figure_extension = '.png'
            time = data['#time']
            pressure = data['static_pressure']

            # time = data['time']
            # pressure = data['pressure']

            n_chunk=4
            lensg=pressure.size
            nperseg=int(lensg/n_chunk)
            nfft=next_greater_power_of_2(nperseg)

            dt=time[1]-time[0]
            fs=1.0/dt

            keep = time > filtered
            pressure = pressure[keep]
            # pdb.set_trace()
            [f_sim,Pxx_sim]=sg.welch(pressure.values,fs=fs,window='hann',nperseg=nperseg,nfft=nfft,scaling=mode)

            f = f_sim
            
            if sim_check == 1:
                image_name = image_name + '_direct'
                plt.plot(f_sim/BPF,10*np.log10(Pxx_sim/4.0e-10),'k',label='LBM Direct Acoustics')

            else:
                pass

            if exp_check == 1:

                exp_file_name = f'{speed}rpm_{i}deg_pk{pk}.h5'
                image_name = image_name + '_exp'
                pressure_exp, dt_exp, f_exp = read_mic_file(os.path.join(exp_data_path,exp_file_name),z)
                time_exp = np.arange(pressure_exp.shape[-1]) * dt_exp


                idt_final = np.where(time_exp >= time.iloc[-1])[0][0]
                idt_start = np.where(time_exp >= time.iloc[0])[0][0]

                pressure_exp = pressure_exp[:,idt_start:idt_final]
                f_exp = f_exp[idt_start:idt_final]

                fmin = 10
                lensg_exp = pressure_exp.size
                #nperseg_exp = int(1/(fmin/5*dt_exp))
                nperseg_exp = lensg_exp/n_chunk
                nfft_exp = next_greater_power_of_2(int(nperseg_exp))
                noverlap_exp = nperseg_exp/2

                if nperseg_exp > lensg_exp:
                    raise RuntimeError('Wrong value for $f_{min}$')

                fs_exp=1.0/dt_exp

                [f_exp,Pxx_exp]=sg.welch(pressure_exp,fs=fs_exp,nperseg=nperseg_exp,nfft=nfft_exp,noverlap=noverlap_exp,scaling=mode)

                plt.plot(f_exp/BPF,10*np.log10(Pxx_exp[0]/4.0e-10),'b',label='Experiment UdeS')

            else:
                pass

            if pa_check == 1:
                labels = ['Rotor','Stator']
                image_name = image_name + '_analogy'
                c = 0
                max_case = 0

                for g,label in zip(geometries,labels):

                    c += 1

                    file_name = f'fwh-nosimp-{g}.mic_{i}_{j}.txt'
                    data_pa = pd.read_csv(os.path.join(sim_data_path, file_name),skiprows=4,delimiter='\s+',header=None)

                    globals()[f'time_pa_{g}'] = data_pa[0]
                    globals()[f'pressure_pa_{g}'] = data_pa[1].values
                    globals()[f'{g}_size'] = globals()[f'pressure_pa_{g}'].size

                    if globals()[f'{g}_size'] >= max_case:

                        max_case = globals()[f'{g}_size']
                        name_max_case = g

                    else:
                        pass

                    n_chunk = 4
                    lensg=globals()[f'pressure_pa_{g}'].size
                    nperseg=lensg/n_chunk
                    nfft=next_greater_power_of_2(int(nperseg))

                    dt=globals()[f'time_pa_{g}'][1]-globals()[f'time_pa_{g}'][0]
                    fs=1.0/dt

                    [globals()[f'f_{g}'],globals()[f'Pxx_{g}']] = sg.welch(globals()[f'pressure_pa_{g}'],fs=fs,window='hann',nperseg=nperseg,nfft=nfft,scaling=mode)

                    plt.plot(globals()[f'f_{g}']/BPF,10*np.log10(globals()[f'Pxx_{g}']/4.0e-10),color = colors[c],label=f'LBM FW-H Acoustics {label}')

                pressure_total = np.zeros(max_case)

                for g in geometries:

                    if g == name_max_case:
                        pressure_total = pressure_total + globals()[f'pressure_pa_{g}']
                    else:
                        globals()[f'pressure_inter_{g}'] = np.interp(globals()[f'time_pa_{name_max_case}'],globals()[f'time_pa_{g}'],globals()[f'pressure_pa_{g}'])
                        pressure_total = pressure_total + globals()[f'pressure_inter_{g}']

                n_chunk = 4
                lensg = pressure_total.size
                nperseg = lensg/n_chunk
                nfft = next_greater_power_of_2(int(nperseg))
                dt = globals()[f'time_pa_{name_max_case}'][1] - globals()[f'time_pa_{name_max_case}'][0]
                fs = 1.0/dt

                [f_total,Pxx_total] = sg.welch(pressure_total,fs=fs,window='hann',nperseg=nperseg,nfft=nfft,scaling=mode)
                f = f_total

                plt.plot(f_total/BPF,10*np.log10(Pxx_total/4.0e-10),color = colors[c+1],label=f'LBM FW-H Acoustics Full Geometry')

            else:
                pass


            if mode == 'density':

                plt.grid(True, which='both', ls='--')
                # plt.legend(prop={'size': 10})
                plt.xlabel('Frequency $\\left[Hz/BPF\\right]$',fontsize=14, style='italic')
                plt.ylabel('PSD $\\left[dB/Hz\\right]$',fontsize=14, style='italic')
                #plt.ticklabel_format(axis='x', style='sci',scilimits=(0,0))
                plt.ylim([10,85])
                plt.xlim([100/BPF,np.max(f_sim)/BPF])
                #plt.title('PSD located ' + title + ' the system')
                ax=plt.gca()
                ax.set_yticks([20,30,40,50,60,70,80])
                ax.set_yticklabels(['','','','','$r_f-10$','$r_f$','$r_f+10$'])
                ax.yaxis.tick_right()
                ax.set_xscale('log')
                ax.xaxis.set_major_formatter(FormatStrFormatter('%1.0f'))
                #ax.yaxis.tick_right()
                plt.xticks(fontsize=12)
                plt.yticks(fontsize=12)
                plt.tight_layout()
                print('Saving confidential figure\n')
                print('---------------------------------------------------------\n')
                plt.savefig(os.path.join(sim_data_path, 'images/psd', image_name + '_confidential_density.png'), dpi=600)
                # plt.show()
                plt.close()
            #     plt.plot(f,10*np.log10(Pxx/4.0e-10),'k',label='LBM simulation')
            #     plt.plot(frequency_exp,10*np.log10(data_exp[z])+cor[z],'b',label='Experiment ECL')
            #     plt.grid(True, which='both', ls='--')
            #     plt.legend()
            #     plt.xlabel('Frequency $\\left[Hz\\right]$',fontsize=12, style='italic')
            #     plt.ylabel('PSD $\\left[\\frac{dB}{Hz}\\right]$',fontsize=12, style='italic')
            #     #plt.ticklabel_format(axis='x', style='sci',scilimits=(0,0))
            #     plt.xlim([100,np.max(f)])
            #     #plt.title('PSD located ' + title + ' the system')
            #     ax=plt.gca()
            #     #ax.set_yticklabels([])
            #     #ax.yaxis.tick_right()
            #     ax.set_xscale('log')
            #     #ax.xaxis.set_major_formatter(FormatStrFormatter('%1.0f'))
            #     plt.xticks(fontsize=12)
            #     plt.yticks(fontsize=12)
            #     plt.tight_layout()
            #     print('Saving figure with levels \n')
            #     print('---------------------------------------------------------\n')
            #     plt.savefig(image_path + image_name + '_density-levels' + figure_extension , dpi=600)
            #     #plt.show()
            #     plt.close()

            # elif mode == 'spectrum':

            #     plt.plot(f,10*np.log10(Pxx/4.0e-10),'k',label='LBM simulation')
            #     plt.plot(frequency_exp,10*np.log10(data_exp[z])+cor[z],'b',label='Experiment ECL')
            #     plt.grid(True, which='both', ls='--')
            #     plt.legend()
            #     plt.xlabel('Frequency $\\left[Hz\\right]$',fontsize=12, style='italic')
            #     plt.ylabel('SPL $\\left[dB\\right]$',fontsize=12, style='italic')
            #     #plt.ticklabel_format(axis='x', style='sci',scilimits=(0,0))
            #     plt.xlim([100,np.max(f)])
            #     #plt.title('SPL located ' + title + ' the system')
            #     ax=plt.gca()
            #     ax.set_xscale('log')
            #     #ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
            #     plt.xticks(fontsize=12)
            #     plt.yticks(fontsize=12)
            #     plt.tight_layout()
            #     print('Saving figure')
            #     print('---------------------------------------------------------\n')
            #     plt.savefig(image_path + image_name + 'spectrum' + figure_extension , dpi=600)
            #     #plt.show()
            #     plt.close()

            #     plt.plot(f/BPF,10*np.log10(Pxx/4.0e-10),'k',label='LBM simulation')
            #     plt.plot(frequency_exp/BPF_exp,10*np.log10(data_exp[z])+cor[z],'b',label='Experiment ECL')
            #     plt.grid(True, which='both', ls='--')
            #     plt.legend()
            #     plt.xlabel('Frequency $\\left[Hz/BPF\\right]$',fontsize=12, style='italic')
            #     plt.ylabel('PSD $\\left[dB\\right]$',fontsize=12, style='italic')
            #     #plt.ticklabel_format(axis='x', style='sci',scilimits=(0,0))
            #     plt.xlim([100/BPF,np.max(f/BPF)])
            #     #plt.title('PSD located ' + title + ' the system')
            #     ax=plt.gca()
            #     #ax.set_yticklabels([])
            #     #ax.yaxis.tick_right()
            #     ax.set_xscale('log')
            #     ax.xaxis.set_major_formatter(FormatStrFormatter('%1.0f'))
            #     plt.xticks(fontsize=12)
            #     plt.yticks(fontsize=12)
            #     plt.tight_layout()
            #     print('Saving confidential figure\n')
            #     print('---------------------------------------------------------\n')
            #     plt.savefig(image_path + image_name + 'spectrum-BPF' + figure_extension , dpi=600)
            #     #plt.show()
            #     plt.close()

            else:
                print('Not a valid mode (density or spectrum)')


def sherfwhwelch(n_blades,speed,mics,data_path,image_path,geometry,case,mode):

    print('The PSD will be calculated for a total of ', len(mics) ,' cases')

    BPF = n_blades * speed / 60

    for i in mics:

        print('---------------------------------------------------------\n')
        print(f'The PSD for mic {i} is being calculated\n')
        print('---------------------------------------------------------\n')

        file_name = f'Mic_{i}'
        image_name = f'fwh-mic_{i}_'

        data = pd.read_csv(data_path + file_name,delimiter='\s+',skiprows=1,names=['iter','time','pressure','nbcontrib'])
        figure_extension = '.png'
        time = data['time']
        pressure = data['pressure']
        contribution = data['nbcontrib']
        filtered = np.where(contribution==contribution[contribution.size-1])[0][0]

        # time = data['time']
        # pressure = data['pressure']

        n_chunk=4
        lensg=pressure.size
        nperseg=lensg/n_chunk
        nfft=next_greater_power_of_2(int(nperseg))

        dt=time[1]-time[0]
        fs=1.0/dt

        time = time[filtered:]
        pressure = pressure[filtered:]
        pressure = pressure - np.mean(pressure)

        [f,Pxx]=sg.welch(pressure,fs=fs,window='hann',nperseg=nperseg,nfft=nfft,scaling=mode)

        figure_extension = '.png'

        if case == 'Plot':

            if mode == 'density':

                plt.plot(f/BPF,10*np.log10(Pxx/4.0e-10),'k--')
                plt.grid(True, which='both', ls='--')
                plt.xlabel('Frequency $\\left[Hz/BPF\\right]$',fontsize=12, style='italic')
                plt.ylabel('PSD $\\left[dB/Hz\\right]$',fontsize=12, style='italic')
                #plt.ticklabel_format(axis='x', style='sci',scilimits=(0,0))
                plt.xlim([100/BPF,np.max(f/BPF)])
                #plt.title('PSD located ' + title + ' the system')
                ax=plt.gca()
                ax.set_yticklabels([])
                ax.yaxis.tick_right()
                ax.set_xscale('log')
                ax.xaxis.set_major_formatter(FormatStrFormatter('%1.0f'))
                plt.xticks(fontsize=12)
                plt.yticks(fontsize=12)
                print('Saving confidential figure\n')
                print('---------------------------------------------------------\n')
                plt.savefig(image_path + image_name + geometry + figure_extension , dpi=600)
                #plt.show()
                plt.close()

                plt.plot(f,10*np.log10(Pxx/4.0e-10),'k--')
                plt.grid(True, which='both', ls='--')
                plt.xlabel('Frequency $\\left[Hz\\right]$',fontsize=12, style='italic')
                plt.ylabel('PSD $\\left[\\frac{dB}{Hz}\\right]$',fontsize=12, style='italic')
                #plt.ticklabel_format(axis='x', style='sci',scilimits=(0,0))
                plt.xlim([100,np.max(f)])
                #plt.title('PSD located ' + title + ' the system')
                ax=plt.gca()
                #ax.set_yticklabels([])
                #ax.yaxis.tick_right()
                ax.set_xscale('log')
                #ax.xaxis.set_major_formatter(FormatStrFormatter('%1.0f'))
                plt.xticks(fontsize=12)
                plt.yticks(fontsize=12)
                print('Saving figure with levels \n')
                print('---------------------------------------------------------\n')
                plt.savefig(image_path + image_name + geometry + '_levels' + figure_extension , dpi=600)
                #plt.show()
                plt.close()

            elif mode == 'spectrum':

                plt.plot(f,10*np.log10(Pxx/4.0e-10),'k--')
                plt.grid(True, which='both', ls='--')
                plt.xlabel('Frequency $\\left[Hz\\right]$',fontsize=12, style='italic')
                plt.ylabel('SPL $\\left[dB\\right]$',fontsize=12, style='italic')
                #plt.ticklabel_format(axis='x', style='sci',scilimits=(0,0))
                plt.xlim([100,np.max(f)])
                #plt.title('SPL located ' + title + ' the system')
                ax=plt.gca()
                ax.set_xscale('log')
                #ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
                plt.xticks(fontsize=12)
                plt.yticks(fontsize=12)
                print('Saving figure')
                print('---------------------------------------------------------\n')
                plt.savefig(image_path + image_name + figure_extension , dpi=1200)
                #plt.show()
                plt.close()

            else:
                print('Not a valid mode (density or spectrum)')

        elif case == 'Data':
            pass

        else:
            print('Not a valid case (Data or Plot)')

        print('Finished process\n')

    return(f,Pxx)

def psd_psnc(n_blades,rpm,variable,first_it,second_it,data_path,image_path,case,transient_rev=5,skiprows=12):

    #probe_type = str(input('Enter if you want to perform the analysis of a surface or experimental probe\n'))
    probe_type = 'surface'

    base_filename = 'probe_in_vanes_{-}pr_{-}z_impigning_side_export.txt' # To compare with input line 121
    print('The base file name is as:', base_filename)
    # dark_colors = ['#000000', '#525252', '#7F7F7F', '#AEABAB', '#D9D9D9']
    # k = 0

    BPF = n_blades * rpm / 60
    rev_time = 60/rpm

    for i in first_it:
        for j in second_it:

            if probe_type == 'surface':

                file_name = f'probe_in_vanes_{i}pr_{j}z_no_impigning_side_export.txt'
                image_name = file_name.split('.')[0]

                data_full = pd.read_csv(os.path.join(data_path,file_name),skiprows=skiprows,delimiter='\s+',header=None)
                time = data_full[0]

                data,y_label = variable_selection(variable,data_full)

            elif probe_type == 'experimental':

                warnings.warn('The value entered has not yet been implemented in this code version')
                sys.exit('Program terminated due to invalid warning')

            else:

                raise ValueError('Not a valid type of probe. Code terminated')

            start_time_i = next(x for x, val in enumerate(time) if val > rev_time * transient_rev)
            time = time[start_time_i:]
            data = data[start_time_i:]

            n_chunk=4
            lensg=data.size
            nperseg=lensg/n_chunk
            nfft=next_greater_power_of_2(int(nperseg))

            dt=time[start_time_i+1]-time[start_time_i]
            fs=1.0/dt

            [f,Pxx]=sg.welch(data,fs=fs,window='hann',nperseg=nperseg,nfft=nfft,scaling='density')

            if case == 'Plot':
                # pdb.set_trace()
                print('--- Generating image for', image_name, '---')

                # plt.plot(f/BPF,10*np.log10(Pxx/4.0e-10),label=f'$r={i}\\%$ $z={j}\\%$')#,color=dark_colors[k])
                # plt.grid(True, which='both', ls='--')
                # plt.xlabel('BPF $\\left[-\\right]$',fontsize=14, style='italic')
                # plt.ylabel('$\\Phi_{pp} \\left[dB/Hz\\right]$',fontsize=14, style='italic')
                # #plt.ticklabel_format(axis='x', style='sci',scilimits=(0,0))
                # plt.xlim([100/BPF,np.max(f/BPF)])
                # plt.legend(loc='upper right', ncol=1, fontsize='small')
                # #plt.title('PSD located ' + title + ' the system')
                # ax=plt.gca()
                # ax.set_yticklabels([])
                # ax.yaxis.tick_left()
                # ax.set_xscale('log')
                # ax.xaxis.set_major_formatter(FormatStrFormatter('%1.0f'))
                # plt.xticks(fontsize=12)
                # plt.yticks(fontsize=12)


                plt.plot(f/BPF,10*np.log10(Pxx/4.0e-10),'k')
                plt.grid(True, which='both', ls='--')
                plt.xlabel('BPF $\\left[-\\right]$',fontsize=14, style='italic')
                plt.ylabel('$\\Phi_{pp} \\left[\\frac{dB}{Hz}\\right]$',fontsize=14, style='italic')
                #plt.ticklabel_format(axis='x', style='sci',scilimits=(0,0))
                plt.xlim([100/BPF,np.max(f/BPF)])
                #plt.title('PSD located ' + title + ' the system')
                ax=plt.gca()
                #ax.set_yticklabels([])
                #ax.yaxis.tick_right()
                ax.set_xscale('log')
                ax.xaxis.set_major_formatter(FormatStrFormatter('%1.0f'))
                plt.xticks(fontsize=12)
                plt.yticks(fontsize=12)
                print('Saving figure with levels \n')
                print('---------------------------------------------------------\n')
                plt.tight_layout()
                plt.savefig(os.path.join(image_path,'PSD/probes/duct/singles',image_name + f'{variable}_PSD.png') , dpi=600)
                #plt.show()
                plt.close()

            elif case == 'Data':
                pass

            else:
                print('Not a valid case (Data or Plot)')

            # k += 1

        print('Saving confidential figure\n')
        print('---------------------------------------------------------\n')

    print('Finished process\n')
    # plt.tight_layout()
    # plt.savefig(os.path.join(image_path,'PSD/probes/duct/comparisons',f'probe_in_duct_{j}deg-comparison-rotor-plane-PSD.png') , dpi=600)
    # #plt.show()
    # plt.close()

    return(f,Pxx)

def psd_comparison(n_blades,speed,pk,angles,mics_sim,mics_exp,sim_data_path,exp_data_path,image_path,case,filtered,mode):

    """

    For .hdf5 exported directly from B&K files and .txt files converted from .pfnc\n
    of PowerFlow

    """

    print('The PSD will be calculated for a total of ', len(angles) + len(mics_exp) ,' comparison cases')

    for i in angles:
        for j,k in zip(mics_sim,mics_exp):

            print('---------------------------------------------------------\n')
            print(f'The PSD for mic {j} in the angle {i} is being calculated\n')
            print('---------------------------------------------------------\n')

            sim_file_name = f'mic_{i}_{j}.txt'
            exp_file_name = f'{speed}rpm_{i}deg_pk{pk}.h5'
            image_name = f'comparison_mic_{j}_{i}'

            data_sim = pd.read_csv(os.path.join(sim_data_path, sim_file_name),delimiter='\s+')
            figure_extension = '.png'
            time_sim = data_sim['#time']
            pressure_sim = data_sim['static_pressure']
            # time = data['time']
            # pressure = data['pressure']

            BPF = n_blades * speed / 60

            n_chunk=4
            lensg_sim=pressure_sim.size
            nperseg_sim=lensg_sim/n_chunk
            nfft_sim=next_greater_power_of_2(int(nperseg_sim))

            dt_sim=time_sim[1]-time_sim[0]
            fs_sim=1.0/dt_sim

            keep = time_sim > filtered
            pressure_sim = pressure_sim[keep]
            sim_len = len(pressure_sim)

            [f_sim,Pxx_sim]=sg.welch(pressure_sim,fs=fs_sim,window='hann',nperseg=nperseg_sim,nfft=nfft_sim,scaling=mode)

            pressure_exp, dt_exp, f_exp = read_mic_file(exp_data_path + exp_file_name,k)
            time_exp = np.arange(pressure_exp.shape[-1]) * dt_exp


            idt_final = np.where(time_exp >= time_sim.iloc[-1])[0][0]
            idt_start = np.where(time_exp >= time_sim.iloc[0])[0][0]

            pressure_exp = pressure_exp[:,idt_start:idt_final]
            f_exp = f_exp[idt_start:idt_final]

            fmin = 10
            lensg_exp = pressure_exp.size
            #nperseg_exp = int(1/(fmin/5*dt_exp))
            nperseg_exp = lensg_exp/n_chunk
            nfft_exp = next_greater_power_of_2(int(nperseg_exp))
            noverlap_exp = nperseg_exp/2

            if nperseg_exp > lensg_exp:
                raise RuntimeError('Wrong value for $f_{min}$')

            fs_exp=1.0/dt_exp
            
            [f_exp,Pxx_exp]=sg.welch(pressure_exp,fs=fs_exp,nperseg=nperseg_exp,nfft=nfft_exp,noverlap=noverlap_exp,scaling=mode)

            figure_extension = '.png'

            if case == 'Plot':

                if mode == 'density':

                    plt.plot(f_sim/BPF,10*np.log10(Pxx_sim/4.0e-10),'k--',label='Simulations')
                    plt.plot(f_exp/BPF,10*np.log10(Pxx_exp[0]/4.0e-10),'b--',label='Experiments')
                    plt.grid(True, which='both', ls='--')
                    plt.legend()
                    plt.xlabel('BPF $\\left[-\\right]$',fontsize=12, style='italic')
                    plt.ylabel('PSD $\\left[dB/Hz\\right]$',fontsize=12, style='italic')
                    #plt.ticklabel_format(axis='x', style='sci',scilimits=(0,0))
                    plt.xlim([100/BPF,np.max(f_sim)/BPF])
                    #plt.ylim([15,1.3*np.max(10*np.log10(Pxx_sim/4.0e-10))])
                    plt.ylim([0,90])
                    #plt.title('PSD located ' + title + ' the system')
                    ax=plt.gca()
                    # ax.set_yticklabels([])
                    # ax.yaxis.tick_right()
                    ax.set_xscale('log')
                    ax.xaxis.set_major_formatter(FormatStrFormatter('%1.0f'))
                    plt.xticks(fontsize=12)
                    plt.yticks([10,25,40,55,70,85],['$r_f-60$','$r_f-45$','$r_f-30$','$r_f-15$','$r_f$','$r_f+15$'],fontsize=12)
                    plt.tight_layout()
                    print('Saving confidential figure\n')
                    print('---------------------------------------------------------\n')
                    plt.savefig(os.path.join(image_path,image_name + figure_extension) , dpi=600)
                    #plt.show()
                    plt.close()

                    plt.plot(f_sim,10*np.log10(Pxx_sim/4.0e-10),'k--',label='Simulations')
                    plt.plot(f_exp,10*np.log10(Pxx_exp[0]/4.0e-10),'b--',label='Experiments')
                    plt.grid(True, which='both', ls='--')
                    plt.legend()
                    plt.xlabel('Frequency $\\left[Hz\\right]$',fontsize=12, style='italic')
                    plt.ylabel('PSD $\\left[\\frac{dB}{Hz}\\right]$',fontsize=12, style='italic')
                    #plt.ticklabel_format(axis='x', style='sci',scilimits=(0,0))
                    plt.xlim([100,np.max(f_sim)])
                    #plt.ylim([15,1.3*np.max(10*np.log10(Pxx_sim/4.0e-10))])
                    plt.ylim([0,90])
                    #plt.title('PSD located ' + title + ' the system')
                    ax=plt.gca()
                    ax.set_xscale('log')
                    #ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
                    plt.xticks(fontsize=12)
                    plt.yticks(fontsize=12)
                    print('Saving figure with levels \n')
                    print('---------------------------------------------------------\n')
                    plt.savefig(os.path.join(image_path, image_name + '_levels' + figure_extension) , dpi=600)
                    #plt.show()
                    plt.close()

                elif mode == 'spectrum':

                    plt.plot(f_sim,10*np.log10(Pxx_sim/4.0e-10),'k--',label='Simulations')
                    plt.plot(f_exp,10*np.log10(Pxx_exp[0]/4.0e-10),'b--',label='Experiments')
                    plt.grid(True, which='both', ls='--')
                    plt.legend()
                    plt.xlabel('Frequency $\\left[Hz\\right]$',fontsize=12, style='italic')
                    plt.ylabel('SPL $\\left[dB\\right]$',fontsize=12, style='italic')
                    #plt.ticklabel_format(axis='x', style='sci',scilimits=(0,0))
                    plt.xlim([100,np.max(f_sim)])
                    plt.ylim([15,1.3*np.max(10*np.log10(Pxx_sim/4.0e-10))])
                    #plt.title('SPL located ' + title + ' the system')
                    ax=plt.gca()
                    ax.set_xscale('log')
                    #ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
                    plt.xticks(fontsize=12)
                    plt.yticks(fontsize=12)
                    print('Saving figure')
                    print('---------------------------------------------------------\n')
                    plt.savefig(image_path + image_name + figure_extension , dpi=600)
                    #plt.show()
                    plt.close()

                else:
                    print('Not a valid mode (density or spectrum)')

            elif case == 'Data':
                pass

            else:
                print('Not a valid case (Data or Plot)')

            print('Finished process\n')

    return(f_sim,Pxx_sim)

def simcomparisonpsd(n_blades,speed,angles,mics_sim1,mics_sim2,sim1_data_path,sim2_data_path,image_path,case,filtered,mode):

    """

    For two .txt files converted from .pfnc format of PowerFlow\n

    """

    print('The PSD will be calculated for a total of ', len(angles) + len(mics_sim1) ,' comparison cases')

    for i in angles:
        for j,k in zip(mics_sim1,mics_sim2):

            print('---------------------------------------------------------\n')
            print(f'The PSD for mic {j} in the angle {i} is being calculated\n')
            print('---------------------------------------------------------\n')

            sim1_file_name = f'mic{j}_{i}.txt'
            sim2_file_name = f'mic{j}_{i}.txt'
            image_name = f'2023-comparison_mic{j}_{i}'



            data_sim1 = pd.read_csv(sim1_data_path + sim1_file_name,delimiter='\s+')
            figure_extension = '.png'
            time_sim1 = data_sim1['#time']
            pressure_sim1 = data_sim1['static_pressure']
            # time = data['time']
            # pressure = data['pressure']

            BPF = n_blades * speed / 60

            n_chunk1=4
            lensg_sim1=pressure_sim1.size
            nperseg_sim1=lensg_sim1/n_chunk1
            nfft_sim1=next_greater_power_of_2(int(nperseg_sim1))

            dt_sim1=time_sim1[1]-time_sim1[0]
            fs_sim1=1.0/dt_sim1

            keep1 = time_sim1 > filtered
            pressure_sim1 = pressure_sim1[keep1]
            sim_len1 = len(pressure_sim1)

            [f_sim1,Pxx_sim1]=sg.welch(pressure_sim1,fs=fs_sim1,window='hann',nperseg=nperseg_sim1,nfft=nfft_sim1,scaling=mode)

            data_sim2 = pd.read_csv(sim2_data_path + sim2_file_name,delimiter='\s+')
            time_sim2 = data_sim2['#time']
            pressure_sim2 = data_sim2['static_pressure']
            # time = data['time']
            # pressure = data['pressure']

            n_chunk2=4
            lensg_sim2=pressure_sim2.size
            nperseg_sim2=lensg_sim2/n_chunk2
            nfft_sim2=next_greater_power_of_2(int(nperseg_sim2))

            dt_sim2=time_sim2[1]-time_sim2[0]
            fs_sim2=1.0/dt_sim2

            keep2 = time_sim2 > filtered
            pressure_sim2 = pressure_sim2[keep2]
            sim_len2 = len(pressure_sim2)

            [f_sim2,Pxx_sim2]=sg.welch(pressure_sim2,fs=fs_sim2,window='hann',nperseg=nperseg_sim2,nfft=nfft_sim2,scaling=mode)


            if case == 'Plot':

                if mode == 'density':

                    plt.plot(f_sim1/BPF,10*np.log10(Pxx_sim1/4.0e-10),'k--',label='2024')
                    plt.plot(f_sim2/BPF,10*np.log10(Pxx_sim2/4.0e-10),'b--',label='2023')
                    plt.grid(True, which='both', ls='--')
                    plt.legend()
                    plt.xlabel('Frequency $\\left[Hz/BPF\\right]$',fontsize=12, style='italic')
                    plt.ylabel('PSD $\\left[dB/Hz\\right]$',fontsize=12, style='italic')
                    #plt.ticklabel_format(axis='x', style='sci',scilimits=(0,0))
                    plt.xlim([100/BPF,np.max(f_sim1)/BPF])
                    plt.ylim([15,1.3*np.max(10*np.log10(Pxx_sim1/4.0e-10))])
                    #plt.title('PSD located ' + title + ' the system')
                    ax=plt.gca()
                    ax.set_yticklabels([])
                    ax.yaxis.tick_right()
                    ax.set_xscale('log')
                    ax.xaxis.set_major_formatter(FormatStrFormatter('%1.0f'))
                    plt.xticks(fontsize=12)
                    plt.yticks(fontsize=12)
                    print('Saving confidential figure\n')
                    print('---------------------------------------------------------\n')
                    plt.savefig(image_path + image_name + figure_extension , dpi=600)
                    #plt.show()
                    plt.close()

                    plt.plot(f_sim1,10*np.log10(Pxx_sim1/4.0e-10),'k--',label='2024')
                    plt.plot(f_sim2,10*np.log10(Pxx_sim2/4.0e-10),'b--',label='2023')
                    plt.grid(True, which='both', ls='--')
                    plt.legend()
                    plt.xlabel('Frequency $\\left[Hz\\right]$',fontsize=12, style='italic')
                    plt.ylabel('PSD $\\left[\\frac{dB}{Hz}\\right]$',fontsize=12, style='italic')
                    #plt.ticklabel_format(axis='x', style='sci',scilimits=(0,0))
                    plt.xlim([100,np.max(f_sim1)])
                    plt.ylim([15,1.3*np.max(10*np.log10(Pxx_sim1/4.0e-10))])
                    #plt.title('PSD located ' + title + ' the system')
                    ax=plt.gca()
                    #ax.set_yticklabels([])
                    #ax.yaxis.tick_right()
                    ax.set_xscale('log')
                    #ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
                    plt.xticks(fontsize=12)
                    plt.yticks(fontsize=12)
                    print('Saving figure with levels \n')
                    print('---------------------------------------------------------\n')
                    plt.savefig(image_path + image_name + '_levels' + figure_extension , dpi=600)
                    #plt.show()
                    plt.close()

                elif mode == 'spectrum':

                    plt.plot(f_sim1,10*np.log10(Pxx_sim1/4.0e-10),'k--',label='Pk 103')
                    plt.plot(f_sim2,10*np.log10(Pxx_sim2/4.0e-10),'b--',label='Pk 109')
                    plt.grid(True, which='both', ls='--')
                    plt.legend()
                    plt.xlabel('Frequency $\\left[Hz\\right]$',fontsize=12, style='italic')
                    plt.ylabel('SPL $\\left[dB\\right]$',fontsize=12, style='italic')
                    #plt.ticklabel_format(axis='x', style='sci',scilimits=(0,0))
                    plt.xlim([100,np.max(f_sim1)])
                    plt.ylim([15,1.3*np.max(10*np.log10(Pxx_sim1/4.0e-10))])
                    #plt.title('SPL located ' + title + ' the system')
                    ax=plt.gca()
                    ax.set_xscale('log')
                    #ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
                    plt.xticks(fontsize=12)
                    plt.yticks(fontsize=12)
                    print('Saving figure')
                    print('---------------------------------------------------------\n')
                    plt.savefig(image_path + image_name + figure_extension , dpi=600)
                    #plt.show()
                    plt.close()

                else:
                    print('Not a valid mode (density or spectrum)')

            elif case == 'Data':
                pass

            else:
                print('Not a valid case (Data or Plot)')

            print('Finished process\n')

    return(f_sim1,Pxx_sim1)

def welch_exp(image_path,path,variable,case,mode):

    data = pd.read_csv(path,delimiter='\s+')
    figure_extension = '.png'
    f = data['#time']
    Pxx = data['static_pressure']

    figure_extension = '.png'

    if case == 'Plot':

        if mode == 'density':

            plt.plot(f,10*np.log10(Pxx/4.0e-10),'k--')
            plt.grid(True, which='both', ls='--')
            plt.xlabel('Frequency $[Hz]$')
            plt.ylabel('PSD $[dB/Hz]$')
            #plt.ticklabel_format(axis='x', style='sci',scilimits=(0,0))
            plt.xlim([100,10000])
            #plt.title('PSD located ' + title + ' the system')
            ax=plt.gca()
            ax.set_yticklabels([])
            ax.yaxis.tick_right()
            ax.set_xscale('log')
            ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
            plt.savefig(image_path + variable + 'max=' + np.str(np.round(np.max(10*np.log10(Pxx/4.0e-10)),decimals=2)) + figure_extension , dpi=1200)
            #plt.show()
            plt.close()

        elif mode == 'spectrum':

            plt.plot(f,10*np.log10(Pxx/4.0e-10),'k--')
            plt.grid(True, which='both', ls='--')
            plt.xlabel('Frequency $[Hz]$')
            plt.ylabel('SPL $[dB]$')
            #plt.ticklabel_format(axis='x', style='sci',scilimits=(0,0))
            plt.xlim([100,10000])
            #plt.title('SPL located ' + title + ' the system')
            ax=plt.gca()
            ax.set_xscale('log')
            ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
            plt.savefig(image_path + variable + figure_extension , dpi=1200)
            #plt.show()
            plt.close()

        else:
            print('Not a valid mode (density or spectrum)')

    elif case == 'Data':
        pass

    else:
        print('Not a valid case (Data or Plot)')

    return(f,Pxx)

def forces_calculation(rpm:int,component:str,axis:str,format:str,data_path:str,file_name:str,case:str,transient_rev=10,n_rev=2,skip=18):

    '''
    This function is used to calculate both torque and thrust from PowerFlow simulations.\n
    The input parameters need to be a .txt file output from the exaritool method.\n
    \n
    - component = upper case either F if the analysis is for thrust or M if it is for torque\n
    - axis = (x,y,z) lower case with the axis aligend with the interest force\n
    - image_folder = with the path where the images will be stored\n
    - data_path = with the path where the data is stored\n
    - file_name = with the name of the file that wants to be analyzed\n
    - case = either Plot or Data depending if the data wants to be ploted or not\n
    - trigger = a number where the filterring wants to be done, if not specified 0.1s will be the lower limit\n
    - skip = the number of lines to skip in the .txt file. 18 is used as default\n

    '''
    
    os.makedirs(os.path.join(os.path.dirname(data_path), f'images/forces/filtered/{component}'), exist_ok=True)
    os.makedirs(os.path.join(os.path.dirname(data_path), f'images/forces/full/{component}'), exist_ok=True)
    
    if component == 'F':
        figure_extension = '_thrust.png'
        y_label = 'Thrust $[N]$'

    elif component == 'M':
        figure_extension = '_torque.png'
        y_label = 'Torque $[Nm]$' 

    else:
        print('Not a valid option for component, use either F for Thrust or M for Torque')

    if format == 'integrated':

        data = pd.read_csv(os.path.join(data_path,file_name),skiprows=skip,delimiter='\s+')
        case_name = file_name.split('.')[0] + f'_{component}{axis}'
        force = component + axis
        time = data['time(s)']
        force_data = (data[force])
        rev_time = 60/rpm

        start_time_i = next(x for x, val in enumerate(time) if val > rev_time * transient_rev)
        end_time_i = next(y for y, val in enumerate(time) if val > rev_time * (transient_rev + n_rev))
        filter_force = force_data[start_time_i:end_time_i]
        filter_time = time[start_time_i:end_time_i]


        if case == 'Plot':

            plt.plot(time,force_data,'k')
            plt.grid()
            plt.xlabel('Time $[s]$',fontsize=14, style='italic')
            plt.ylabel(y_label,fontsize=14, style='italic')
            #plt.title('Thrust of the system')
            plt.xticks(fontsize=14)
            plt.yticks(fontsize=14)
            #ax=plt.gca()
            #ax.set_yticklabels([])
            #ax.yaxis.tick_right()
            plt.tight_layout()
            plt.savefig(os.path.join(os.path.dirname(data_path), f'images/forces/full/{component}', case_name + figure_extension), dpi=600)
            #plt.show()
            plt.close()

            plt.plot(filter_time,filter_force,'k')
            plt.grid()
            plt.xlabel('Time $[s]$',fontsize=14, style='italic')
            plt.ylabel(y_label,fontsize=14, style='italic')
            #plt.title('Filtered Thrust of the system')
            plt.xticks(fontsize=14)
            plt.yticks(fontsize=14)
            #ax=plt.gca()
            #ax.set_yticklabels([])
            #ax.yaxis.tick_right()
            plt.tight_layout()
            plt.savefig(os.path.join(os.path.dirname(data_path), f'images/forces/filtered/{component}', case_name + f'_{n_rev}revn_filtered' + figure_extension), dpi=600)
            #plt.show()
            plt.close()

        elif case == 'Data':
            pass

        else:
            print('Not a valid case (Data or Plot)')

        print('The mean', y_label, 'after the transient data is: ', np.mean(filter_force))
        
    elif format == 'strips':
        
        filename = f'forces_component_F{axis}.txt'
        data = pd.read_csv(os.path.join(data_path,filename),delimiter=',')
        
        time = data['time(s)']
        plt.figure()
        
        for i in range(1,data.shape[1]):
            plt.plot(time,data.iloc[:,i],label=f'Strip {i}')
        
        plt.grid()
        plt.xlabel('Time $[s]$',fontsize=14, style='italic')
        plt.ylabel(y_label,fontsize=14, style='italic')
        #plt.title('Filtered Thrust of the system')
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        #ax=plt.gca()
        #ax.set_yticklabels([])
        #ax.yaxis.tick_right()
        plt.tight_layout()
        plt.savefig(os.path.join(data_path, f'images/forces/filtered/{component}', filename.split('.')[0] + '_Strips' + figure_extension), dpi=600)
        #plt.show()
        plt.close()
        

    return(time,force_data,filter_time,filter_force)

def polar_forces_plot(rpm:int,data_path:str,file_name:str,case:str,transient_rev=5,n_rev=2,skip=18):
    
    '''
    This function is used to calculate the polar plot of the forces from PowerFlow simulations.\n
    The input parameters need to be a .txt file output from the exaritool method.\n
    This function will output three different images associated to Fx(Fy), Fx(Fz) and Fz(Fy)\n
    \n
    - rpm = the rotational speed of the system in revolutions per minute\n
    - data_path = with the path where the data is stored\n
    - image_folder = with the path where the images will be stored\n
    - file_name = with the name of the file that wants to be analyzed\n
    - case = either Plot or Data depending if the data wants to be ploted or not\n
    - transient_rev = the starting number where the filterring wants to be done. If not specified 5rev will be the lower limit\n
    - n_rev = the number of revolutions to consider for the filtering after the transient. If not specified the setup value is 5
    - skip = the number of lines to skip in the .txt file. 18 is used as default\n
    '''
    
    os.makedirs(os.path.join(os.path.dirname(data_path), 'images/polar_forces'), exist_ok=True)
    
    data = pd.read_csv(os.path.join(data_path,file_name),skiprows=skip,delimiter='\s+')
    case_name = file_name.split('.')[0]
    forces_components = ['Fx','Fy','Fz']
    combination_list = []
    
    for i in itertools.combinations(forces_components,2):
        combination_list.append(i)
    
    time = data['time(s)']
    globals()['force_Fx'] = (data[forces_components[0]])
    globals()['force_Fy'] = (data[forces_components[1]])
    globals()['force_Fz'] = (data[forces_components[2]])
    rev_time = 60/rpm
    
    time_filter, globals()['force_Fx_filter'], index = data_filtering(time,globals()['force_Fx'],rev_time*transient_rev,rev_time * (transient_rev + n_rev))
    time_filter, globals()['force_Fy_filter'], index = data_filtering(time,globals()['force_Fy'],rev_time*transient_rev,rev_time * (transient_rev + n_rev))
    time_filter, globals()['force_Fz_filter'], index = data_filtering(time,globals()['force_Fz'],rev_time*transient_rev,rev_time * (transient_rev + n_rev))
    
    if case == 'Plot':
        
        for i in combination_list:
            
            force_1 = globals()[f'force_{i[0]}']
            force_2 = globals()[f'force_{i[1]}']
            filter_force_1 = globals()[f'force_{i[0]}_filter']
            filter_force_2 = globals()[f'force_{i[1]}_filter']
            
            plt.plot(force_1,force_2,'k')
            plt.grid()
            plt.xlabel(f'{i[0]} $[N]$',fontsize=14, style='italic')
            plt.ylabel(f'{i[1]} $[N]$',fontsize=14, style='italic')
            #plt.title('Thrust of the system')
            plt.xticks(fontsize=14)
            plt.yticks(fontsize=14)
            #ax=plt.gca()
            #ax.set_yticklabels([])
            #ax.yaxis.tick_right()
            plt.tight_layout()
            print('Saving non filtered polar forces')
            print(40*'-')
            plt.savefig(os.path.join(os.path.dirname(data_path), 'images/polar_forces', case_name + f'-{i[0]}-{i[1]}.png'), dpi=600)
            #plt.show()
            plt.close()

            plt.plot(filter_force_1,filter_force_2,linestyle='-', color='k')
            plt.grid()
            plt.xlabel(f'{i[0]} $[N]$',fontsize=14, style='italic')
            plt.ylabel(f'{i[1]} $[N]$',fontsize=14, style='italic')
            #plt.title('Filtered Thrust of the system')
            plt.xticks(fontsize=14)
            plt.yticks(fontsize=14)
            # plt.ylim([1.25,3.25])
            # plt.xlim([-4.75,-1.75])
            ax=plt.gca()
            ax.set_yticklabels([])
            ax.set_xticklabels([])
            #ax.yaxis.tick_right()
            plt.tight_layout()
            print('Saving filtered polar forces')
            print(40*'-')
            plt.savefig(os.path.join(os.path.dirname(data_path), 'images/polar_forces', case_name + f'-{i[0]}-{i[1]}-{n_rev}revn_filtered.png'), dpi=600)
            #plt.show()
            plt.close()
    
    return

def induced_velocity_fluidline(root_c:float,tip_c:float,n_strips:int,filename:str,data_path:str,image_path:str,units_p='m',units_v='m/sec'):
    
    '''
    Function to evaluate the induced velocity from a fluid line of PF
    '''

    data = pd.read_csv(os.path.join(data_path,filename),delimiter=',')
    position = pd.to_numeric(data[f'Position[Length:{units_p}]'], errors='coerce')
    velocity = pd.to_numeric(data[f'Value[Velocity:{units_v}]'], errors='coerce')
    
    start_i = np.where((position > root_c) & (position < tip_c))[0][0]
    finish_i = np.where((position > root_c) & (position < tip_c))[0][-1]
    
    position = position[start_i:finish_i]
    velocity = velocity[start_i:finish_i]
    strip_step = int(np.floor(position.size/n_strips))
    
    position_avg = np.zeros(n_strips)
    velocity_avg = np.zeros(n_strips)
    
    with open(os.path.join(data_path,filename.split('.')[0]+f'_{n_strips}strip_avg.txt'),'w') as file:
        file.write(f'position[{units_p}] velocity{units_v}\n')
        for i in range(1,n_strips+1):
            
            position_avg[i-1] = np.mean(position.loc[start_i+(i-1)*strip_step:start_i+i*strip_step])
            velocity_avg[i-1] = np.mean(velocity.loc[start_i+(i-1)*strip_step:start_i*i*strip_step])
            
            file.write(f'{position_avg[i-1]} {velocity_avg[i-1]}\n')
                    
    plt.figure()
    plt.plot(position,velocity,'k')
    plt.grid()
    plt.xlabel(f'Position $[{units_p}]$',fontsize=14, style='italic')
    plt.ylabel(f'Axial Velocity $[{units_v}]$',fontsize=14, style='italic')
    #plt.title('Thrust of the system')
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    #ax=plt.gca()
    #ax.set_yticklabels([])
    #ax.yaxis.tick_right()
    plt.tight_layout()
    print('-------Saving non-avg figure---------\n')
    plt.savefig(os.path.join(image_path, 'exports/induced_velocity/non-avg', filename.split('.')[0] + '_non-avg.png'), dpi=600)
    plt.close()

    plt.figure()
    plt.plot(position_avg,velocity_avg,color='k',marker='o')
    plt.grid()
    plt.xlabel(f'Position $[{units_p}]$',fontsize=14, style='italic')
    plt.ylabel(f'Axial Velocity $[{units_v}]$',fontsize=14, style='italic')
    #plt.title('Thrust of the system')
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    #ax=plt.gca()
    #ax.set_yticklabels([])
    #ax.yaxis.tick_right()
    plt.tight_layout()
    print('-------Saving avg figure---------\n')
    plt.savefig(os.path.join(image_path, 'exports/induced_velocity/avg', filename.split('.')[0] + '_avg.png'), dpi=600)
    plt.close()

    

def psnc_surface_avg(rpm:int,first_it:list,second_it:list,variable:str,data_path:str,image_path:str,case:str,transient_rev=10,n_rev=2,skiprows=12):

    '''

    Avg value one per time for all the different locations

    '''

    base_filename = 'probe_in_duct_{-}pr_{-}deg_export.txt' # To compare with input line 121
    print('The base file name is as:', base_filename)
    figure_extension = '.png'
    rev_time = 60/rpm

    first_file = 'probe_in_duct_0pr_0deg_export.txt' # Input to change
    image_name = variable + '-mean-' + first_file.split('_')[0] + '-' + first_file.split('_')[1] + '-' + first_file.split('_')[2]
    first_data = pd.read_csv(os.path.join(data_path,first_file),skiprows=skiprows,delimiter='\s+',header=None)
    
    time = first_data[0]
    mean_data = np.zeros(time.size)
    time_filter, data_filter, index = data_filtering(time,mean_data,rev_time*transient_rev,rev_time * (transient_rev + n_rev))
    mean_data = np.zeros(time_filter.size)
    
    z_iter = 0
    for z in range(index[0],index[1]):
        # pdb.set_trace()
        inst_data = np.zeros(len(first_it) * len(second_it))
        i_iter = 0
        
        for i in first_it:
            # pdb.set_trace()
            for k in second_it:

                file_name = f'probe_in_duct_{i}pr_{k}deg_export.txt' # Input to change
                data_full = pd.read_csv(os.path.join(data_path,file_name),skiprows=skiprows,delimiter='\s+',header=None)
                
                # print(f'doing variable selection for 1st it {i} and 2nd it {k}')
                data,y_label = variable_selection(variable,data_full)
                
                inst_data[i_iter] = data[z]
                i_iter += 1
        
        print(f'doing variable selection for timsetep {z}/{index[1]}')
        mean_data[z_iter] = np.mean(inst_data)
        z_iter += 1

    with open(os.path.join(data_path,image_name+f'-{n_rev}rev-data.txt'),'w') as file:
        for i,j in zip(time_filter,mean_data):
            file.write(f'{i} {j}\n')

    # time_filter, data_filter = data_filtering(time,mean_data,rev_time*transient_rev,rev_time * (transient_rev + n_rev))

    # pdb.set_trace()
    if case == 'Plot':

        print('--- Generating image for', image_name, '---')

        # plt.plot(time,mean_data,'k')
        # # plt.axvline(initialization_t, color = 'r')
        # plt.grid()
        # plt.xlabel('Time $[s]$',fontsize=14, style='italic')
        # plt.ylabel(y_label,fontsize=14, style='italic')
        # plt.xticks(fontsize=14)
        # plt.yticks(fontsize=14)
        # #plt.title('Time trace of ' + title)
        # plt.tight_layout()
        # plt.savefig(os.path.join(image_path, 'timetrace/probes/vanes/impigning/full', image_name + figure_extension), dpi=600)
        # #plt.show()
        # plt.close()

        plt.plot(time_filter,mean_data,'k')
        # plt.axvline(initialization_t, color = 'r')
        plt.grid()
        plt.xlabel('Time $[s]$',fontsize=14, style='italic')
        plt.ylabel(y_label,fontsize=14, style='italic')
        plt.xticks(fontsize=14)
        plt.yticks(fontsize=14)
        #plt.title('Time trace of ' + title)
        plt.tight_layout()
        plt.savefig(os.path.join(image_path, 'timetrace/probes/vanes/impigning/filtered', image_name + f'-{n_rev}rev-filtered' + figure_extension), dpi=600)
        #plt.show()
        plt.close()

    elif case == 'Data':
        pass
    else:
        print('Not a valid case (Data or Plot)')



    return mean_data,time_filter

def check_monotonicity(*parameters):
    diff_parameters = np.diff(parameters)

    if np.all(diff_parameters > 0):
        pass
    else:
        print('The given parameters are not strictly monotonic and the gci analysis may induce wrong results')

def gci(image_path,grid_number,variable_ylabel,case,*parameters):

    check_monotonicity(parameters)
    r = 2 # Value defined due to PowerFlow meshing technique

    if grid_number == 3:
        coarse_parameter = parameters[0]
        medium_parameter = parameters[1]
        fine_parameter = parameters[2]

        parameters_array = np.array([fine_parameter,medium_parameter,coarse_parameter])


        # Order of convergence (p). To check the theoretical value

        p = (np.log((coarse_parameter - medium_parameter)/(medium_parameter - fine_parameter)))/(np.log(r))

        print('---------------------------------------')
        print('The obtained value of p is', np.round(p,decimals=2))
        print('---------------------------------------')

        # Define the value of Fs based on the amount of grids

        Fs = 1.25

        gci_12 = ((Fs*np.abs((fine_parameter-medium_parameter)/fine_parameter))/((r**p) - 1)) * 100
        gci_23 = ((Fs*np.abs((medium_parameter-coarse_parameter)/medium_parameter))/((r**p) - 1)) * 100

        gci = (gci_23)/((r**p)*gci_12)

        print('---------------------------------------')
        print('The obtained value of GCI is', np.abs(np.round(gci_23,decimals=2)), '%')
        print('The ratio of the two grids is', np.round(gci,decimals=2))
        print('value suggested to be close to 1')
        print('---------------------------------------')

    elif grid_number == 2:
        medium_parameter = parameters[0]
        fine_parameter = parameters[1]

        parameters_array = np.array([fine_parameter,medium_parameter])

        # Due to amount of grid used the theoretical parameter must be used (assumed to be 2)

        p = 2

        # Define the value of Fs based on the amount of grids

        Fs = 3

        gci_23 = ((Fs*np.abs((fine_parameter-medium_parameter)/fine_parameter))/((r**p) - 1)) * 100

        print('---------------------------------------')
        print('The obtained value of GCI is', np.abs(np.round(gci_23,decimals=2)), '%')
        print('---------------------------------------')

    else:
        print('The convergence study does not use the minimum amount of grid (2)')



    # Expected value for the richardson extrapolation (f_0) based on the fine case data

    f_0 = fine_parameter + ((medium_parameter - fine_parameter) * (r**p))/(r**p - 1)
    grid_array = np.arange(1,grid_number+1,dtype=int)


    plt.plot(grid_array,parameters_array,'ko-',label='Simulation parameters')
    plt.plot(0,f_0,'kx',label='Richardson extrapolation parameter')
    plt.xlabel('Normalized grid spacing')
    plt.ylabel(variable_ylabel)
    plt.legend()
    #plt.title('Grid convergence study')
    ax=plt.gca()
    ax.yaxis.set_major_formatter(ticker.StrMethodFormatter("{x:.2f}"))
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.savefig(image_path + case + '.png', dpi=1200)
    #plt.show()
    plt.close()

    #pdb.set_trace()

    print('---------------------------------------')
    print('The value for the', variable_ylabel, 'on the fine grid is')
    print(fine_parameter, 'and expected to be', np.round(f_0,decimals=2))
    print('with an error band of', np.abs(np.round(gci_23,decimals=2)), '%')
    print('assuming a zero grid spacing')
    print('---------------------------------------')

def ospl(freq,pxx,*f_lim):

    if np.size(f_lim) == 0:

        ospl = 10*np.log10(np.trapz(pxx,freq)/4.0e-10)

    else:

        freq_f = freq[np.where((freq>f_lim[0])&(freq<f_lim[1]))[0]]
        pxx_f = pxx[np.where((freq>f_lim[0])&(freq<f_lim[1]))[0]]
        ospl = 10*np.log10(np.trapz(pxx_f,freq_f)/4.0e-10)

    print('The obtained OSPL value is ',np.round(ospl,decimals=2),'dB')
    return(np.round(ospl,decimals=2))

def ospl_correction(freq,pxx,correction,*f_lim):

    if np.size(f_lim) == 0:

        ospl = 10*np.log10(np.trapz(pxx,freq))+correction

    else:

        freq_f = freq[np.where((freq>f_lim[0])&(freq<f_lim[1]))]
        pxx_f = pxx[np.where((freq>f_lim[0])&(freq<f_lim[1]))]
        ospl = 10*np.log10(np.trapz(pxx_f,freq_f))+correction

    print('The obtained OSPL value is ',np.round(ospl,decimals=2),'dB')
    return(np.round(ospl,decimals=2))

def extract_data(freq,pxx,target,case):

    if case == 'exact':
        freq_f = freq[np.where(freq>=target)[0][0]]
        pxx_f = pxx[np.where(freq>=target)[0][0]]
    elif case == 'average':
        freq_f = np.average([freq[np.where(freq>target)[0][0]],freq[[np.where(freq<target)[0][-1]]]])
        pxx_f = np.average([pxx[np.where(freq>target)[0][0]],pxx[np.where(freq<target)[0][-1]]])


    print(np.round(freq_f,decimals=2),'Hz is the closest frequency to target')

    print('The value at the mentioned frequency is', np.round(10*np.log10(pxx_f/4.0e-10),decimals=2),'dB')
    return(np.round(10*np.log10(pxx_f/4.0e-10)))

def directivity(variable_name,image_path,angles,cases,legend,*parameters):

    max_val = np.zeros(cases)

    for i in range(cases):
        print('parameter',i)
        plt.polar(angles*np.pi/180,np.transpose(parameters[i]),'o-')
        max_val[i] = np.max(parameters[i])
        i += 1

    plt.grid(True, which='both', ls='--')
    plt.legend(legend,ncol=2)
    plt.ylim([0,1.05*np.max(max_val)])
    ax=plt.gca()
    ax.set_yticklabels([])
    ax.yaxis.tick_right()
    plt.savefig(image_path + 'directivity_' + variable_name + 'max=' + np.str(np.round(np.max(max_val),decimals=2)) + '.png' , dpi=600)
    #plt.show()
    plt.close()
    return()

def spectrogram(n_blades,rpm,first_it,second_it,data_path,image_path,type_probe,transient_rev,skiprows=12):

    #type_probe = str(input('Enter if you want to perform the analysis of a fluid, surface or experimental probe\n'))
    base_filename = 'probe_in_duct_{-}pr_{-}z_no_impigning_side_export.txt' # To compare with input line 121
    print('The base file name is as:', base_filename)

    BPF = n_blades * rpm / 60
    rev_time = 60/rpm

    for i in first_it:
        for j in second_it:

            if type_probe == 'fluid':
                
                file_name = f'mic{i}_{j}.txt'
                image_name = file_name.split('.')[0]
                data_full = pd.read_csv(os.path.join(data_path,file_name),delimiter='\s+')
                time = data_full['#time']
                data = data_full['static_pressure']

            elif type_probe == 'surface':

                file_name = f'probe_in_vanes_{i}pr_{j}z_impigning_side_export.txt' # Input to change
                image_name = file_name.split('.')[0]
                data_full = pd.read_csv(os.path.join(data_path,file_name),skiprows=skiprows,delimiter='\s+',header=None)
                time = data_full[0]
                data = data_full[1]

            elif type_probe == 'experimental':

                warnings.warn('The value entered has not yet been implemented in this code version')
                sys.exit('Program terminated due to invalid warning')
            else:

                raise ValueError('Not a valid type of probe. Code terminated')

            start_time_i = next(x for x, val in enumerate(time) if val > rev_time * transient_rev)
            time = time[start_time_i:]
            data = data[start_time_i:]

            data = data - np.mean(data)

            n_chunk=4
            lensg=data.size
            nperseg=lensg/n_chunk
            dt = time[start_time_i+1]-time[start_time_i]
            fs = 1.0/dt
            nfft=next_greater_power_of_2(int(nperseg))

            f, t, pxx = sg.spectrogram(data,fs=fs,window='hann',nperseg=int(nperseg),nfft=nfft,scaling='density')

            cmap = plt.get_cmap('hot_r')
            fig = plt.pcolormesh(t, f/BPF, 10*np.log10(pxx/4.0e-10), shading='gouraud',cmap=cmap,vmin=10,vmax=40)
            plt.ylim([100/BPF,np.max(f)/BPF])
            ax = plt.gca()
            ax.set_yscale('log')
            ax.yaxis.set_major_formatter(FormatStrFormatter('%1.0f'))
            ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
            plt.ylabel('BPF $[-]$',fontsize=14, style='italic')
            plt.xlabel('Time $(s)$',fontsize=14, style='italic')
            cbar = plt.colorbar()
            cbar.set_label('PSD $[dB/Hz]$',fontsize=14, fontstyle='italic')
            plt.xticks(fontsize=14)
            plt.yticks(fontsize=14)
            plt.tight_layout()
            print(f'Saving spectrogram for data with iterating {i}-{j}')
            plt.savefig(os.path.join(image_path, 'spectrograms',image_name + '_spectrogram-2.png'), dpi=600)
            plt.close()


    return()

def coherence(nblades,speed,image_path,variable,path_0,name_0,path_1,name_1):

    BPF = nblades * speed / 60

    data_0 = pd.read_csv(path_0+name_0,delimiter='\s+')
    time_0 = data_0['#time']
    pressure_0 = data_0['static_pressure']

    pressure_0 = pressure_0 - np.mean(pressure_0)

    data_1 = pd.read_csv(path_1+name_1,delimiter='\s+')
    time_1 = data_1['#time']
    pressure_1 = data_1['static_pressure']

    pressure_1 = pressure_1 - np.mean(pressure_1)

    dt = time_0[1]-time_0[0]

    if len(time_0) != len(time_1):
        print('WARNING: A cut of data need to be done due to different array size')
        if len(time_0) > len(time_1):
            time_0 = time_0[-len(time_1):]
            pressure_0 = pressure_0[-len(time_0):]
        elif len(time_1) > len(time_0):
            time_1 = time_1[-len(time_0):]
            pressure_1 = pressure_1[-len(time_0):]


    n_chunk=4
    lensg=pressure_0.size
    nperseg=lensg/n_chunk
    fs = 1.0/dt
    nfft=next_greater_power_of_2(int(nperseg))

    [f, cor] = sg.coherence(pressure_0,pressure_1, fs=fs, window='hann', nperseg=nperseg, nfft=nfft, detrend='constant')

    plt.figure()
    plt.plot(f/BPF,cor,'k--',)
    plt.grid(True, which='both', ls='--')
    plt.xlabel('Frequency $\\left[Hz/BPF\\right]$',fontsize=12, style='italic')
    plt.ylabel('Coherence $\\left[-\\right]$',fontsize=12, style='italic')
    plt.xlim([100/BPF,np.max(f/BPF)])
    #plt.title('The mean coherence is' + str(np.mean(cor)))
    ax=plt.gca()
    ax.set_xscale('log')
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.1f'))
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.savefig(image_path + variable +'_coherence.png', dpi=600)
    #plt.show()

    return(f,cor)


