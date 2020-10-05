## Dictionary to store all parameters
simParameters = {}

## Simulation type choices
simParameters['record_data'] = True
simParameters['record_data_full'] = False
simParameters['plot_figures'] = True
simParameters['compare_results'] = False
simParameters['end_hydrolysis'] = True
simParameters['frame_rate_aim'] = float(0.25) # in seconds
simParameters['record_cap'] = False
simParameters['record_length'] = False
simParameters['steady_state_analysis'] = False

# Key simulation parameters
simParameters['EB'] = 0
#EB_growth_speed = np.array([[0, 1.68],[0.02, 2.79],[0.05, 2.79],[0.1, 3.36]])
if simParameters['EB'] == 0: #set growth speed (only used for our special case... remove later)    
    kBC = 0.1
    D_tip = 2000   
    growth_speed = 1.8
elif simParameters['EB'] == 0.02:
    growth_speed = 2.79 
    kBC = 0.20
    D_tip = 3600  
elif simParameters['EB'] == 0.05:
    growth_speed = 2.79 
    kBC = 0.24
    D_tip = 3420      
elif simParameters['EB'] == 0.1:
    kBC = 0.39
    D_tip = 3100
    #growth_speed = 3.36 #µm/min
    growth_speed = 3.72
print('growth speed set to: ', growth_speed)
    
simParameters['growth_speed'] = growth_speed
#simParameters['growth_speed'] = 1.5 #µm/min
simParameters['dL_dimer'] = 0.008/13 # Eukaryotic tubulin dimer length in um
#simParameters['dL_dimer'] = 0.008/5 # Prokaryotic tubulin dimer length in um
simParameters['growth_rate_one'] = growth_speed/(60*simParameters['dL_dimer']) #rate of one dimer per s
## Best values for kBC and D_tip:
# 0 nM EB:      0.08 and 2900 OR 0.07 and 3720 || 0.07 and 3720
# 20 nM EB:     0.18 and 3500 OR 0.19 and 3580 || 0.20 and 3400
# 50 nM EB:     0.28 and 3000 OR 0.23 and 3840 || 0.24 and 3420
# 100 nM EB:    0.50 and 1700 OR 0.53 and 2100 || 0.39 and 3100

N_unstable = int(15);

simParameters['kBC'] = kBC #s^-1 
simParameters['D_tip'] = D_tip #tip diffusion nm^2 /s


## tip noise parameters (relative noise if different when in A, B, C state...)
#noise_STD_A = float(1) #currently not used anyore --> TODO: remove from simulation!!
#noise_STD_B = float(1)
#noise_STD_C = float(1)
#noise_STD_seed = float(0) #no noise within seed
#simParameters['tip_noise_relative'] = [noise_STD_A, 
#      noise_STD_B, noise_STD_C, noise_STD_C, noise_STD_seed]


simParameters['unstable_cap_criteria'] = N_unstable 
simParameters['seed_resolution'] = 0.005 #spatial resolution at seed (to distinguish catastrophes)

# Parameters for time-dependent cap criterium
simParameters['unstable_cap_time'] = False # Enable to simulate a time-dependent Cap threshold
simParameters['unstable_cap_start'] = N_unstable*2  # N_unstable*1.5
simParameters['unstable_cap_end'] = N_unstable      # N_unstable*(1/1.5)
simParameters['unstable_cap_rate'] = 0.007

# Parameters for time-dependent tip diffusion
simParameters['D_tip_time'] = False # Enable to simulate a time-dependent D_tip
simParameters['D_tip_start'] =  D_tip*(2/3)    # 1) D_tip*(1/2), 2)
simParameters['D_tip_end'] = D_tip            # 1) D_tip*(5/4), 2)
simParameters['D_tip_rate_T'] = 0.02   #1) 0.025

# Parameters for length-dependent tip diffusion
simParameters['D_tip_length'] = False # Enable to simulate a length-dependent D_tip
simParameters['D_tip_start'] =  D_tip*(1/3)    # 1) D_tip*(1/2), 2) 0
simParameters['D_tip_end'] = D_tip*(5/4)            # 1) D_tip*(5/4), 2)
simParameters['D_tip_rate_L'] = 0.2 # rate per added dimer

# Parameters for time-dependent hydrolysis rate
simParameters['kBC_time'] = False # Enable to simulate a time-dependent D_tip
simParameters['kBC_start'] = kBC*(1/2)      # kBC*(1/1.5)
simParameters['kBC_end'] = kBC              # kBC*1.5
simParameters['kBC_rate'] = 0.005       #1) 0.01

# if barrier, then give distance in dimers, otherwise set to False
simParameters['barrier'] = False #550*13 (= 3.4 um) or False
simParameters['DXbarrier'] = 0.005 #0.01

simParameters['CAP_threshold'] = float(0)

## Parameters for simulation run time and accuracy
# Criteria for early stopping (if run too long...):
simParameters['no_cat'] = 500 #aim: number of catastrophes to simulate (if run not too long)
simParameters['max_MT_growth_time'] = 2000
simParameters['too_stable_check'] = 5 #if more than this #MTs grow longer than max_time
simParameters['total_max_time'] = 500*simParameters['no_cat'] # stop if simulated time > total_max_time

simParameters['P_max'] = 0.05 # maximum probability of hydrolysis event during one time step

simParameters['tip_window'] = 1.0 #1um window for comet

# Nucleation correction
simParameters['nucleation_threshold'] = 400 #400 subunits is equal to 250 nm (~ 2 pixels)

#memory and imaging options:
simParameters['show_fraction'] = 40 #for example figure --> must be <= min_length_run
simParameters['min_length_run'] = 40
simParameters['min_length_begin'] = 100 # duration of growth from seed
simParameters['take_EB_profiles'] = True

simParameters['washout'] = False #True
simParameters['washout_time'] = 160 # in seconds

#screening options:
simParameters['shift_ks_contact'] = 0.01
simParameters['shift_ks_cat'] = 0.01

# EB binding characteristics
simParameters['EB_KD'] = 122 # nM (Roth et al., 2018)
