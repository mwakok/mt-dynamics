"""
Function to run and initialize the MT simulation

@author: FlorianHuber
"""

import simulation_parameters as simpa
from parameters import ParameterSet
import csv
import json
 


def dict_to_json(mydict, file_json): 
    # save dictionary as json file
    with open(file_json, 'w') as outfile:  
        json.dump(mydict, outfile)
        

def json_to_dict(file_json): 
    # create dictionary from json file
    with open(file_json) as infile:  
        mydict = json.load(infile)

    return mydict
    


def dict_to_csv(mydict, file_csv):
    # Writes dictionary to csvfile
    csvfile = csv.writer(open(file_csv, 'w', newline=''))#, quoting=csv.QUOTE_NONNUMERIC)
    for key in mydict.keys():
        csvfile.writerow((key,mydict[key], type(mydict[key])))



def csv_to_dict(file_csv):
    # Load csv file and create dictionary 
    mydict = {}
    reader = csv.reader(open(file_csv, 'r', newline=''))#, quoting=csv.QUOTE_NONNUMERIC)
    for rows in reader:
        print(rows)
        #if row[2] ==
        mydict[rows[0]] = rows[1]
            
    return mydict


def load_parameters():
    # add all parameters to one dictionary

#    if np.where(simpa.EB_growth_speed[:,0] == 0.2)[0].shape[0] == 1:
#        growth_speed = simpa.EB_growth_speed[np.where(simpa.EB_growth_speed[:,0] == simpa.EB),1]
#    else:
#        print('No growth speed found for given EB concentration')
#    sPARAM['growth_speed'] = growth_speed
#    print('Growth speed set to: ', growth_speed)

    simParameters = simpa.simParameters
    simParameters['growth_rate_one'] = growth_speed/(60*simpa.dL_dimer) #rate of one dimer per s
    
    
def analyse_EB_signal(simPa, EB_comet_sum, barrier_contact_times):
    
    # Select valid runs
    b = []
    for a in range(0,len(EB_comet_sum)):
        b.append(len(EB_comet_sum[a]) * simPa.frame_rate_actual - barrier_contact_times[a]) 
    
    valid_runs = np.where(np.array(b) > simPa.min_length_run)[0]  
    
    max_barrier_contact_frames = int(round(np.max(barrier_contact_times/simPa.frame_rate_actual),0)) 
    min_length_run_frames = int(min_length_run/frame_rate_actual)
    frame_window = min_length_run_frames + max_barrier_contact_frames
     
    EB_signal = np.zeros((len(valid_runs), frame_window+1)) #min_length_run+1+max_barrier_contact)) #put individual runs into one np.array
    normalize_EB_signal = np.zeros(frame_window+1) #min_length_run+1+max_barrier_contact)
    for a in range(0,len(valid_runs)):
        #frame_barrier_contact = int(round(barrier_contact_times[valid_runs[a]]/frame_rate_actual,0))
        frame_barrier_contact = int(np.round(barrier_contact_times[valid_runs[a]]/simPa.frame_rate_actual,0))
        EB_signal[a][(max_barrier_contact_frames-frame_barrier_contact):frame_window] \
        = np.array(EB_comet_sum[valid_runs[a]])[0:(min_length_run_frames+frame_barrier_contact)]
        normalize_EB_signal[(max_barrier_contact_frames-frame_barrier_contact):frame_window] +=1
        
    EB_signal_average = np.sum(EB_signal, axis=0)
    EB_signal_average = EB_signal_average/normalize_EB_signal
    
    return valid_runs, EB_signal, EB_signal_average