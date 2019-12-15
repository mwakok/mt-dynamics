#
# mt-dynamics
#
# Copyright 2019 Florian Huber
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

import simulation_parameters as simpa
from parameters import ParameterSet
from itertools import islice
import csv
import json
import numpy as np


def dict_to_json(mydict, file_json): 
    # Save dictionary as json file
    with open(file_json, 'w') as outfile:  
        json.dump(mydict, outfile)
        

def json_to_dict(file_json): 
    # Create dictionary from json file
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


def load_parameters(simpa,
                    growth_speed):
    """ Add all parameters to one dictionary
    """

    simParameters = simpa.simParameters
    simParameters['growth_rate_one'] = growth_speed/(60*simpa.dL_dimer) #rate of one dimer per s
    return simParameters
    
    
def analyse_EB_signal(simPa, EB_comet_sum, barrier_contact_times):
    # TODO: function redundant with function in plotting_functions.py
    
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


def analyse_EB_profile(simPa, 
                       MT_length_full, 
                       EB_profiles, 
                       dt, 
                       w_size):
    # TODO: function redundant with function in plotting_functions.py
    """ Calculate the mean GTP/GDP-Pi profile at the microtubule end during steady-state growth
    
    Args:
    -------
    
    """
       
    # Initialize arays
    v_mean = []
    EB_mean = []
    
    # Microscope parameters
    resolution = 0.25 #um
    
    # Number of dimers in PSF
    resolution_dimers = int(np.ceil(resolution/simPa.dL_dimer))  
    
    # Loop over number of simulated microtubules
    for num_run in range(len(MT_length_full)):
        
        # Obtain the time and length arrays
        time = np.arange(0, len(MT_length_full[num_run]), 1) * dt
        MT_length = (np.asarray(MT_length_full[num_run]) - MT_length_full[num_run][0]) * simPa.dL_dimer *1000
            
        # Find the local mean growth speeds in order to exclude pausing state from the profile analysis    
        if len(MT_length) > w_size:
            di = 0
            v_fit = np.zeros(len(MT_length) - w_size + 1)
            for i in window(MT_length, w_size):
                v_fit[di] = np.polyfit(np.linspace(0, dt*w_size-1, w_size), i, 1)[0]
                di = di + 1
        else: v_fit = [] 
        
        # Calculate mean growth speed
        v = np.polyfit(time, MT_length, 1)
        v_mean.append(v[0])
        v_thres = 0.6
        
        # Identify steady-state growth events in the trace
        matches = [i for i, x in enumerate(v_fit) if x > v_thres*v_mean[-1]]
        matches = np.asarray(matches) + w_size//2
        
        if matches.size > 0:
            for mm in matches:                
                # Extend the EB profile array 
                EB_new = np.append(EB_profiles[num_run][mm], np.zeros(resolution_dimers))
                EB_mean.append(EB_new)
    
    
    return EB_mean


# -----------------------------------------------------------------------------
# ----------------------------- Helper functions ------------------------------
# -----------------------------------------------------------------------------
    

def gaussian(x, mu, sig):
    # Define Gaussian distribution
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))


def window(seq, n):
    """ # Define a sliding window.
    Returns a sliding window (of width n) over data from the iterable
    s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...
    """   
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result