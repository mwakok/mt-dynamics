#
# mt-dynamics
#
# Copyright 2019 Florian Huber, Maurits Kok
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