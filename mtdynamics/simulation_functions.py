import mtdynamics.simulation_parameters as simpa
from parameters import ParameterSet
from itertools import islice
import csv
import json
import numpy as np


def dict_to_json(mydict: dict, file_json: str):
    """ Save dictionary as json file."""
    with open(file_json, 'w') as outfile:
        json.dump(mydict, outfile)


def json_to_dict(file_json: str):
    """ Create dictionary from json file."""
    with open(file_json) as infile:
        mydict = json.load(infile)

    return mydict


def dict_to_csv(mydict: dict, file_csv: str):
    """ Writes dictionary to csvfile."""
    csvfile = csv.writer(open(file_csv, 'w', newline=''))#, quoting=csv.QUOTE_NONNUMERIC)
    for key in mydict.keys():
        csvfile.writerow((key,mydict[key], type(mydict[key])))


def csv_to_dict(file_csv: str):
    """ Load csv file and create dictionary."""
    mydict = {}
    reader = csv.reader(open(file_csv, 'r', newline=''))#, quoting=csv.QUOTE_NONNUMERIC)
    for rows in reader:
        print(rows)
        #if row[2] ==
        mydict[rows[0]] = rows[1]

    return mydict


def load_parameters(simpa, growth_speed: float):
    """ Add all parameters to one dictionary
    """

    simParameters = simpa.simParameters
    simParameters['growth_rate_one'] = growth_speed/(60*simpa.dL_dimer) #rate of one dimer per s
    return simParameters


# -----------------------------------------------------------------------------
# ----------------------------- Helper functions ------------------------------
# -----------------------------------------------------------------------------

def gaussian(x: float, mu: float, sig: float):
    """ Define Gaussian distribution."""
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))


def window(seq, n):
    """ Define a sliding window.
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
