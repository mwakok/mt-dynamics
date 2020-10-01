import os
import sys
import numpy as np
from parameters import ParameterSet
import pytest

path_root = os.path.dirname(os.getcwd())
sys.path.insert(0, os.path.join(path_root, "mtdynamics"))
print(os.path.join(path_root, "mtdynamics"))

#from . import simulation_main as simulation_main
from mtdynamics.simulation_main import mt_run
import mtdynamics.simulation_parameters as simParameters #import simulation parameters


def test_basic_run():
    np.random.seed(1000)
    
    # Load parameters from simulation_parameters.py
    simPa = ParameterSet(simParameters.simParameters)
    
    # Change some parameters
    number_of_catastrophes = 5
    simPa.no_cat = number_of_catastrophes
    simPa.D_tip = 3000
    simPa.plot_figures = False
    simPa.record_data = False
    
    dt, MT_length_sum, MT_length_full, CATASTROPHE_TIMES, CATASTROPHE_LENGTH, \
    barrier_contact_times, EB_comet_sum, MT_cap, Cap_threshold, frame_rate_actual, \
    EB_profiles, washout_times, catastrophe_washout, Cap_length_sum  = mt_run(simPa)
    
    assert dt == pytest.approx(0.2461538, 1e-6), "Expected different dt"
    assert frame_rate_actual == pytest.approx(0.2461538, 1e-6), "Expected different frame rate"
    assert len(MT_cap) == number_of_catastrophes
    assert barrier_contact_times.shape == (0,), "Expected no barrier contact times"
