import os
import sys
import numpy as np
from parameters import ParameterSet
import pytest

path_root = os.path.dirname(os.getcwd())
sys.path.insert(0, os.path.join(path_root, "code"))

import simulation_functions as sFUNC #import simulation functions
import simulation_parameters as simParameters #import simulation parameters
import simulation_main as sMAIN


def test_basic_run():
    np.random.seed(1000)
    
    # Load parameters from simulation_parameters.py
    simPa = ParameterSet(simParameters.simParameters)
    
    # Change some parameters
    number_of_catastrophes = 5
    simPa.no_cat = number_of_catastrophes
    simPa.D_tip = 3000
    
    dt, MT_length_sum, MT_length_full, CATASTROPHE_TIMES, CATASTROPHE_LENGTH, \
    barrier_contact_times, EB_comet_sum, MT_cap, Cap_threshold, frame_rate_actual, \
    EB_profiles, washout_times, catastrophe_washout, Cap_length_sum  = sMAIN.MT_RUN(simPa)
    
    assert dt == pytest.approx(0.2461538, 1e-6), "Expected different dt"
    assert frame_rate_actual == pytest.approx(0.2461538, 1e-6), "Expected different frame rate"
    assert len(MT_cap) == number_of_catastrophes
    assert barrier_contact_times.shape == (0,), "Expected no barrier contact times"
