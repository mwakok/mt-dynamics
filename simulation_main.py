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

# Import packages
import numpy as np

# When switching to Object oriented?
#from parameters import ParameterSet
#
#class MT_run(object):
#    """
#    Main simulation class.
#    MT_run object represents one simulation run
#    """
#    def __init__(self, sim_parameters):
#        self.simulation_done = False
#        self.dt = 0
#        self.MT_length_sum = 0  
#        self.CATASTROPHE_TIMES = 0, 
#        self.CATASTROPHE_LENGTH = 0 
#        self.barrier_contact_times = 0 
#        self.EB_comet_sum = 0 
#        self.cap_end = 0 
#        self.frame_rate_actual = 0 
#        self.EB_profiles = 0
#        
#        # Import parameter set
#        self.simPa = ParameterSet(sim_parameters)
#
#
#    def initialize(self):     

# -----------------------------------------------------------------------------
# ------------ MAIN FUNCTIONS FOR MT SIMULATION ("maurer model") --------------
# -----------------------------------------------------------------------------

def MT_RUN(simPa):
    """
    Main simulation function.
    
    MTs are modeled as a 1D string of tubulin dimers. 
    Those are represented by a 1D numpy array where the different dimer states
    are described by different integers.
    
    5 -- stable seed dimers (no hydrolysis, no disassembly)
    1 -- 
    2 --
    
    
    Args:
    --------
    simPa: parameter set 
        Simulation parameters in "ParameterSet" format.
    """  
    
    # Set time step length:
    dt_max = simPa.P_max/simPa.kBC  # Make sure no transition probability is > P_max
    dt_max = min(simPa.frame_rate_aim, dt_max)  # Do not allow time step to be > than frame rate
    
    dimers_per_step = max(1, int(simPa.growth_rate_one * dt_max))  # See how many dimers can grow in one step
    dt = dimers_per_step / simPa.growth_rate_one
    tip_noise = (2*dt*simPa.D_tip)**0.5/(simPa.dL_dimer*1000)
    
    print('EB: ', simPa.EB, '   kBC: ', simPa.kBC, '   D_tip: ', simPa.D_tip)
    print('--- time step dt = %.3f s' %dt) 
    print('--- MTs growth ', dimers_per_step, ' dimers per step (D_tip = ', simPa.D_tip, ' -> noise STD= ', tip_noise, ' )' )
    frame_taking = round(simPa.frame_rate_aim/dt)
    frame_rate_actual = frame_taking*dt
    
    # Initialize arrays etc.:
    f = np.ones(simPa.unstable_cap_criteria)  # Filter for end-of-cap testing
    
    L_seed = int(np.ceil(simPa.tip_window/simPa.dL_dimer)) #seed 
    L = L_seed
    MT = np.ones(L, dtype=np.int) # Python! so MT[L-1] = MT at position L
    current_cap_end = L_seed
    MT[0:L_seed] = 5  # Set to seed state

    # Calculate different growth rates (NO GOOD EXPERIMENTAL DATA --> =parameters)
    P_growth = min(1, simPa.growth_rate_one*dt) #can be removed????
    
    # Calculate transition rates based on B-C model:   
    P_BC = simPa.kBC * dt
    
    deterministic_growth = np.zeros(1)

    # -------------------------------------------------------------------------
    # ------------------------- Start catastrophe loop ------------------------
    # -------------------------------------------------------------------------
    N_max = int(np.floor(simPa.max_MT_growth_time/dt)) # max number of timesteps
    
    # Initialize variables
    CATASTROPHE_TIMES = []
    CATASTROPHE_LENGTH = []
    catastrophe_at_barrier =[]
    EB_comet_sum = []
    EB_profile_intermediate = np.zeros(L_seed, dtype=np.int)
    EB_profiles = []
    show_fraction = simPa.show_fraction
    show_fraction_frames = int(show_fraction/frame_rate_actual)
    MT_length_sum = np.zeros((simPa.no_cat, show_fraction_frames))
    barrier_contact_times = np.zeros(0);
    
    too_long = 0 # reset
    run_time = 0 # reset
    
    for cat_number in range(0,simPa.no_cat):
        # Reset everything:
        L = L_seed
        Lnew = L
        MT = np.ones(L, dtype=np.int)
        current_cap_end = L_seed
        MT[0:L_seed] = 5 #set state Seed
        
        growth = np.zeros(1)
        #catastrophes = np.zeros(0)
        cap_end = [0] 
        #cap_length = [0] 
        MT_length = [L] 
    
        #barrier_contact_count = 0    
        barrier_contact = 0 #set contact = off:
        
        #washout experiment:
        washout = 0 #reset
        
        #EB monitoring:
        EB_comet = [0]
        EB_profiles_wholerun = []
        
        #interrupt run if too long (e.g. when MT too stable)
        if too_long > simPa.too_stable_check:
            break;
        if run_time > simPa.total_max_time:
            print("run stopped - was longer than time limit!")
            break;
        
        # ---------------------------------------------------------------------
        # ------------------------ Start time loop ----------------------------
        # ---------------------------------------------------------------------
        for i in range(1,N_max+1):  # N_max timesteps
            
            if i == N_max: # Check if MT grows for too long
                run_time += N_max*dt
                too_long += 1
                print("Too long growth event! (", simPa.EB, simPa.kBC, simPa.D_tip, ")")
            
            
            # Tip growth 
            # -----------------------------------------------------------------
            
            deterministic_growth = deterministic_growth + P_growth #here: only to time "growth event"
            
            if deterministic_growth >= 1: #growth event:
                deterministic_growth = deterministic_growth-1
                growth = np.round(dimers_per_step + tip_noise*np.random.randn(1))
                '''if MT[L-1] == 1:
                    growth = round(1+ simPa.noise_STD_A*np.random.randn(1)) #%= random number with std= 1; Round to draw discrete number
                elif MT[L-1] == 4:
                    growth = round(P_growth_C + simPa.noise_STD_C*np.random.randn(1))
                elif MT[L-1] == 2:
                    growth = round(P_growth_BE + simPa.noise_STD_BE*np.random.randn(1))
                else: # for MT[L-1] == 3
                    growth = round(P_growth_B + simPa.noise_STD_B*np.random.randn(1))   '''
                    
                if washout == 1: #ONLY FOR TESTING
                    if growth > 0:
                        growth = 0 #do not allow growth after tubulin washout
                    
                if abs(growth) > 0: #growth process
                    Lnew = L + int(growth)
    #                if hasattr(FACT, 'barrier'): #if variable barrier present
    #                    if Lnew > (barrier - np.ceil(simPa.DXbarrier/simPa.dL_dimer)) and barrier_contact == 0:
    #                        barrier_contact_count += 1
    #                    if barrier_contact_count > 5 and barrier_contact == 0:
    #                        barrier_contact = (i-5)*dt                
                    
                    
                    if Lnew < L:  # MT shrinks
                        if Lnew < L_seed:
                            Lnew = L_seed # +1 to forbid catastrophes at seed
                        else:
                            MT[Lnew:L+1] = 0 
                    elif Lnew > L:  # MT grows
                        
                        #  Change size of MT array when necessary
                        if Lnew > len(MT):  
                            MT = np.append(MT,np.zeros(Lnew-L + 500, dtype=np.int))                
                        
                        if simPa.barrier:  # If variable barrier present  
        
                            if Lnew > (simPa.barrier - np.ceil(simPa.DXbarrier/simPa.dL_dimer)) and barrier_contact == 0: #contact if closer than 0.1um, 0.2µm, 0.3µm?
                                barrier_contact = i*dt 
                                #barrier_contact = (i // frame_taking) * frame_taking * dt #time of first 'contact' with wall                    
                                
                            if Lnew > simPa.barrier:
                                
                                # NO brownian ratchet-like stalling:
                                Lnew = simPa.barrier #still possible to grow, however only till barrier
#                                if barrier_contact == 0: #if first moment of "contact"
#                                    barrier_contact = (i // frame_taking) * frame_taking * dt #time of first 'contact' with wall                    
#                                    #barrier_contact = i * dt #time of first 'contact' with wall   
#                                
                                # Brownian ratchet-like stalling:
                                #F_barrier = simPa.k_barrier*(Lnew - barrier)  # Spring like response at barrier
                                #Lnew = L + (Lnew-L)*(np.random.rand(1)< np.exp(- F_barrier*1000*simPa.dL_dimer/4.1))[0]  # Mogilner: V=Vmax*exp(force*dl/(k_b*T)), k_bT=4.1pN/nm divide by 8nm/dimer
                                if Lnew > L:
                                    MT[L:Lnew] = 2 #start with "B" state                   
                            else:
                                MT[L:Lnew] = 2 #start with "B" state                
                        else:
                            MT[L:Lnew] = 2 #start with "B" state  
                    L = Lnew
                    L = max(L, L_seed) #avoid L < L_seed
                    
                    if simPa.washout_experiment:
                        if barrier_contact > 0 and washout ==0:
                            washout = 1
                    
            # Hydrolysis --> Multi-step reaction ~ Maurer 2014 A->B (->BE) -> C
            # -----------------------------------------------------------------
            
            #old: elements_in_A = np.where(MT == 1)[0] 
            elements_in_B = np.where(MT == 2)[0] 
            #old: from maurer model: elements_in_BE = np.where(MT == 3)[0] 
            #elements_in_C = np.where(MT == 4)[0]   
            
            # now: shift between different stages:           
            # 1- throw dices:
            #old: randomA = np.random.rand(len(elements_in_A))
            randomB = np.random.rand(len(elements_in_B))
            #old: from maurer model:randomBE = np.random.rand(len(elements_in_BE))
            #randomC = np.random.rand(len(elements_in_C))
            
            # 2- change state of elements based on dices:
            #old: MT[elements_in_A[np.where(randomA<P_AB)[0]]] += 1 #go from 1=A to 2=B
            #old: from maurer model: MT[elements_in_A[np.where(randomA<P_ABE)[0]]] += 1 # and from 2=B to 3=BE
            #old: from maurer model: MT[elements_in_A[np.where(randomA>(1-P_AC))[0]]] += 3 #from 1=A to 4=C 
            
            #old: from maurer model: MT[elements_in_B[np.where(randomB<P_BBE)[0]]] += 1 #from 2=B to 3=BE
            MT[elements_in_B[np.where(randomB<P_BC)[0]]] += 2 #from 2=B to 4=C
            
            #old: from maurer model: MT[elements_in_BE[np.where(randomBE<P_BEB)[0]]] -= 1 #from 3=BE to 2=B
            #old: from maurer model: MT[elements_in_BE[np.where(randomBE>(1-P_BEC))[0]]] += 1 #from 3=BE to 4=C
            
               
            # Update end of cap position
            # -----------------------------------------------------------------
            
            old_current_cap_end = current_cap_end
            if simPa.unstable_cap_criteria > 1:
                half_testbox = int(simPa.unstable_cap_criteria/2) # int() here same as int(np.floor())
                #f = filter to check if N=unstable_cap_criteria neighboring dimer rings are unstable
                Mtest = np.convolve(MT[:]!=4,f[:],'same') #convolution to find # of places where X neighbors are all hydrolysed enough
                Mtest[0:(L_seed-half_testbox-1)] = simPa.unstable_cap_criteria                
                #wrong?:Mtest[0:(L_seed-2*half_testbox+simPa.CAP_threshold)] = unstable_cap_criteria                
                current_cap_end = np.where(np.concatenate(([0], Mtest[1:(L-half_testbox+1)])) <= simPa.CAP_threshold)[0][-1] + half_testbox
                #wrong:current_cap_end = np.where(np.concatenate(([0], Mtest[1:(L-2*half_testbox+simPa.CAP_threshold)])) <= simPa.CAP_threshold)[0][-1] + half_testbox               
                if current_cap_end < (L_seed-1): 
                    current_cap_end = L_seed - 1
            else:  # Go to faster way
                current_cap_end = L_seed + np.where(np.concatenate(([4], MT[L_seed:L])) == 4)[0][-1] # end at last element = 'C' (=4)
            
            if current_cap_end < old_current_cap_end:  # Hydrolysis cannot go back!
                current_cap_end = old_current_cap_end
                
            
            EB_profile_intermediate = EB_profile_intermediate + np.array(MT[(L-L_seed):L] == 2)

            # Update positions to be saved    
            if i % frame_taking == 0:              
                        #if barrier: #if variable barrier present
                            #if Lnew > (barrier - np.ceil(simPa.DXbarrier/simPa.dL_dimer)) and barrier_contact == 0:
                            #    barrier_contact = i*dt
                #add new data point:
                MT_length.append(L) 
                cap_end.append(current_cap_end)
                EB_comet.append(len(np.where(MT == 2)[0])) #len(elements_in_B)) #BE CAREFUL: NOW it is elements in B (not BE like in maurer model), is not 100% correct measure  
#                EB_comet.append(len(np.where(MT[(L-L_seed):L] == 2)[0])) #len(elements_in_B)) #BE CAREFUL: NOW it is elements in B (not BE like in maurer model), is not 100% correct measure  
                
                if simPa.take_EB_profiles:
                    EB_profiles_wholerun.append(np.array(EB_profile_intermediate)/frame_taking)
                    EB_profile_intermediate[0:L_seed] = 0 #reset
                
#                #check if barrier contact:
#                if barrier:
#                    if barrier_contact == 0:
#                        if np.mean(MT_length[-(1+int(1/frame_taking)):-1]) > (barrier - np.ceil(simPa.DXbarrier/simPa.dL_dimer)):
#                            barrier_contact = i*dt - np.random.rand() * frame_taking * dt 
#                
            
            # Check if catastrophe occured
            # -----------------------------------------------------------------
            
            if (L - current_cap_end) <= 0 and L > L_seed: 
                
                # Record data from run before break      
                # Without Barrier: always record ||| with barrier: only record when in contact
                
                if (simPa.barrier != 0) == (barrier_contact != 0): # =IF true/true or false/false
                    CATASTROPHE_TIMES.append(i*dt)
                    CATASTROPHE_LENGTH.append(L-L_seed)        
                    MT_length[-1] = L
                    cap_end[-1] = MT_length[-1] #set cap length to zero in case it gets lost due to binning (if frame_rate > dt)        
                    EB_comet[-1] = len(np.where(MT == 2)[0]) # NOT REAL EB, but "B" state !!              
                    
                    EB_comet_sum.append(list(reversed(EB_comet)))
                    
                    if simPa.take_EB_profiles:
                        EB_profiles.append(EB_profiles_wholerun)  # [-int(simPa.show_fraction/frame_rate_actual):-1])
                        
                    if len(MT_length) > show_fraction_frames:  # Save last piece of MT position
                        MT_length_sum[cat_number][0:show_fraction_frames] = MT_length[-(1+show_fraction_frames):-1]
                else:  # Meaning no contact with barrier
                    CATASTROPHE_TIMES.append(0)
                    CATASTROPHE_LENGTH.append(0)  
                
                    
                if simPa.barrier and barrier_contact != 0:
                    #barrier_contact_times = np.append(barrier_contact_times,[i*dt - barrier_contact])
                    #match contact time to last "frame"
                    #barrier_contact_times = np.append(barrier_contact_times,[(i // frame_taking)*frame_taking*dt - barrier_contact])
                    barrier_contact_times = np.append(barrier_contact_times,[i * dt - barrier_contact])
                    catastrophe_at_barrier.append(cat_number)
                    
                #barrier_contact_count = 0        
                barrier_contact = 0 #set contact = off:    
                
                # Output progress:
                print('\r', 'catastrophes: ', np.size(CATASTROPHE_TIMES), ' of ', simPa.no_cat, end="")
                run_time += i*dt
                
                break;
    
    return dt, MT_length_sum, CATASTROPHE_TIMES, CATASTROPHE_LENGTH, barrier_contact_times, EB_comet_sum, cap_end, frame_rate_actual, EB_profiles

        
   



        
    
               



