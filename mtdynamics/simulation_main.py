import numpy as np


# -----------------------------------------------------------------------------
# ------------ MAIN FUNCTIONS FOR MT SIMULATION ("maurer model") --------------
# -----------------------------------------------------------------------------

def mt_run(simPa):
    """
    Main simulation function.

    MTs are modeled as a 1D string of tubulin dimers.
    Those are represented by a 1D numpy array where the different dimer states
    are described by different integers.

    5 -- stable seed dimers (no hydrolysis, no disassembly)
    2 -- 'B' state dimer
    4 -- 'C' state dimer


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
    print('--- MTs growth ', dimers_per_step, ' dimers per step (D_tip = ',
          simPa.D_tip, ' -> noise STD= ', tip_noise, ' )')
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
    catastrophe_at_barrier = []
    catastrophe_washout = []
    EB_comet_sum = []
    EB_profile_intermediate = np.zeros(L_seed, dtype=np.int)
    EB_profiles = []
    show_fraction = simPa.show_fraction
    #show_begin = simPa.min_length_begin
    show_fraction_frames = int(show_fraction/frame_rate_actual)
    MT_length_sum = np.zeros((simPa.no_cat, show_fraction_frames))
    Cap_length_sum = np.zeros((simPa.no_cat, show_fraction_frames))
    MT_length_full = []
    barrier_contact_times = np.zeros(0)
    washout_times = np.zeros(0)
    MT_cap = []

    too_long = 0 # reset
    run_time = 0 # reset

    for cat_number in range(0, simPa.no_cat):
        # Reset everything:
        L = L_seed
        Lnew = L
        MT = np.ones(L, dtype=np.int)
        current_cap_end = L_seed
        MT[0:L_seed] = 5 # set state Seed
        growth = np.zeros(1)
        cap_end = [0]
        MT_length = [L]
        barrier_contact = 0 #set contact = off:

        #washout experiment:
        washout = 0 #reset

        #EB monitoring:
        EB_comet = [0]
        EB_profiles_wholerun = []

        Cap_threshold = []
        cap_crit = []

        # Interrupt run if too long (e.g. when MT is too stable)
        if too_long > simPa.too_stable_check:
            break
        if run_time > simPa.total_max_time:
            print("run stopped - was longer than time limit!")
            break

        # ---------------------------------------------------------------------
        # ------------------------ Start time loop ----------------------------
        # ---------------------------------------------------------------------
        for i in range(1, N_max + 1):  # N_max timesteps

            if i == N_max: # Check if MT grows for too long
                run_time += N_max * dt
                too_long += 1
                print("Too long growth event! (", simPa.EB, simPa.kBC, simPa.D_tip, ")")

            # -----------------------------------------------------------------
            # Option to include 'MT aging' ------------------------------------
            # -----------------------------------------------------------------
            # --> Make cap_threshold time-dependent
            if simPa.unstable_cap_time:
                criteria_new = (simPa. unstable_cap_start - simPa.unstable_cap_end) * np.exp(-1*simPa.unstable_cap_rate*(i*dt)) + simPa.unstable_cap_end
                simPa.unstable_cap_criteria = int(np.rint(criteria_new))
                cap_crit.append(simPa.unstable_cap_criteria)
                f = np.ones(simPa.unstable_cap_criteria) #filter for end-of-cap testing
            else:
                cap_crit.append(simPa.unstable_cap_criteria)
                f = np.ones(simPa.unstable_cap_criteria) #filter for end-of-cap testing

            # --> Make D_tip time-dependent
            if simPa.D_tip_time:
                D_tip_new = (simPa.D_tip_end - simPa.D_tip_start) * (1 - np.exp(-1*simPa.D_tip_rate_T*(i*dt))) + simPa.D_tip_start
                tip_noise = (2*dt*D_tip_new)**0.5/(simPa.dL_dimer*1000)

            # --> Make D_tip length-dependent
            if simPa.D_tip_length:
                D_tip_new = (simPa.D_tip_end - simPa.D_tip_start) * (1 - np.exp(-1*simPa.D_tip_rate_L*L)) + simPa.D_tip_start
                tip_noise = (2*dt*D_tip_new)**0.5/(simPa.dL_dimer*1000)

            # --> Make k_BC time-dependent
            if simPa.kBC_time:
                kBC_new = (simPa.kBC_end - simPa.kBC_start) * (1 - np.exp(-1*simPa.kBC_rate*(i*dt))) + simPa.kBC_start
                P_BC = kBC_new * dt

            # -----------------------------------------------------------------
            # Tip growth ------------------------------------------------------
            # -----------------------------------------------------------------
            deterministic_growth = deterministic_growth + P_growth #here: only to time "growth event"

            if deterministic_growth >= 1: # if growth event:
                deterministic_growth = deterministic_growth-1
                growth = np.round(dimers_per_step + tip_noise*np.random.randn(1))
                """ Option to include dimer-type dependent growth dynamics
                if MT[L-1] == 1:
                    growth = round(1+ simPa.noise_STD_A*np.random.randn(1)) #%= random number with std= 1; Round to draw discrete number
                elif MT[L-1] == 4:
                    growth = round(P_growth_C + simPa.noise_STD_C*np.random.randn(1))
                elif MT[L-1] == 2:
                    growth = round(P_growth_BE + simPa.noise_STD_BE*np.random.randn(1))
                else: # for MT[L-1] == 3
                    growth = round(P_growth_B + simPa.noise_STD_B*np.random.randn(1))"""

                # Option to simulate 'washout' experiment
                if washout == 1:
                    if growth > 0:
                        growth = 0 #do not allow growth after tubulin washout

                # Actual growth --> addition of int(growth) dimers to MT
                if abs(growth) > 0:
                    Lnew = L + int(growth)

                    if Lnew < L:  # MT shrinks
                        if Lnew < L_seed:
                            Lnew = L_seed # +1 to forbid catastrophes at seed
                        else:
                            MT[Lnew:L+1] = 0
                    elif Lnew > L:  # MT grows

                        #  Change size of MT array when necessary
                        if Lnew > len(MT):
                            MT = np.append(MT, np.zeros(Lnew - L + 500, dtype=np.int))

                        # Option to simlate growth against rigid barrier
                        if simPa.barrier:
                            if Lnew > (simPa.barrier - np.ceil(simPa.DXbarrier/simPa.dL_dimer)) and barrier_contact == 0: #contact if closer than 0.1um, 0.2µm, 0.3µm?
                                barrier_contact = i*dt
                                #barrier_contact = (i // frame_taking) * frame_taking * dt #time of first 'contact' with wall

                            if Lnew > simPa.barrier:

                                # Default: NO brownian ratchet-like stalling:
                                # Rigid barrier --> Still possible to grow, however only till barrier.
                                Lnew = simPa.barrier

                                # Alternative: Brownian ratchet-like stalling:
                                # Spring like response at barrier
                                #F_barrier = simPa.k_barrier*(Lnew - barrier)
                                # Mogilner: V=Vmax*exp(force*dl/(k_b*T)), k_bT=4.1pN/nm divide by 8nm/dimer
                                #Lnew = L + (Lnew-L)*(np.random.rand(1)< np.exp(- F_barrier*1000*simPa.dL_dimer/4.1))[0]
                                if Lnew > L:
                                    MT[L:Lnew] = 2 #start with "B" state
                            else:
                                MT[L:Lnew] = 2 #start with "B" state
                        else:
                            MT[L:Lnew] = 2 #start with "B" state
                    L = Lnew
                    L = max(L, L_seed) #avoid L < L_seed

                    if simPa.washout:
                        if (i * dt) >= simPa.washout_time and washout == 0:
                            washout = 1
                            washout_time = i * dt

            # -----------------------------------------------------------------
            # Hydrolysis ------------------------------------------------------
            # Based on ~ Maurer 2014 A->B (->BE) -> C
            # Here: Simplified multi-step reaction, ignoring fast A -> B:
            # B -> C
            # -----------------------------------------------------------------
            elements_in_B = np.where(MT == 2)[0]
            randomB = np.random.rand(len(elements_in_B))

            # Transition from dimer state 2=B to 4=C
            MT[elements_in_B[np.where(randomB < P_BC)[0]]] += 2

            """
            # -----------------------------------------------------------------
            # Possible alternative:
            # Full multi-step reaction ~ Maurer 2014 A->B (->BE) -> C
            # -----------------------------------------------------------------

            elements_in_A = np.where(MT == 1)[0]
            elements_in_B = np.where(MT == 2)[0]
            # based on maurer model:
            elements_in_BE = np.where(MT == 3)[0]
            elements_in_C = np.where(MT == 4)[0]

            # Shift between different stages:
            # 1- throw dices:
            randomA = np.random.rand(len(elements_in_A))
            randomB = np.random.rand(len(elements_in_B))
            # based on maurer model:
            randomBE = np.random.rand(len(elements_in_BE))
            randomC = np.random.rand(len(elements_in_C))

            # 2- change state of elements based on dices:
            MT[elements_in_A[np.where(randomA<P_AB)[0]]] += 1 #go from 1=A to 2=B
            MT[elements_in_A[np.where(randomA<P_ABE)[0]]] += 1 # and from 2=B to 3=BE
            MT[elements_in_A[np.where(randomA>(1-P_AC))[0]]] += 3 #from 1=A to 4=C
            MT[elements_in_B[np.where(randomB<P_BBE)[0]]] += 1 #from 2=B to 3=BE
            MT[elements_in_B[np.where(randomB<P_BC)[0]]] += 2 #from 2=B to 4=C
            MT[elements_in_BE[np.where(randomBE<P_BEB)[0]]] -= 1 #from 3=BE to 2=B
            MT[elements_in_BE[np.where(randomBE>(1-P_BEC))[0]]] += 1 #from 3=BE to 4=C
            """


            # -----------------------------------------------------------------
            # Update end of cap position
            # -----------------------------------------------------------------
            old_current_cap_end = current_cap_end
            if simPa.unstable_cap_criteria > 1:
                half_testbox = int(simPa.unstable_cap_criteria / 2)

                # Convolution to find num of places where X neighbors are all hydrolysed enough
                Mtest = np.convolve(MT[:] != 4, f[:], 'same')
                Mtest[0:(L_seed-half_testbox - 1)] = simPa.unstable_cap_criteria

                current_cap_end = np.where(np.concatenate(([0], Mtest[1:(L - half_testbox + 1)])) \
                                           <= simPa.CAP_threshold)[0][-1] + half_testbox
                if current_cap_end < (L_seed - 1):
                    current_cap_end = L_seed - 1
            else:  # Go to faster way
                current_cap_end = L_seed + np.where(np.concatenate(([4], MT[L_seed:L])) == 4)[0][-1] # end at last element = 'C' (=4)

            if current_cap_end < old_current_cap_end:  # Hydrolysis cannot go back!
                current_cap_end = old_current_cap_end

            EB_profile_intermediate = EB_profile_intermediate + np.array(MT[(L-L_seed):L] == 2)

            # Update positions to be saved
            if i % frame_taking == 0:
                # Add new data point:
                MT_length.append(L)
                cap_end.append(current_cap_end)
                EB_comet.append(len(np.where(MT == 2)[0])) #BE CAREFUL: here it is elements in B (not BE like in maurer model), is not 100% correct measure for EB

                if simPa.take_EB_profiles:
                    EB_profiles_wholerun.append(np.array(EB_profile_intermediate)/frame_taking)
                    EB_profile_intermediate[0:L_seed] = 0 #reset

            # -----------------------------------------------------------------
            # Check if catastrophe occured
            # -----------------------------------------------------------------
            if (L - current_cap_end) <= 0 and L > L_seed:

                # Record data from run before break
                # Without Barrier: always record ||| with barrier: only record when in contact

                if (simPa.barrier != 0) == (barrier_contact != 0): # = IF true/true or false/false
                    CATASTROPHE_TIMES.append(i*dt)
                    CATASTROPHE_LENGTH.append(L-L_seed)
                    MT_length[-1] = L
                    cap_end[-1] = MT_length[-1] #set cap length to zero in case it gets lost due to binning (if frame_rate > dt)
                    EB_comet[-1] = len(np.where(MT == 2)[0]) # NOT REAL EB, but "B" state !!

                    #  Collect cap position and EB comet size
                    MT_length_full.append(MT_length)
                    MT_cap.append(cap_end)
                    EB_comet_sum.append(list(reversed(EB_comet)))

                    # Collect current cap criterium
                    Cap_threshold.append(cap_crit)

                    if simPa.take_EB_profiles:
                        EB_profiles.append(EB_profiles_wholerun)

                    # Save last piece of MT position
                    if len(MT_length) > show_fraction_frames:
                        MT_length_sum[cat_number][0:show_fraction_frames] = MT_length[-(1+show_fraction_frames):-1]
                        Cap_length_sum[cat_number][0:show_fraction_frames] = cap_end[-(1+show_fraction_frames):-1]

                else:  # Meaning no contact with barrier
                    CATASTROPHE_TIMES.append(0)
                    CATASTROPHE_LENGTH.append(0)

                if simPa.barrier and barrier_contact != 0:
                    barrier_contact_times = np.append(barrier_contact_times, [i * dt - barrier_contact])
                    catastrophe_at_barrier.append(cat_number)

                if simPa.washout and washout != 0:
                    washout_times = np.append(washout_times, [i*dt-washout_time])
                    catastrophe_washout.append(cat_number)

                barrier_contact = 0 # reset contact = off:

                # Output progress:
                print('\r', 'catastrophes: ', np.size(CATASTROPHE_TIMES), ' of ', simPa.no_cat, end="")
                run_time += i*dt

                break

    return dt, MT_length_sum, MT_length_full, CATASTROPHE_TIMES, CATASTROPHE_LENGTH, \
        barrier_contact_times, EB_comet_sum, MT_cap, Cap_threshold, frame_rate_actual, \
        EB_profiles, washout_times, catastrophe_washout, Cap_length_sum
