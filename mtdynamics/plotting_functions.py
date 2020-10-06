import numpy as np
import random
import time
from matplotlib import pyplot as plt
import pandas as pd
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from scipy.special import gammainc
from itertools import islice

# Initialize random number generator:
np.random.seed(int(100*time.perf_counter()))

## ----------------------------------------------------------------------------
## --------------------------- Plotting functions -----------------------------
## ----------------------------------------------------------------------------

def KMG_analysis(events, times):
    """ Kaplan-Meier survivial funciton with greenwood's forumla to estimate variance.

    Args:
    -------
    events:
    times:
    """
    S = np.ones(len(times)+1) #survival probability
    S[0] = 1
    V = np.ones(len(times)+1) #variance of S (Greenwood's formula)
    V_cumulative = np.zeros(len(times)+1)
    V[0] = 0

    num_of_events = np.sum(events)
    for i, times in enumerate(times):
        S[i] = S[i-1] * (1 - events[i-1]/num_of_events)
        V_cumulative[i] = V_cumulative[i-1] + events[i-1]/(num_of_events*(num_of_events-events[i-1]))
        V[i] = S[i]**2 * V_cumulative[i]

    return S, V



def distribution_alternatives(distribution, num_alter, overlap):
    """ include "new" error function for cumulative distributions

    Args:
    ------
    distribution: list, array
        Original data to compare to.
    num_alter: int
        Number of alternative distributions to generate
    overlap: int, float
        What fraction to (randomly) draw from original one.
    """
    if overlap < 1:
        # Assume that fraction between 0 and 1 was given
        overlap = int(np.ceil(overlap*len(distribution)))
    else:
        # Assume that overlap was given as number of desired elements (not percentage!)
        print("Overlap given will not be interpreted as percentage!")
        overlap = int(overlap)

    num_initial = len(distribution)

    distribution_alternatives = np.zeros((overlap, num_alter))

    if distribution[0] == 0: #if (cumulative) distribution sorted and startins from 0
        for i in range(0, num_alter):
            random_index = random.sample(range(1, num_initial), overlap-1)
            random_index = np.append(random_index, 0)
            random_index = np.sort(random_index)
            distribution_alternatives[:,i] = distribution[random_index]

    else:
        for i in range(0, num_alter):
            random_index = random.sample(range(0, num_initial), overlap)
            random_index = np.sort(random_index)
            distribution_alternatives[:,i] = distribution[random_index]

    return distribution_alternatives


def distribution_compare(Cum_hist1,
                         Cum_hist2,
                         num_interpol=10000):
    """
    Function to compare to input distributions.

    Args:
    --------
    Cum_hist1: list, array #TODO: check!
    Cum_hist2: list, array #TODO: check!
    num_interpol: int, optional
        Number of interpolation bins. Default = 10000.
    """
    y1 = 1/(len(Cum_hist1)-1) * np.arange(0, len(Cum_hist1) , 1)
    y2 = 1/(len(Cum_hist2)-1) * np.arange(0, len(Cum_hist2) , 1)
    fit1 = interp1d(Cum_hist1, y1, kind='nearest')
    fit2 = interp1d(Cum_hist2, y2, kind='nearest')
    xnew = np.linspace(0,min(max(Cum_hist1),max(Cum_hist2)), num=num_interpol) # only look at first 95% (ignore weird end)
    return (np.mean((fit1(xnew) - fit2(xnew))**2))


def valid_EB_runs(simPa,
                  EB_comet_sum,
                  barrier_contact_times = []):
    """ Function to select valid runs (runs longer than min_length_run).

    Args:
    -------
    simPa: parameter set
        Simulation parameters in "ParameterSet" format.
    EB_comet_sum: list
        List containing EB counts in comet. #TODO: check
    barrier_contact_times: list
        List containing barrier contact times.
    """

    # Select valid runs
    b = []
    if simPa.barrier: # if barrier present
        for a in range(0, len(EB_comet_sum)):
            b.append(len(EB_comet_sum[a]) * simPa.frame_rate_actual - barrier_contact_times[a])
    else:
        for a in range(0, len(EB_comet_sum)):
            b.append(len(EB_comet_sum[a]) * simPa.frame_rate_actual)

    valid_runs = np.where(np.array(b) > simPa.min_length_run)[0]

    return valid_runs


def analyse_EB_signal(simPa,
                      EB_comet_sum,
                      barrier_contact_times):
    """ Function to analyse EB signal

    Args:
    -------
    simPa: parameter set
        Simulation parameters in "ParameterSet" format.
    EB_comet_sum: list
        List containing EB counts in comet. #TODO: check
    barrier_contact_times: list
        List containing barrier contact times.
    """

    # Select valid runs
    valid_runs = valid_EB_runs(simPa, EB_comet_sum, barrier_contact_times)

    max_barrier_contact_frames = int(round(np.max(barrier_contact_times/simPa.frame_rate_actual),0))
    min_length_run_frames = int(simPa.min_length_run/simPa.frame_rate_actual)
    frame_window = min_length_run_frames + max_barrier_contact_frames

    EB_signal = np.zeros((len(valid_runs), frame_window+1)) #simPa.min_length_run+1+max_barrier_contact)) #put individual runs into one np.array
    normalize_EB_signal = np.zeros(frame_window+1) #simPa.min_length_run+1+max_barrier_contact)

    for a in range(0,len(valid_runs)):
        frame_barrier_contact = int(np.round(barrier_contact_times[valid_runs[a]]/simPa.frame_rate_actual,0))
        EB_signal[a][(max_barrier_contact_frames-frame_barrier_contact):frame_window] \
        = np.array(EB_comet_sum[valid_runs[a]])[0:(min_length_run_frames+frame_barrier_contact)]
        normalize_EB_signal[(max_barrier_contact_frames-frame_barrier_contact):frame_window] +=1

    EB_signal_average = np.sum(EB_signal, axis=0)
    EB_signal_average = EB_signal_average/normalize_EB_signal

    return EB_signal, EB_signal_average, max_barrier_contact_frames, min_length_run_frames, frame_window


def analyse_EB_profile(simPa,
                       MT_length_full,
                       EB_profiles,
                       w_size):
    """ Calculate the mean GTP/GDP-Pi (or EB) profile at the microtubule end during steady-state growth.
    The resulting profile is not convolved wit a Gaussian and thus represents the theoretical profile.

    Args:
    ------
    simPa: parameter set
        Simulation parameters in "ParameterSet" format.
    TODO: add documentation
    ...
    """

    # Initialize arays
    v_mean = []
    EB_mean = []

    # Microscope parameters
    resolution = 0.25 #um

    # Number of dimers in PSF
    resolution_dimers = int(np.ceil(resolution/simPa.dL_dimer))

    # Loop over number of simulated microtubules
    count = 0
    for num_run in range(len(MT_length_full)):

        # Obtain the time and length arrays
        time = np.arange(0, len(MT_length_full[num_run]), 1) * simPa.frame_rate_actual
        MT_length = (np.asarray(MT_length_full[num_run]) - MT_length_full[num_run][0]) * simPa.dL_dimer *1000

        # Calculate mean growth speed
        v = np.polyfit(time, MT_length, 1)
        v_mean.append(v[0])

        if simPa.steady_state_analysis:
            # Find the local mean growth speeds in order to exclude pausing state from the profile analysis
            if len(MT_length) > w_size:
                di = 0
                v_fit = np.zeros(len(MT_length) - w_size + 1)
                for i in window(MT_length, w_size):
                    v_fit[di] = np.polyfit(np.linspace(0, simPa.frame_rate_actual*w_size-1, w_size), i, 1)[0]
                    di = di + 1
            else: v_fit = []

            # Set velocity threshold
            v_thres = 0.6

            # Identify steady-state growth events in the trace
            matches = [i for i, x in enumerate(v_fit) if x > v_thres*v_mean[-1]]
            matches = np.asarray(matches) + w_size//2

            if matches.size > 0:
                for mm in matches:
                     # Extend the EB profile array
                     EB_new = np.append(EB_profiles[num_run][mm], np.zeros(resolution_dimers))
                     if count > 0:
                         EB_mean = EB_mean*(count/(count+1)) + EB_new*(1/(count+1))
                         count += 1
                     else:
                        EB_mean = EB_new
                        count += 1

        else: # Include the complete growth trajectory to calculate the mean profile

            for mm in range(len(EB_profiles[num_run])-1):
                # Extend the EB profile array
                EB_new = np.append(EB_profiles[num_run][mm], np.zeros(resolution_dimers))
                if count > 0:
                    EB_mean = EB_mean*(count/(count+1)) + EB_new*(1/(count+1))
                    count += 1
                else:
                    EB_mean = EB_new
                    count += 1

    return EB_mean, v_mean



# -----------------------------------------------------------------------------
# -------------------------- Small helper functions ---------------------------
# -----------------------------------------------------------------------------

def frange(start, stop, step):
    """ Function as alternative for "range, since "range" does not support floats.
    """
    i = start
    while i < stop:
        yield i
        i += step


def list_dim(lst):
    """ Function to return the dimension of a list (e.g. nested list).
    """
    if not type(lst) == list:
        return 0
    return len(lst) + list_dim(lst[0])


# Define Gaussian distribution
def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))

def gamma_cdf(x, n, r):
    return gammainc(n, r*x)

def exp_cdf(x, k):
    return 1 - np.exp(-k*x)

# Define a sliding window
def window(seq, n):
    """ Generater that returns a sliding window (of width n) over data from the iterable/
    s -> (s0,s1,...s[n-1]), (s1,s2,...,sn), ...
    """
    it = iter(seq)
    result = tuple(islice(it, n))
    if len(result) == n:
        yield result
    for elem in it:
        result = result[1:] + (elem,)
        yield result

# -----------------------------------------------------------------------------
# ----------------------------- Figure functions ------------------------------
# -----------------------------------------------------------------------------

# Figure styles
from matplotlib.font_manager import FontProperties
font = FontProperties()
font.set_family('sans-serif')
font.set_style('normal')
font.set_weight('light')


def fig_sim_verification(simPa,
                         file_figure,
                         num_fig,
                         MT_length_full,
                         cap_end,
                         d_steps=1):
    """ Compare the fixed parameter v_g and D_tip with the simulated results.

    d_steps = 2 # in time_steps

    Args:
    ------
    simPa: parameter set
        Simulation parameters in "ParameterSet" format
    file_figure: str
        Folder for storing figures and data
    num_fig: int
        Figure number
    MT_length_full: numpy array ? #TODO:check
        ...
    cap_end
        ...
    d_steps: ...
    """

    # Calculate the growth fluctuations
    v_fluc = []
    c_fluc = []

    # Remove cap=seed position
    L_seed = int(np.ceil(simPa.tip_window/simPa.dL_dimer))
    for i in range(len(cap_end)):
        cap_temp = cap_end
        index = np.argmax(np.asarray(cap_end[i]) > L_seed)
        del cap_temp[i][:index]

    cap_end_ss = cap_temp

    for i in range(len(MT_length_full)):
        v_fluc.extend(np.diff(MT_length_full[i]))
        c_fluc.extend(np.diff(cap_end_ss[i]))

    sample_size = len(c_fluc)

    c_fluc = np.sort(c_fluc)
    index = np.argmax(np.asarray(c_fluc) > 0)
    c_fluc = c_fluc[int(index):]

    v_fluc = np.asarray(v_fluc)*(simPa.dL_dimer*1000)
    c_fluc = np.asarray(c_fluc)*(simPa.dL_dimer*1000)

    # Calculate the growth distribution based on the fixed parameters
    mu = simPa.growth_rate_one*simPa.frame_rate_actual # mean growth rate in dimers/frame
    sig = ((2*simPa.D_tip*(simPa.frame_rate_actual*d_steps))**0.5)
    x = np.arange(mu-5*sig, mu+5*sig, 1)
    G = gaussian(x, mu, sig)
    G = G / np.sum(G) # normalize gaussian

    # Plot the results
    fig , (ax1, ax2) = plt.subplots(1,2, figsize=(12, 7))

    ax1.hist(v_fluc, bins = 60, density = True, color = "skyblue", label = "simulated data")
    ax1.plot(x, G, 'r', label = "theoretical distribution")

    ax2.hist(c_fluc, bins = 60, density = True, color = "skyblue", label = "simulated data")

    move =  len(c_fluc)/sample_size
    pause = 1 - move
    step_mean = np.mean(c_fluc)
    step_std =  np.std(c_fluc)

    if simPa.record_data:
            filename = file_figure + '_fig' + str(int(num_fig))
            plt.savefig(filename + '.eps', format='eps', dpi=1000)
            plt.savefig(filename + '.png', format='png', dpi=200)

    plt.show()

    print('Pausing probability: %.2f' %pause)
    print('Step probability: %.2f' %float(1-pause))
    print('Mean step size: %.1f +- %.1f nm' %(step_mean, step_std))
    print('Mean pausing duration: %.2f sec ' %float(step_mean/(simPa.growth_speed*1000/60)))


def fig_cat_dist(simPa,
                 file_figure,
                 num_fig,
                 catastrophe_times,
                 Cum_dist_compare):
    """ Catastrophe distribution compared to data

    Args:
    ------
    simPa: parameter set
        Simulation parameters in "ParameterSet" format
    file_figure: str
        Folder for storing figures and data
    num_fig: int
        Figure number
    catastrophe_times: numpy array
        Array of catastrophe times
    Cum_hist_compare: list, array #TODO: check
        Cumulative catastrophe time distribution for comparison
    """

    if not isinstance(catastrophe_times, np.ndarray):
        print('Catastrophe times input format must be numpy array!')
        catastrophe_times = np.zeros(0)

    if catastrophe_times.shape[0] > 1:
        tau_c = np.mean(catastrophe_times) #i*dt/len(catastrophe_times)
        print('Mean catastrophe time: %.2f s' %tau_c)

        n_bins = int(np.ceil(len(catastrophe_times)/10))

        fig = plt.figure(3)
        plt.clf()

        ## Compare to data
        bins=np.histogram(np.hstack((Cum_dist_compare,catastrophe_times)), bins=n_bins)[1] #get the bin edges
        Hist_exp, edges_exp = np.histogram(Cum_dist_compare, bins = bins)
        bin_width = edges_exp[1]
        plt.bar((edges_exp[:-1] + bin_width/2) , np.float_(Hist_exp)/(sum(Hist_exp)), bin_width, alpha=0.5, color='gray')
        Hist, edges = np.histogram(catastrophe_times, bins = bins)
        plt.plot((edges[1:] -edges[1]/2), np.float_(Hist)/(sum(Hist)),'r-', linewidth=1.0)

        #plt.title('Catastrophe distribution')
        plt.xlabel('time [s]')
        plt.ylabel('fraction of event')

        fig.suptitle('Catastrophe distribution', fontsize=14, fontweight='bold')
        plt.ax = fig.add_subplot(111)
        fig.subplots_adjust(top=0.9)

        #Add parameters to figure
        figtext = ['$v_{g} = %.2f \mu m/min$' %float(simPa.growth_rate_one*(60*simPa.dL_dimer))]
        figtext = ['$EB = %.2f \mu M$' %float(simPa.EB)]
        figtext.append('$D_{tip} = %.2f nm^2/s$)' %simPa.D_tip)
        figtext.append('Cap unstable when in state "C" ')
        figtext.append('in %r out of %r dimer layers.' %(int(simPa.unstable_cap_criteria-simPa.CAP_threshold),int(simPa.unstable_cap_criteria)))
        figtext.append('Tip states:B->C with the rates:' )
        figtext.append('$k_{hyd} = %.3f s^{-1}$' %(simPa.kBC ))
        figtext.append('Results (n = %d) -------------------------------' %len(catastrophe_times))
        figtext.append(r'$\tau_{C} = %.2f s$' %tau_c)

        figDX = 0.045
        for m in range(len(figtext)):
            plt.ax.text(0.4, 0.9-m*figDX, figtext[m], fontproperties=font,
                verticalalignment='bottom', horizontalalignment='left',
                transform=plt.ax.transAxes, color='black', fontsize=8)

        if simPa.record_data:
            filename = file_figure + '_fig' + str(int(num_fig))
            plt.savefig(filename+'.eps', format='eps', dpi=1000)
            plt.savefig(filename+'.png', format='png', dpi=200)

        plt.show()

    else:
        print('No proper input found.')


def fig_cat_cumulative(simPa, file_figure, num_fig, Cum_dist, Cum_dist_compare = [0]):
    """ Plot cumulative catastrophe distribution (or barrier contact time distribution).
    Compare to (experimental) data if given.

    Args:
    -------
    simPa: parameter set
        Simulation parameters in "ParameterSet" format.
    file_figure: str
        Folder for storing figures and data.
    num_fig: int
        Figure number.
    Cum_dist: list, array #TODO:check
        Cumulative catastrophe time distribution.
    Cum_hist_compare: list, array #TODO:check
        Cumulative catastrophe time distribution for comparison.
    """

    fig = plt.figure(1, figsize=(12, 7))
    plt.clf()

    # Check input cumulative distribution
    if isinstance(Cum_dist, list) and list_dim(Cum_dist) > 1 and list_dim(Cum_dist[0]) > 0:
        if isinstance(Cum_dist[0], np.ndarray):
            print(list_dim(Cum_dist), ' different cumulative distributions found. ')
        else:
            print('Error: Input cumulative distributions must be numpy arrays or lists of numpy arrays.' )
    elif isinstance(Cum_dist, list) and list_dim(Cum_dist) == 1 and isinstance(Cum_dist[0], np.ndarray):
        pass;
    elif isinstance(Cum_dist, np.ndarray):
        Cum_dist = [Cum_dist] #put numpy array into list
    else:
        print('Error: Input cumulative distributions must be numpy arrays or lists of numpy arrays.' )

    if len(Cum_dist_compare) > 1: # i.e.if comparison data is given
        if isinstance(Cum_dist_compare, list):
            if list_dim(Cum_dist_compare) == 1:
                comparing_index = np.zeros(list_dim(Cum_dist))
            elif list_dim(Cum_dist_compare) == list_dim(Cum_dist):
                #Assume that one comparison distribution given for each Cum_dist + same ordering
                print('Function assumes same pairing of distributions: 1-1, 2-2, ... ')
                comparing_index = np.arange(0, list_dim(Cum_dist_compare))
            else:
                print('Error: Dimension of comparison distribution(s) does not match.' )
                comparing_index = []
        elif isinstance(Cum_dist_compare, np.ndarray):
            Cum_dist_compare = [Cum_dist_compare]
            comparing_index = np.zeros(list_dim(Cum_dist))
        else:
            print('Error: Input distributions must be numpy arrays or lists of numpy arrays.' )
            comparing_index = []

    if list_dim(Cum_dist) > 1:
        c_range = 1/(list_dim(Cum_dist)-1)
    else:
        c_range = 1
    print(c_range)
    for i, Cum_dist in enumerate(Cum_dist):
        print((0.95*(i+1)*c_range, 0.1, 0.1))
        plt.step(Cum_dist, 1/(len(Cum_dist)-1) * np.arange(0, len(Cum_dist) , 1),
                 where='post', color=(0.95-0.7*(i)*c_range, 0.1, 0.1 + 0.8*(i)*c_range), linewidth=1.5, label='model results')

        if len(Cum_dist_compare) > 1:
            Cum_dist_compare_selected = Cum_dist_compare[int(comparing_index[i])]

            #generate and draw distributions of same length as experimental data
            print(Cum_dist_compare_selected.shape)
            print(comparing_index)
            overlap = Cum_dist_compare_selected.shape[0]
            print('overlap: ', overlap)
            num_distributions = 100
            if overlap < len(Cum_dist): #needed: more simulation data points than experimental ones
                Cum_dist_variants = distribution_alternatives(Cum_dist, num_distributions, overlap)

                for m in range(0, num_distributions):
                    plt.step(Cum_dist_variants[:,m], 1/(overlap-1) * np.arange(0, overlap , 1),
                             where='post', color=(0.95-0.7*(i)*c_range, 0.3, 0.1 +0.8*(i)*c_range), alpha=0.25, linewidth=1.0)
            plt.step(Cum_dist_compare_selected, 1/(len(Cum_dist_compare_selected)-1) * np.arange(0, len(Cum_dist_compare_selected) , 1),
                         where='post', color='black', linewidth=1.5, label='experimental data')

    if simPa.barrier:
        plt.title('Cumulative contact-time distribution')
    else:
        plt.title('Cumulative catastrophe time distribution')

    plt.xlabel('time [s]')
    plt.legend(fontsize=14)

    plt.ax = fig.add_subplot(111)
    fig.subplots_adjust(top=0.9)

    figtext = ['$v_{g} = %.2f \mu m/min$' %float(simPa.growth_rate_one*(60*simPa.dL_dimer))]
    figtext.append('$EB = %.2f \mu M$' %float(simPa.EB))
    figtext.append('$D_{tip} = %.2f nm^2/s$)' %simPa.D_tip)
    figtext.append('Cap unstable when in state "C" ')
    figtext.append('in %r out of %r dimer layers.' %(int(simPa.unstable_cap_criteria-simPa.CAP_threshold),int(simPa.unstable_cap_criteria)))
    figtext.append('Tip states: B->C with the rates:' )
    figtext.append('$k_{hyd} = %.3f s^{-1}$' %simPa.kBC)
    figtext.append('dt = %.2f s  ||  V = %.2f um/s' %(simPa.dt, simPa.growth_rate_one*60*simPa.dL_dimer))
    figtext.append('actual frame rate = %.2f /s' %simPa.frame_rate_actual)

    figDX = 0.045
    for m in range(len(figtext)):
        plt.ax.text(0.6, 0.65-m*figDX, figtext[m], fontproperties=font,
            verticalalignment='bottom', horizontalalignment='left',
            transform=plt.ax.transAxes, color='black', fontsize=10)

    if simPa.record_data:
        file_figure = file_figure + '_fig' + str(int(num_fig))
        plt.savefig(file_figure +'.pdf', format='pdf', dpi=1000) #, transparent=True)
        plt.savefig(file_figure +'.png', format='png', dpi=200) #, transparent=True )
        file_csv = file_figure[:-10] + "EB" + str(simPa.EB*1000)[:-2] + "_" + str(simPa.kBC) + "_" + str(simPa.D_tip) + ".csv"
        Cum_dist_pd = pd.DataFrame(np.round(Cum_dist,2))
        Cum_dist_pd.to_csv(file_csv, header=None, index=None)

    plt.show()


def fig_EB_at_barrier(simPa,
                      file_figure,
                      num_fig,
                      EB_comet_sum,
                      barrier_contact_times):
    """ Plot EB intensity (here = elements in state "B") before and at barrier contact.

    Args:
    -------
    simPa: parameter set
        Simulation parameters in "ParameterSet" format.
    file_figure: str
        Folder for storing figures and data.
    num_fig: int
        Figure number.
    EB_comet_sum: list, array #TODO:check
        Number of "B"s during a time window before catastophe.
    barrier_contact_times: list, array #TODO:check
        List/array containing barrier contact times.
    """

    # Select valid runs
    valid_runs = valid_EB_runs(simPa, EB_comet_sum, barrier_contact_times)

    EB_signal, EB_signal_average, max_barrier_contact_frames, min_length_run_frames, frame_window = analyse_EB_signal(simPa,
                                                                                 EB_comet_sum, barrier_contact_times)

    def func(x, a, b, c):
        return a * np.exp(-b * x) + c

    xdata = simPa.frame_rate_actual * np.arange(0,int(max_barrier_contact_frames/2.5))
    ydata = EB_signal_average[(max_barrier_contact_frames-int(max_barrier_contact_frames/2.5)):max_barrier_contact_frames][::-1]
    # Make fit:
    popt, pcov = curve_fit(func, xdata, ydata, p0=(np.max(ydata), simPa.kBC, 1),maxfev=1000)

    norm = np.mean(EB_signal_average[-50:-2])

    # Start plotting
    if num_fig != 0:
        plt.figure(num_fig)
        plt.clf()
        for a in range(0,len(valid_runs)):
            plt.plot(simPa.frame_rate_actual * np.arange(-min_length_run_frames+1,max_barrier_contact_frames+1),
                     EB_signal[a][0:frame_window][::-1]/norm,color=plt.cm.Reds(0.3+0.7*a/len(valid_runs)))
        plt.plot(simPa.frame_rate_actual * np.arange(-min_length_run_frames+1,max_barrier_contact_frames+1),
                 EB_signal_average[0:frame_window][::-1]/norm,'black', linewidth=3.0)

        plt.title("number of B's")
        plt.xlabel('time before barrier contact [s]');
        plt.ylabel('EB comet intensity');
        plt.xlim(-10, 20)
        print(popt)
        plt.plot(xdata, func(xdata, *popt)/norm, 'c--',)
        plt.text(8,1.5*max(ydata/norm),'decay rate (exp. fit): %.2f' %popt[1])

        if simPa.record_data:
                filename = file_figure + '_fig' + str(int(num_fig))
                plt.savefig(filename+'.eps', format='eps', dpi=1000)
                plt.savefig(filename+'.png', format='png', dpi=200)

        plt.show()

    return popt[1]


def fig_EB_before_cat(simPa, file_figure, num_fig, EB_comet_sum, barrier_contact_times=[]):
    """ Plot EB intensity (here = elements in state "B") before catastrophe.

    Args:
    -------
    simPa: parameter set
        Simulation parameters in "ParameterSet" format.
    file_figure: str
        Folder for storing figures and data.
    num_fig: int
        Figure number.
    EB_comet_sum: list, array #TODO:check
        Number of "B"s during a time window before catastophe.
    barrier_contact_times: list, array #TODO:check
        List/array containing barrier contact times.
        Needed if barrier not False
    """

    EB_output = []
    min_length_run_frames = int(simPa.min_length_run/simPa.frame_rate_actual)

    # Select valid runs
    valid_runs = valid_EB_runs(simPa, EB_comet_sum, barrier_contact_times)

    EB_signal_before_cat = np.zeros((len(valid_runs), min_length_run_frames+1)) # put individual runs into one np.array
    for a in range(0,len(valid_runs)):
        EB_signal_before_cat[a][0:min_length_run_frames] = \
        np.array(EB_comet_sum[valid_runs[a]])[0:min_length_run_frames]

    EB_signal_mean = np.mean(EB_signal_before_cat, axis=0)
    EB_signal_std = np.std(EB_signal_before_cat, axis=0)

    # 95% confidence intervals
    CI_upper = EB_signal_mean + 1.96*(EB_signal_std/np.sqrt(len(valid_runs)))
    CI_lower = EB_signal_mean - 1.96*(EB_signal_std/np.sqrt(len(valid_runs)))

    # Normalize the EB signal
    norm = np.mean(EB_signal_mean[-50:])

    # Combine time and mean EB value in single output array
    EB_output = np.vstack((CI_upper,CI_lower))
    EB_output = np.vstack((EB_signal_mean, EB_output))
    EB_output = np.flip(np.delete(EB_output,-1,1), 1)/norm
    EB_output = np.vstack((simPa.frame_rate_actual*np.arange(-min_length_run_frames,0),EB_output)).T

    # Start plotting
    if num_fig != 0:
        plt.figure(num_fig)
        plt.clf()
        plt.plot(simPa.frame_rate_actual * np.arange(-min_length_run_frames,0),
                 EB_signal_mean[0:min_length_run_frames][::-1]/norm,'red', linewidth=2.0)

        for a in range(0,len(valid_runs)):
            plt.plot(simPa.frame_rate_actual * np.arange(-min_length_run_frames,0),
                     EB_signal_before_cat[a][0:min_length_run_frames][::-1]/norm,color=plt.cm.Reds(0.3+0.7*a/len(valid_runs)))
            plt.plot(simPa.frame_rate_actual * np.arange(-min_length_run_frames,0),
                     EB_signal_mean[0:min_length_run_frames][::-1]/norm,'black', linewidth=3.0)
        plt.title("Mean GTP/GDP-Pi prior to catastrophe")
        plt.xlabel('time before catastrophe [s]');
        plt.ylabel('GTP/GDP-Pi');

        if simPa.record_data:
                filename = file_figure + '_fig' + str(int(num_fig))
                plt.savefig(filename+'.eps', format='eps', dpi=1000)
                plt.savefig(filename+'.png', format='png', dpi=200)

        plt.show()

    return EB_output


def fig_MT_before_cat(simPa,
                      file_figure,
                      num_fig,
                      MT_length_sum,
                      cap_end_sum):
    """ Figure showing tip and cap position before catastrope.

    Args:
    -------
    simPa: parameter set
        Simulation parameters in "ParameterSet" format.
    file_figure: str
        Folder for storing figures and data.
    num_fig: int
        Figure number.
    MT_length_sum:
        TODO
    cap_end_sum:
        TODO
    """
    # Remove invalid growth events from dataset
    index = np.sum(MT_length_sum, axis=1)
    index = np.where(index == 0)[0]
    MT_length = np.delete(MT_length_sum, index, 0)
    C_length = np.delete(cap_end_sum, index, 0)

    # Normalize length to max position
    MT_length_norm = (MT_length.T - np.max(MT_length, axis=1)).T
    C_length_norm = (C_length.T - np.max(C_length, axis=1)).T

    # Calculate mean and error MT tip position
    MT_length_mean = np.mean(MT_length_norm, axis=0)*(simPa.dL_dimer*1000)
    MT_length_mean = MT_length_mean + abs(np.max(MT_length_mean))
    MT_length_std = np.std(MT_length_norm, axis=0)

    # Calculate mean and error cap end position
    C_length_mean = np.mean(C_length_norm, axis=0)*(simPa.dL_dimer*1000)
    C_length_mean = C_length_mean + abs(np.max(C_length_mean))
    C_length_std = np.std(C_length_norm, axis=0)

    # Ensure that tip and cap coincide at moment of catastrophe
    C_length_mean = C_length_mean + MT_length_mean[-1]

    # Calculate 95% confidence intervals of tip position
    MT_CI_upper = MT_length_mean + (1.96 * (MT_length_std/np.sqrt(MT_length.shape[0])))
    MT_CI_lower = MT_length_mean - (1.96 * (MT_length_std/np.sqrt(MT_length.shape[0])))

#    MT_CI_upper = MT_length_mean + MT_length_std
#    MT_CI_lower = MT_length_mean - MT_length_std


    # Calculate 95% confidence intervals of cap end position
    C_CI_upper = C_length_mean + (1.96 * (C_length_std/np.sqrt(C_length.shape[0])))
    C_CI_lower = C_length_mean - (1.96 * (C_length_std/np.sqrt(C_length.shape[0])))

#    C_CI_upper = C_length_mean + C_length_std
#    C_CI_lower = C_length_mean - C_length_std

    # Calculate pausing duration prior to catastrophe (def: < 10% mean growth speed)
    MT_length_diff = np.diff(MT_length_mean) / simPa.frame_rate_actual
    Pausing = np.where(MT_length_diff < 0.1*simPa.growth_speed*1000/60)[0][0]
    Pausing = (MT_length.shape[1] - Pausing)*simPa.frame_rate_actual

    # Tip position output
    MT_length_output = np.vstack((MT_CI_upper, MT_CI_lower))
    MT_length_output = np.vstack((MT_length_mean, MT_length_output))
    MT_length_output = np.vstack((simPa.frame_rate_actual * np.arange(-len(MT_length_mean),0), MT_length_output)).T

    # Cap position output
    C_length_output = np.vstack((C_CI_upper, C_CI_lower))
    C_length_output = np.vstack((C_length_mean, C_length_output))
    C_length_output = np.vstack((simPa.frame_rate_actual * np.arange(-len(C_length_mean),0), C_length_output)).T

    #Start plotting
    if num_fig != 0:

        plt.figure(num_fig)
        plt.clf()

        # Plot mean tip position
        plt.plot(simPa.frame_rate_actual * np.arange(-len(MT_length_mean),0),
                 MT_length_mean,'k', linewidth=2.0)

        # Plot upper 95% confidence interval of tip position
        plt.plot(simPa.frame_rate_actual * np.arange(-len(MT_length_mean),0),
                 MT_CI_upper,'k--', linewidth=1.0)

        # Plot lower 95% confidence interval of tip position
        plt.plot(simPa.frame_rate_actual * np.arange(-len(MT_length_mean),0),
                 MT_CI_lower,'k--', linewidth=1.0)

        # Plot mean cap position
        plt.plot(simPa.frame_rate_actual * np.arange(-len(C_length_mean),0),
                 C_length_mean,'r', linewidth=2.0)

        # Plot upper 95% confidence interval of cap position
        plt.plot(simPa.frame_rate_actual * np.arange(-len(C_length_mean),0),
                 C_CI_upper,'r--', linewidth=1.0)

        # Plot lower 95% confidence interval of cap position
        plt.plot(simPa.frame_rate_actual * np.arange(-len(C_length_mean),0),
                 C_CI_lower,'r--', linewidth=1.0)

        plt.title("MT length prior to catastrophe")
        plt.xlabel('time before catastrophe [s]');
        plt.ylabel('Mean MT position [nm]');

        #Add parameters to figure
        figtext = ['Shrinkage = %.2f nm'  %abs(MT_length_mean[-1])]
        figtext.append('Pausing = %.2f s'  %Pausing)

        ## Figure styles
        from matplotlib.font_manager import FontProperties
        font = FontProperties()
        font.set_family('sans-serif')
        font.set_style('normal')
        font.set_weight('light')

        figDX = 0.045
        for m in range(len(figtext)):
            plt.ax.text(0.1,0.82-m*figDX, figtext[m], fontproperties=font,
                verticalalignment='bottom', horizontalalignment='left',
                transform=plt.ax.transAxes, color='black', fontsize=11)

        plt.xlim(-30, 2)

        if simPa.record_data:
                filename = file_figure + '_fig' + str(int(num_fig))
                plt.savefig(filename+'.eps', format='eps', dpi=1000)
                plt.savefig(filename+'.png', format='png', dpi=200)

        plt.show()

        print('Mean shrinkage before catastrophe = %.2f nm'  %abs(MT_length_mean[-1]))
        print('Mean pause duration before catastrophe = %.2f s'  %Pausing)

    return MT_length_output, C_length_output


def fig_EB_cat_hist(simPa,
                    file_figure,
                    num_fig,
                    EB_comet_sum,
                    barrier_contact_times,
                    EB_average_frames = 2):
    """ Have a look at EB intensity at catastrophe...

    Args:
    -------
    simPa: parameter set
        Simulation parameters in "ParameterSet" format.
    file_figure: str
        Folder for storing figures and data.
    num_fig: int
        Figure number.
    EB_comet_sum: list, array #TODO:check
        Number of "B"s during a time window before catastophe.
    barrier_contact_times: list, array #TODO:check
        List/array containing barrier contact times.
    EB_average_frames: int
        Number of frames to average over. Default = 2.
    """
    EB_intensity_before_cat = []
    EB_intensity_at_barrier = []
    EB_mean = []

    # Select valid runs
    valid_runs = valid_EB_runs(simPa, EB_comet_sum, barrier_contact_times)

    EB_signal, EB_signal_average, max_barrier_contact_frames, min_length_run_frames, frame_window = analyse_EB_signal(simPa,
                                                                                 EB_comet_sum, barrier_contact_times)

    for a in range(0,len(valid_runs)):
        EB_intensity_before_cat.append(np.mean(np.array(EB_comet_sum[a])[0:(EB_average_frames+1)])) # :-1]))
        barrier_contact_frame = int(round(barrier_contact_times[valid_runs[a]]/simPa.frame_rate_actual,0))
        EB_intensity_at_barrier.append(np.mean(np.array(EB_comet_sum[a])[barrier_contact_frame:(barrier_contact_frame+EB_average_frames+1)]))

        EB_mean.append(np.mean(EB_signal[a][max_barrier_contact_frames:frame_window]))

    fig, ax = plt.subplots(figsize=(8, 8)) #figure(9, figsize=(8, 8))
    plt.clf()
    map = plt.scatter(EB_intensity_before_cat/np.mean(EB_mean), EB_intensity_at_barrier/np.mean(EB_mean),
                c = barrier_contact_times[valid_runs], alpha=0.5, cmap='CMRmap')
    plt.xlim(xmax=1)
    fig.colorbar(map, ax = ax, label = 'barrier contact time [s]')

    plt.title('EB intensity before catastrophe (%.2f nM EB)' %(simPa.EB*1000))
    plt.xlabel('EB intensity right before catastrophe (last %.0f frames), relative to mean' %EB_average_frames)
    plt.ylabel('EB intensity right before barrier contact, relative to mean')
    plt.legend(fontsize=14)

    if simPa.record_data:
        filename = file_figure + '_fig' + str(num_fig) + '_relative'
        plt.savefig(filename+'.eps', format='eps', dpi=1000)
        plt.savefig(filename+'.png', format='png', dpi=200)

    plt.show()

    fig, ax = plt.subplots(figsize=(8, 8)) #figure(9, figsize=(8, 8))
    plt.clf()

    hist_data, hist_bins = np.histogram(EB_intensity_before_cat/np.mean(EB_mean), np.arange(0,1.1,0.1))
    bin_width = hist_bins[1]

    plt.bar((hist_bins[:-1] + bin_width/2) , np.float_(hist_data)/(np.sum(hist_data)), 0.9*bin_width, alpha=0.8)
    plt.title('Relative B-state intensity at catastrophe')
    plt.xlabel('Relative B-state intensity (#of elements in state "B" div. by mean)')
    plt.ylabel('Probability')

    if simPa.record_data:
        filename = file_figure + '_fig' + str(num_fig) + '_histogram'
        plt.savefig(filename+'.eps', format='eps', dpi=1000)
        plt.savefig(filename+'.png', format='png', dpi=200)

    #PROBLEM with plt.hist --> normed=1 and density=1 don't work properly
    #plt.hist(EB_intensity_before_cat/np.mean(EB_mean), np.arange(0,1.1,0.1), density=True, histtype='bar', rwidth=0.8)
    plt.show()

    return EB_intensity_before_cat/np.mean(EB_mean)


def fig_display_examples(simPa,
                         file_figure,
                         num_fig,
                         MT_length_sum,
                         catastrophe_times,
                         EB_comet_sum,
                         barrier_contact_times=[]):
    """ Show selection of examples (tip position + EB intensity)

    Args:
    -------
    simPa: parameter set
        Simulation parameters in "ParameterSet" format.
    file_figure: str
        Folder for storing figures and data.
    num_fig: int
        Figure number.
    MT_length_sum:
        TODO
    catastrophe_times:
        TODO
    EB_comet_sum: list, array #TODO:check
        Number of "B"s during a time window before catastophe.
    barrier_contact_times: list, array #TODO:check
        List/array containing barrier contact times.
    """
    min_length_run_frames = int(simPa.min_length_run/simPa.frame_rate_actual)

    # Select valid runs
    valid_runs = valid_EB_runs(simPa, EB_comet_sum, barrier_contact_times)

    EB_signal_before_cat = np.zeros((len(valid_runs), min_length_run_frames+1)) #put individual runs into one np.array
    for a in range(0,len(valid_runs)):
        EB_signal_before_cat[a][0:min_length_run_frames] = \
        np.array(EB_comet_sum[valid_runs[a]])[0:min_length_run_frames]

    EB_signal_average = np.sum(EB_signal_before_cat, axis=0)
    EB_signal_average = EB_signal_average/len(valid_runs)

    show_fraction_frames = int(simPa.show_fraction/simPa.frame_rate_actual)
    valid_runs = np.where(catastrophe_times > simPa.show_fraction)[0]

    plt.figure(num_fig, figsize=(15, 10))
    plt.clf()
    f, axarr = plt.subplots(nrows=5, ncols=5, sharey=True, sharex=True, figsize=(15, 10))
    for m in range(0,5):
        for n in range(0,5):
            skip = 0
            axarr[m, n].plot(simPa.frame_rate_actual * np.arange(-min_length_run_frames,0) ,MT_length_sum[valid_runs[skip+m+5*n]][0::], 'black')
            axarr[m, n].set_title('catastrophe %.0f' %(skip+m+5*n))
            axarr[m, n]
            plt.plot(simPa.frame_rate_actual * np.arange(-min_length_run_frames,0) ,EB_signal_before_cat[skip+m+5*n][0:show_fraction_frames][::-1],'red')

    if simPa.record_data:
            filename = file_figure + '_fig' + str(int(num_fig))
            plt.savefig(filename+'.eps', format='eps', dpi=1000)
            plt.savefig(filename+'.png', format='png', dpi=200)

    plt.show()


def fig_EB_profile(simPa, file_figure, num_fig, EB_profiles, MT_length_full, w_size):
    """ Figure to display EB profile.

    Args:
    -------
    simPa: parameter set
        Simulation parameters in "ParameterSet" format.
    file_figure: str
        Folder for storing figures and data.
    num_fig: int
        Figure number.
    EB_profiles:
        TODO
    MT_length_full:
        TODO
    w_size:
        TODO
    """
    # Analyse the simulated EB profiles
    EB_mean, v_mean = analyse_EB_profile(simPa, MT_length_full, EB_profiles, w_size)

    # Calculate the mean EB profile
#    x = np.arange(0, len(EB_mean[0]), 1) * simPa.dL_dimer *1000
#    y = np.mean(EB_mean, axis=0)

    x = np.arange(0, len(EB_mean), 1) * simPa.dL_dimer *1000
    y = EB_mean
    #plt.plot(x, y)

    # Define exponential function
    def exponenial_func(x, a, b):
        return a*np.exp(b*x)

    # Calculate the maturation rate (Duellberg, 2016), i.e. the hydrolysis rate
    ind = np.argmax(y)
    popt, pcov = curve_fit(exponenial_func, x[0:ind], y[0:ind], p0=(1e-2, 1e-3))

    xx= np.linspace(0, x[ind], 1000)
    yy = exponenial_func(xx, *popt)

    # Align profile from left to right and set tip position to zero
    x = -x
    x += np.round(np.argmax(y)*(simPa.dL_dimer*1000))
    xx = -xx
    xx += np.round(np.argmax(y)*(simPa.dL_dimer*1000))

    fig = plt.figure(1, figsize=(12, 7))
    plt.clf()

    plt.plot(x,y,'k.', xx, yy, '--r')

    plt.title('Mean GTP/GDP-Pi profile', fontsize=14)
    plt.xlabel('Position [nm]', fontsize=12)
    plt.ylabel('Intensity [a.u.]', fontsize = 12)

    plt.ax = fig.add_subplot(111)
    fig.subplots_adjust(top=0.9)

    ## Figure styles
    from matplotlib.font_manager import FontProperties
    font = FontProperties()
    font.set_family('sans-serif')
    font.set_style('normal')
    font.set_weight('light')

    #Add parameters to figure
    figtext = ['Simulation parameters:']
    figtext.append('$N_{sim} = %.0f$' %len(MT_length_full))
    figtext.append('$v_{g} = %.2f$ $nm/s$' %float(simPa.growth_rate_one*(simPa.dL_dimer*1000)))
    figtext.append('$k_{hyd} = %.2f$ $s^{-1}$' %(simPa.kBC))
    figtext.append('$D_{tip} = %.0f$ $nm^{2}/s$' %simPa.D_tip)
    figtext.append('')
    figtext.append('Measured values:')
    figtext.append('$N_{profiles} = %.0f$' %len(EB_mean))
    figtext.append('$v_{g} = %.2f$ $nm/s$' %np.mean(v_mean))
    figtext.append('$L_{comet} = %.0f$ $nm$' %float(1/popt[1]))
    figtext.append('$k_{m} = %.2f$ $s^{-1}$' %float(np.mean(v_mean)*popt[1]))
    figtext.append('$I_{max} = %.2f$ $a.u.$' %float(np.max(y)))

    figDX = 0.045
    for m in range(len(figtext)):
        plt.ax.text(0.75, 0.9-m*figDX, figtext[m], fontproperties=font,
            verticalalignment='bottom', horizontalalignment='left',
            transform=plt.ax.transAxes, color='black', fontsize=12)

    if simPa.record_data:
            filename = file_figure + '_fig' + str(int(num_fig))
            plt.savefig(filename+'.eps', format='eps', dpi=1000)
            plt.savefig(filename+'.png', format='png', dpi=200)

    plt.show()

def fig_MT_ageing(simPa, file_figure, num_fig, c_times):
    """ Calculate the age-dependent microtubule catastrophe frequency

    Args:
    -------
    simPa: parameter set
        Simulation parameters in "ParameterSet" format.
    file_figure: str
        Folder for storing figures and data.
    num_fig: int
        Figure number.
    c_times: numpy.array
        Array of catastrophe times.
    """

    X_dist = np.sort(c_times)
    Y_dist = np.cumsum(np.ones(len(c_times)))/len(c_times)

    f = (Y_dist/X_dist)/(1-Y_dist)

    C_freq = np.vstack((X_dist,Y_dist)).T

    plt.figure(1, figsize=(12, 7))
    plt.clf()

    plt.plot(X_dist, f)
    plt.ylim(0, 0.02)
    plt.xlim(0, 400)
    plt.title('Microtubule ageing', fontsize=14)
    plt.xlabel('Microtubule age [s]', fontsize=12)
    plt.ylabel('Catastrophe frequency [$s^{-1}$]', fontsize = 12)

    if simPa.record_data:
            filename = file_figure + '_fig' + str(int(num_fig))
            plt.savefig(filename+'.eps', format='eps', dpi=1000)
            plt.savefig(filename+'.png', format='png', dpi=200)
    plt.show()

    return C_freq


def fig_dist_fit(simPa,
                 file_figure,
                 num_fig,
                 Cum_dist):
    """

    Args:
    -------
    simPa: parameter set
        Simulation parameters in "ParameterSet" format.
    file_figure: str
        Folder for storing figures and data.
    num_fig: int
        Figure number.
    Cum_dist:
        TODO
    """

    fig, ax1 = plt.subplots(figsize=(12, 7))
    plt.clf()

    if isinstance(Cum_dist, list) and list_dim(Cum_dist) > 1 and list_dim(Cum_dist[0]) > 0:
        if isinstance(Cum_dist[0], np.ndarray):
            print(list_dim(Cum_dist), ' different cumulative distributions found. ')
        else:
            print('Error: Input cumulative distributions must be numpy arrays or lists of numpy arrays.' )
    elif isinstance(Cum_dist, list) and list_dim(Cum_dist) == 1 and isinstance(Cum_dist[0], np.ndarray):
        pass;
    elif isinstance(Cum_dist, np.ndarray):
       Cum_dist = [Cum_dist] #put numpy array into list
    else:
        print('Error: Input cumulative distributions must be numpy arrays or lists of numpy arrays.' )


    x = Cum_dist[0]
    x_fit = np.linspace(0,x[-1],1000)
    y = np.linspace(0,1,len(x))

    # Fit cumulative distribution to the Gamma function
    popt1, pcov1 = curve_fit(gamma_cdf,x, y, p0=(1, 1e-2))
    y1 = gamma_cdf(x_fit, *popt1)
    print(popt1)
    print(pcov1)

    # Fit cumulative distribution to an exponential
    popt2, pcov2 = curve_fit(exp_cdf, x, y, p0=(1e-2))
    y2 = exp_cdf(x_fit, *popt2)

    if list_dim(Cum_dist) > 1:
        c_range = 1/(list_dim(Cum_dist)-1)
    else:
       c_range = 1

    for i, Cum_dist in enumerate(Cum_dist):
        plt.step(Cum_dist, 1/(len(Cum_dist)-1) * np.arange(0, len(Cum_dist) , 1),
                 where='post', color=(0.95-0.7*(i)*c_range, 0.1, 0.1 + 0.8*(i)*c_range), linewidth=1.5, label='Simulation')

    plt.plot(x_fit, y1, 'k--', linewidth=1.5, label='Gamma fit')
    plt.plot(x_fit, y2, 'k:', linewidth=1.5, label='Exponential fit')

    plt.title('Microtubule lifetime distribution', fontsize=14)
    plt.xlabel('time [s]')
    plt.ylabel('Cumulative fraction')

    plt.ax = fig.add_subplot(111)
    figtext = ['Simulation parameters:']
    figtext.append('$v_{g} = %.1f$ $nm/s$' %float(simPa.growth_rate_one*(simPa.dL_dimer*1000)))
    figtext.append('$k_{hyd} = %.2f$ $s^{-1}$' %(simPa.kBC))
    if simPa.D_tip_time:
        figtext.append('$k_{ageing} = %.3f$ $s^{-1}$' %(simPa.D_tip_rate_T))
    elif simPa.D_tip_length:
        figtext.append('$k_{ageing} = %.3f$ $s^{-1}$' %(simPa.D_tip_rate_L))

    figtext.append('Gamma fit parameters:')
    figtext.append('$steps = %.2f$' %popt1[0])
    figtext.append('$rate = %.3f$ $s^{-1}$' %popt1[1])

    figDX = 0.045
    for m in range(len(figtext)):
        plt.ax.text(0.7, 0.82-m*figDX, figtext[m], fontproperties=font,
            verticalalignment='bottom', horizontalalignment='left',
            transform=plt.ax.transAxes, color='black', fontsize=11)

    plt.legend(fontsize=12)

    # Plot insert with microtubule ageing parameters
    left, bottom, width, height = [0.6, 0.21, 0.25, 0.25]
    ax2 = fig.add_axes([left, bottom, width, height])

    # Calculate the evolution of tip diffusion
    D_tip = []
    if simPa.D_tip_time:
        for i in range(len(x_fit)):
            D_tip.append((simPa.D_tip_end - simPa.D_tip_start) * (1 - np.exp(-1*simPa.D_tip_rate_T*(x_fit[i]))) + simPa.D_tip_start)
    elif simPa.D_tip_length:
        x_fit = np.linspace(0,50000,50000)*simPa.dL_dimer
        for i in range(len(x_fit)):
            D_tip.append((simPa.D_tip_end - simPa.D_tip_start) * (1 - np.exp(-1*simPa.D_tip_rate_L*(x_fit[i]))) + simPa.D_tip_start)
    else:
        D_tip = np.ones(len(x_fit))*simPa.D_tip

    # Plot tip diffusion
    ax2.plot(x_fit, D_tip, 'k')
    ax2.set_ylim(0, 1.1*np.max(D_tip))
    if simPa.D_tip_time:
        ax2.set_xlabel('time [s]')
    elif simPa.D_tip_length:
        ax2.set_xlabel('length [um]')
    ax2.set_ylabel('Tip diffusion [$nm^{2}s^{-1}$]')
    ax2.set_title('Microtubule ageing')

    if simPa.record_data:
        filename = file_figure + '_fig' + str(int(num_fig))
        plt.savefig(filename+'.eps', format='eps', dpi=1000)
        plt.savefig(filename+'.png', format='png', dpi=200)

    plt.show()


def fig_nucleation(simPa,
                   file_figure,
                   num_fig,
                   nucleation_times):
    """ Figure showing nucleation time histogram.

    Args:
    -------
    simPa: parameter set
        Simulation parameters in "ParameterSet" format.
    file_figure: str
        Folder for storing figures and data.
    num_fig: int
        Figure number.
    nucleation_times: list, array
        List, array of nucleation times.
    """
    plt.subplots(figsize=(12, 7))
    plt.clf()
    plt.hist(nucleation_times, bins='auto')
    plt.xlabel('time [s]')
    plt.ylabel('Counts')
    plt.title('Nucleation time', fontsize=14)
    plt.show()


def fig_washout(simPa,
                file_figure,
                num_fig,
                washout_times,
                catastrophe_washout,
                MT_length_sum):
    """ Figure showing MT behavior for washout experiment.

    Args:
    -------
    simPa: parameter set
        Simulation parameters in "ParameterSet" format.
    file_figure: str
        Folder for storing figures and data.
    num_fig: int
        Figure number.
    washout_times: numpy.array, list
        List/array of washout time points.
    catastrophe_washout:
        TODO: Maurits --> add description
    """
    # Plot the washout delay times
    fig, axs = plt.subplots(1,3, figsize=(14,7))
    axs[0].hist(washout_times, bins='auto')
    axs[0].set_xlabel('delay time [s]')
    axs[0].set_ylabel('counts')
    axs[0].set_xlim(left=0)

    # Calculate the shrinkage length
    shrinkage_length = []
    delay_steps = np.round(washout_times / simPa.dt)
    delay_steps = delay_steps.astype(np.int16)
    for i in range(len(catastrophe_washout)):
        shrinkage_length.append((MT_length_sum[catastrophe_washout[i],-delay_steps[i]] - MT_length_sum[catastrophe_washout[i],-1])*simPa.dL_dimer*1000)

    # Plot the shrinkage length
    axs[1].hist(shrinkage_length, bins='auto')
    axs[1].set_xlabel('shrinkage length [nm]')
    axs[1].set_ylabel('counts')
    axs[1].set_xlim(left=0)

    # Calculate the shrinkage speeds
    axs[2].hist(shrinkage_length/washout_times,bins='auto')

    # Plot the shrinkage speed
    axs[2].set_xlabel('shrinkage speed [nm/s]')
    axs[2].set_ylabel('counts')

    if simPa.record_data:
        filename = file_figure + '_fig' + str(int(num_fig))
        plt.savefig(filename+'.eps', format='eps', dpi=1000)
        plt.savefig(filename+'.png', format='png', dpi=200)


def fig_cap_size_before_cat(simPa,
                            file_figure,
                            num_fig,
                            MT_length_full,
                            cap_end):
    """ Figure showing cap size before catastrophe.

    Args:
    -------
    simPa: parameter set
        Simulation parameters in "ParameterSet" format.
    file_figure: str
        Folder for storing figures and data.
    num_fig: int
        Figure number.
    MT_length_full:
        TODO
    cap_end:
        TODO
    """
    min_length_run_frames = int(simPa.min_length_run/simPa.frame_rate_actual)

    MT_run_length = []
    for i in range(len(cap_end)):
        MT_run_length.append(len(cap_end[i]))

    valid_runs = np.where(np.array(MT_run_length) > min_length_run_frames)[0]

    cap_before_cat = np.zeros((len(valid_runs), min_length_run_frames+1)) #put individual runs into one np.array
    for a in range(0,len(valid_runs)):
        cap_before_cat[a][0:min_length_run_frames] = \
        np.array(MT_length_full[valid_runs[a]])[-min_length_run_frames:] - np.array(cap_end[valid_runs[a]])[-min_length_run_frames:]

    cap_average = np.sum(cap_before_cat, axis=0)
    cap_average = cap_average/len(cap_end)
    cap_average = cap_average*simPa.dL_dimer*1000
    cap_before_cat = cap_before_cat*simPa.dL_dimer*1000

    plt.figure(num_fig)
    plt.clf()

    plt.plot(simPa.frame_rate_actual * np.arange(-min_length_run_frames,0),
             cap_average[0:min_length_run_frames][::1],'red', linewidth=2.0)

    for a in range(0,len(valid_runs)):
        plt.plot(simPa.frame_rate_actual * np.arange(-min_length_run_frames,0),
                 cap_before_cat[a][0:min_length_run_frames][::1],color=plt.cm.Reds(0.3+0.7*a/len(valid_runs)))
        plt.plot(simPa.frame_rate_actual * np.arange(-min_length_run_frames,0),
                cap_average[0:min_length_run_frames][::1],'black', linewidth=3.0)
    plt.title("Mean Cap size prior to catastrophe")
    plt.xlabel('time before catastrophe [s]');
    plt.ylabel('Length of cap [nm]');

    if simPa.record_data:
            filename = file_figure + '_fig' + str(int(num_fig))
            plt.savefig(filename+'.eps', format='eps', dpi=1000)
            plt.savefig(filename+'.png', format='png', dpi=200)

    plt.show()
