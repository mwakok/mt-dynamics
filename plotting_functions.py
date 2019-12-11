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
import random
from time import clock #, strftime
from matplotlib import pyplot as plt
import pandas as pd 
from scipy.interpolate import interp1d
from scipy.special import gammainc
from itertools import islice

# Initialize random number generator:
np.random.seed(int(100*clock()))

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
    
    ax1.hold(True)
    ax1.hist(v_fluc, bins = 60, density = True, color = "skyblue", label = "simulated data")
    ax1.plot(x, G, 'r', label = "theoretical distribution")
    ax1.hold(False)
    
    ax2.hist(c_fluc, bins = 60, density = True, color = "skyblue", label = "simulated data")

    move =  len(c_fluc)/sample_size
    pause = 1 - move
    step_mean = np.mean(c_fluc)
    step_std =  np.std(c_fluc)              
    
    if simPa.record_data:        
            filename = file_figure + '_fig' + str(int(num_fig))
            plt.savefig(filename+'.eps', format='eps', dpi=1000)
            plt.savefig(filename+'.png', format='png', dpi=200) 
            
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
        
        if len(Cum_dist_compare) > 0:     
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
        file_csv = file_figure[:-10] + "EB" + str(simPa.EB*1000) + "_cumhist.csv"
        Cum_dist_pd = pd.DataFrame(np.round(Cum_dist,2))
        Cum_dist_pd.to_csv(file_csv, header=None, index=None)

    plt.show()  
    

def fig_EB_at_barrier(simPa, file_figure, num_fig, EB_comet_sum, barrier_contact_times):
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
        Barrier contact times.
    """
    
    # Select valid runs
    valid_runs = valid_EB_runs(simPa, EB_comet_sum, barrier_contact_times)
    
    EB_signal, EB_signal_average, max_barrier_contact_frames, min_length_run_frames, frame_window = analyse_EB_signal(simPa, 
                                                                                                             EB_comet_sum, 
                                                                                                             barrier_contact_times)
    
    plt.figure(num_fig)
    plt.clf()
    for a in range(0,len(valid_runs)):
        plt.plot(simPa.frame_rate_actual * np.arange(-min_length_run_frames+1,max_barrier_contact_frames+1), 
                 EB_signal[a][0:frame_window][::-1],color=plt.cm.Reds(0.3+0.7*a/len(valid_runs))) 
    
        plt.plot(simPa.frame_rate_actual * np.arange(-min_length_run_frames+1,max_barrier_contact_frames+1),
             EB_signal_average[0:frame_window][::-1],'black', linewidth=3.0) 
    plt.title("number of B's")
    plt.xlabel('time before barrier contact [s]');
    plt.ylabel('EB comet intensity');
    
    from scipy.optimize import curve_fit
    def func(x, a, b, c):
        return a * np.exp(-b * x) + c
#   alternative tests...
#    def func(x, a, b, c, d, e):
#        return a * np.exp(-b * x) + c*np.exp(-d * x) + e
#    def func(x, a, b, c, d, e):
#        return a * np.exp(-b * x) + c * x + d
    xdata = simPa.frame_rate_actual * np.arange(0,int(max_barrier_contact_frames/2.5))
    ydata = EB_signal_average[(max_barrier_contact_frames-int(max_barrier_contact_frames/2.5)):max_barrier_contact_frames][::-1]
#    
#    def func(x, a, b, c, d, e):
#        return a * np.exp(-b * x) + d * np.exp(-e * x)+ c
#    xdata = simPa.frame_rate_actual * np.arange(0,int(max_barrier_contact_frames/3))
#    ydata = EB_signal_average[(max_barrier_contact_frames-int(max_barrier_contact_frames/3)):max_barrier_contact_frames][::-1]
##    
    popt, pcov = curve_fit(func, xdata, ydata, p0=(np.max(ydata), simPa.kBC, 1),maxfev=1000)
    print(popt)
    plt.plot(xdata, func(xdata, *popt), 'c--',)
    
    plt.text(2,2*max(ydata),'decay rate (exp. fit): %.2f' %popt[1])
    
    if simPa.record_data:        
            filename = file_figure + '_fig' + str(int(num_fig))
            plt.savefig(filename+'.eps', format='eps', dpi=1000)
            plt.savefig(filename+'.png', format='png', dpi=200) 
    
    plt.show()
    
    
def fig_EB_before_cat(simPa, file_figure, num_fig, EB_comet_sum, barrier_contact_times=[]):
    ## Plot EB intensity (here = elements in state "B") before catastrophe  
    # simPa - simulation parameters in "ParameterSet" format
    # file_figure - folder for storing figures and data
    # num_fig - figure number
    # EB_comet_sum - 
    # barrier_contact_times - needed if barrier not False

#    min_length_run = 2* simPa.min_length_run # only consider runs longer than xx frames
    min_length_run_frames = int(simPa.min_length_run/simPa.frame_rate_actual)  
        
    # Select valid runs
    valid_runs = valid_EB_runs(simPa, EB_comet_sum, barrier_contact_times)
    
    EB_signal_before_cat = np.zeros((len(valid_runs), min_length_run_frames+1)) #put individual runs into one np.array
    for a in range(0,len(valid_runs)):
        EB_signal_before_cat[a][0:min_length_run_frames] = \
        np.array(EB_comet_sum[valid_runs[a]])[0:min_length_run_frames]
        
    EB_signal_average = np.sum(EB_signal_before_cat, axis=0)
    EB_signal_average = EB_signal_average/len(valid_runs)
    
    plt.figure(num_fig)
    plt.clf()
    #simPa.frame_rate_actual * np.arange(-min_length_run_frames,0)
    plt.plot(simPa.frame_rate_actual * np.arange(-min_length_run_frames,0), 
             EB_signal_average[0:min_length_run_frames][::-1],'red', linewidth=2.0) 
    
    #plt.plot(range(-simPa.show_fraction+1,1) ,EB_signal_average[0:simPa.show_fraction][::-1],'red', linewidth=2.0) 
    for a in range(0,len(valid_runs)):
        plt.plot(simPa.frame_rate_actual * np.arange(-min_length_run_frames,0), 
                 EB_signal_before_cat[a][0:min_length_run_frames][::-1],color=plt.cm.Reds(0.3+0.7*a/len(valid_runs))) 
        plt.plot(simPa.frame_rate_actual * np.arange(-min_length_run_frames,0), 
                 EB_signal_average[0:min_length_run_frames][::-1],'black', linewidth=3.0)
    plt.title("number of B's")
    plt.xlabel('time before catastrophe [s]');
    plt.ylabel('EB comet intensity');
    
    if simPa.record_data:        
            filename = file_figure + '_fig' + str(int(num_fig))
            plt.savefig(filename+'.eps', format='eps', dpi=1000)
            plt.savefig(filename+'.png', format='png', dpi=200) 
    
    plt.show()
    
    
def fig_EB_cat_hist(simPa, file_figure, num_fig, EB_comet_sum, barrier_contact_times, EB_average_frames = 2):
    ## Have a look at EB intensity at catastrophe...
    # simPa - simulation parameters in "ParameterSet" format
    # file_figure - folder for storing figures and data
    # num_fig - figure number
    # EB_comet_sum - num of "B"s during a time window before catastophhe
    # barrier_contact_times
    #
    # EB_average_frames = 2 (average over X frames to get B intensity at catastrophe)
    
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
    #plt.ylabel(' No. of dimers in "B-state" before catastrophe (averaged over %.2f s)' %(EB_average_frames * simPa.frame_rate_actual))
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
    
    
def fig_display_examples(simPa, file_figure, num_fig, MT_length_sum, catastrophe_times, EB_comet_sum, barrier_contact_times=[]):
    ## Show selection of examples (tip position + EB intensity)
    # simPa - simulation parameters in "ParameterSet" format
    # file_figure - folder for storing figures and data
    # num_fig - figure number
    
#    # only consider runs longer than xx frames
#    
#    # only consider runs longer than xx frames
#    min_length_run = 2* simPa.min_length_run # only consider runs longer than xx frames
#    #if FACT.show_fraction > min_length_run:
#    #    FACT.show_fraction = min_length_run
#    min_length_run_frames = int(simPa.min_length_run/simPa.frame_rate_actual)  
#        
#    b = []
#    if simPa.barrier: #if barrier present  
#        for a in range(0,len(EB_comet_sum)):
#            b.append(len(EB_comet_sum[a]) * simPa.frame_rate_actual - barrier_contact_times[a])
#    else:
#        for a in range(0,len(EB_comet_sum)):
#            b.append(len(EB_comet_sum[a]) * simPa.frame_rate_actual)
#            
#    test = np.where(np.array(b) > min_length_run)[0]  
#    
#    
#    EB_signal_before_cat = np.zeros((len(test), min_length_run_frames+1)) #put individual runs into one np.array
#    for a in range(0,len(test)):
#        EB_signal_before_cat[a][0:min_length_run_frames] = \
#        np.array(EB_comet_sum[test[a]])[0:min_length_run_frames]
    
    
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
    MT_length_average = np.sum(MT_length_sum[valid_runs][0::], axis = 0)/len(valid_runs)
    
    plt.figure(num_fig, figsize=(15, 10))
    plt.clf()
    f, axarr = plt.subplots(nrows=5, ncols=5, sharey=True, sharex=True, figsize=(15, 10))
    axarr2 = []
    for m in range(0,5):
        for n in range(0,5):
            axarr[m, n].plot(simPa.frame_rate_actual * np.arange(-min_length_run_frames,0) ,MT_length_sum[valid_runs[m+5*n]][0::], 'black')
            axarr[m, n].set_title('catastrophe %.2f' %float(m+5*n))
            axarr[m, n]
            axarr2 = axarr[m,n].twinx()
            plt.plot(simPa.frame_rate_actual * np.arange(-min_length_run_frames,0) ,EB_signal_before_cat[m+5*n][0:show_fraction_frames][::-1],'red') 
    #
    #
    #fig, ax1 = plt.subplots()
    #ax1.plot(range(-show_fraction+1,1) ,MT_length_average,'Black', linewidth=2.0) 
    #for a in range(0,len(b)):
    #    ax1.plot(range(-show_fraction+1,1) ,MT_length_sum[b[a]][0::], color=plt.cm.Greys(0.3+0.7*a/len(b))) 
    ## plt.title('MT length before catastrophe')
    #ax1.set_xlabel('time before catastrophe [s]');
    #ax1.set_ylabel('MT length');
    #
    #ax2 = ax1.twinx()
    #ax2.plot(range(-show_fraction+1,1) ,EB_signal_average[0:show_fraction][::-1],'red', linewidth=2.0) 
    if simPa.record_data:        
            filename = file_figure + '_fig' + str(int(num_fig))
            plt.savefig(filename+'.eps', format='eps', dpi=1000)
            plt.savefig(filename+'.png', format='png', dpi=200) 
            
    plt.show()