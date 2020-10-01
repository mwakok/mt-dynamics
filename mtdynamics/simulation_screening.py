from parameters import ParameterSet
import numpy as np
from matplotlib import pyplot as plt

# Add multi core parallelization
from concurrent.futures import ThreadPoolExecutor, as_completed

# Import main simulation function
import simulation_main as sMAIN #import main simulation function

           
class simulation_screen(object):
    """
    Run screen of simulations over 1 to 3 parameters 
    
    Parameters
    ----------
    sim_parameters: dict
        Dictionary including all simulation parameters
    param1_name: str
        Name of parameter 1 to screen (use name according to sim_parameter naming).
    params1: list
        List of values to use for screening of parameter 1.
    param2_name: str, optional
        Name of parameter 2 to screen (use name according to sim_parameter naming).
    params1: list, optional
        List of values to use for screening of parameter 2.
    param3_name: str, optional
        Name of parameter 3 to screen (use name according to sim_parameter naming).
    params3: list, optional
        List of values to use for screening of parameter 3.
    """
    def __init__(self, sim_parameters, param1_name, params1, 
                 param2_name=None, params2=[0], 
                 param3_name=None, params3=[0]):
        self.param1_name = param1_name
        self.param2_name = param2_name
        self.param3_name = param3_name
        self.params1 = params1
        self.params2 = params2
        self.params3 = params3
        self.sim_parameters = sim_parameters
        self.dist_compare_cat = []
        self.dist_compare_contact = []

        # Initialize output variables: 
        dim1 = len(self.params1)
        dim2 = len(self.params2)
        dim3 = len(self.params3)
        # Output as numpy arrays: 
        self.SCAN_dt = np.zeros((dim1,dim2,dim3))
        self.SCAN_mean_cat_times = np.zeros((dim1,dim2,dim3))
        self.SCAN_mean_contact_times = np.zeros((dim1,dim2, dim3))
        # Output as lists:
        self.SCAN_EB_comet_sum = [[[0]*dim3 for i in range(dim2)] for j in range(dim1)]
        self.SCAN_EB_comet_sum_barrier = [[[0]*dim3 for i in range(dim2)] for j in range(dim1)]
        self.SCAN_cat_times = [[[0]*dim3 for i in range(dim2)] for j in range(dim1)]
        self.SCAN_cat_length = [[[0]*dim3 for i in range(dim2)] for j in range(dim1)]
        self.SCAN_barrier_contact_times = [[[0]*dim3 for i in range(dim2)] for j in range(dim1)]
        self.SCAN_cap_end = [[[0]*dim3 for i in range(dim2)] for j in range(dim1)]
        
        # Initialize screening parameter sets: 
        parameter_collection = []
        for i, param1 in enumerate(self.params1):
            if self.param2_name == None: # screen 1 parameter
                sim_parameters = self.sim_parameters
                sim_parameters[self.param1_name] = param1
                print(sim_parameters)
                parameter_collection.append([i, 0, 0, sim_parameters.copy()])
            elif self.param3_name == None: # screen 2 parameters
                for j, param2 in enumerate(self.params2):
                    sim_parameters = self.sim_parameters
                    sim_parameters[self.param1_name] = param1
                    sim_parameters[self.param2_name] = param2  
                    parameter_collection.append([i, j, 0, sim_parameters.copy()])
            else: # screen 3 parameter
                for j, param2 in enumerate(self.params2):
                    for k, param3 in enumerate(self.params3):           
                        sim_parameters = self.sim_parameters
                        sim_parameters[self.param1_name] = param1
                        sim_parameters[self.param2_name] = param2
                        sim_parameters[self.param3_name] = param3  
                        parameter_collection.append([i, j, k, sim_parameters.copy()])

        self.parameter_collection = parameter_collection
        


    def run_screen_multicore(self, barrier_or_free = 'both', num_workers=4):
        """
        Run simulation parameter scan on multiple threads.

        Parameters
        ----------
        barrier_or_free: str
            Must be "barrier', 'free', or 'both'.
        num_workers: int
            Number of threads for parallelization.
            Usually the number of CPUs or 
        """
        
        #Create a pool of processes. By default, one is created for each CPU in your machine.
        #with concurrent.futures.ProcessPoolExecutor() as executor:
        print("Run parameter scan on ", num_workers, "number of workers.")
        with ThreadPoolExecutor(max_workers=num_workers) as executor:
            futures = [executor.submit(self.single_run, X, barrier_or_free=barrier_or_free) for X in self.parameter_collection]
                   
        for future in as_completed(futures):
            print("Single simulation run done.")
            #print(future.result())


    def run_screen(self, barrier_or_free = 'both'):
        """
        Run simulation parameter scan without parallelization.

        Parameters
        ----------
        barrier_or_free: str
            Must be "barrier', 'free', or 'both'.
        """
        
        print("Run parameter scan without parallelization.")
        for X in self.parameter_collection:
            self.single_run(X, barrier_or_free=barrier_or_free)




    def single_run(self, parameter_collect, barrier_or_free = 'both'):
        """
        Run simulations for one single parameter set

        Parameters
        ----------
        parameter_collect: list
            List including the position along the screening parameters,
            and the respective parameter set 
            --> [i,j,k, sim_parameters] 
            with sim_parameters a dictionary containing the parameters
        barrier_or_free: str
            Must be "barrier', 'free', or 'both'.
        """
        
        i, j, k, sim_parameters = parameter_collect    
        print('\n', '==='*5, 'Start run ', i,' - ',j, ' - ', k, ) 
        
        if (barrier_or_free == 'free') or (barrier_or_free == 'both'): 
            if sim_parameters['barrier']:
                barrier = sim_parameters['barrier']
            else:
                print("No barrier distance given in parameter file")
        
        if (barrier_or_free == 'free') or (barrier_or_free == 'both'): 
            sim_parameters['barrier'] = False
            simPa = ParameterSet(sim_parameters)
            
            dt, MT_length_sum, CATASTROPHE_TIMES, CATASTROPHE_LENGTH, barrier_contact_times, \
            EB_comet_sum, cap_end, frame_rate_actual, EB_profiles  = sMAIN.MT_RUN(simPa)
        
            self.SCAN_dt[i,j,k] = dt
            self.SCAN_EB_comet_sum[i][j][k] = EB_comet_sum
            self.SCAN_cat_times[i][j][k] = CATASTROPHE_TIMES
            self.SCAN_cat_length[i][j][k] = CATASTROPHE_LENGTH
            if len(CATASTROPHE_TIMES) > 10: # needs to have minimum 10 events to count !    
                tau_c = np.mean(np.array(CATASTROPHE_TIMES)) 
                print(len(CATASTROPHE_TIMES),' catastrophes | catastrophe time: %.2f s' %tau_c)
            else:
                tau_c = 2000
            self.SCAN_mean_cat_times[i,j,k] = tau_c    
            
        if (barrier_or_free == 'barrier') or (barrier_or_free == 'both'): 
            print('barrier at: ', barrier)
            
            sim_parameters['barrier'] = barrier
            simPa = ParameterSet(sim_parameters)
            dt, MT_length_sum, CATASTROPHE_TIMES, CATASTROPHE_LENGTH, barrier_contact_times, \
            EB_comet_sum, cap_end, frame_rate_actual, EB_profiles  = sMAIN.MT_RUN(simPa)
            
            if (barrier_or_free == 'barrier'):
                self.SCAN_dt[i,j,k] = dt
            self.SCAN_barrier_contact_times[i][j][k] = barrier_contact_times
        
            if barrier_contact_times.shape[0] > 0:     
                mean_contact = np.mean(barrier_contact_times) 
                print('Barrier contact time: %.2f s' %mean_contact)
            else:
                mean_contact = 0
            
            self.SCAN_mean_contact_times[i,j,k] = mean_contact
            self.SCAN_EB_comet_sum_barrier[i][j][k] = EB_comet_sum    

           
    def compare_to_results(self, cum_hist_exp, type = 'free'):
        """ Compare results from simulation screening to empirical data
        
        Parameters
        ----------
        cum_hist_exp: numpy array
            Cumulative distribution of contact or catastrophe times.
        type: str
            Type of empirical data, freely growing (type='free') or stalled ('barrier').
        """
        dim1 = len(self.params1)
        dim2 = len(self.params2)
        dim3 = len(self.params3)      
        dist_compare_cat = np.zeros((dim1,dim2,dim3))
        dist_compare_contact = np.zeros((dim1,dim2,dim3))
        
        # Loop through all screening simulations
        for i in range(dim1):
            for j in range(dim2):
                for k in range(dim3):
                    if type == 'free':
                        # freely growing MTs ------------------------------------- 
                        cum_hist = np.array(self.SCAN_cat_times[i][j][k]) # load respective model curve
                        cum_hist = np.append(cum_hist, 0)
                        cum_hist = np.sort(cum_hist)
                        
                        if cum_hist.shape[0] < 5:
                                 dist_compare_cat[i,j,k] = 0        
                        else:
            #                Previously used KS score:
            #                stats_shifted = np.zeros(20+1)
            #                for m in range(0,20+1):
            #                    #stats_shifted[m] = stats.ks_2samp(cum_hist_exp, (1 + FACT.shift_ks_cat - FACT.shift_ks_cat/10*m)*cum_hist)[1]
            #                    stats_shifted[m] = distribution_compare(cum_hist_exp, (1 + FACT.shift_ks_cat - FACT.shift_ks_cat/10*m)*cum_hist)
            #                
            #                dist_compare_cat[i,j,k] = min(stats_shifted)
                             dist_compare_cat[i,j,k] = distribution_compare(cum_hist_exp, cum_hist)

                    elif type == 'barrier':
                        # barrier stalled MTs -------------------------------------  
                        cum_hist = np.array(self.SCAN_barrier_contact_times[i][j][k]) #load respective model curve
                        cum_hist = np.append(cum_hist, 0)
                        cum_hist = np.sort(cum_hist)
                        
                        if cum_hist.shape[0] < 5:
                                dist_compare_contact[i,j,k] = 0        
                        else:
            #                Previously used KS score:                
            #                stats_shifted = np.zeros(20+1)
            #                for m in range(0,20+1):
            #                    #stats_shifted[m] = stats.ks_2samp(cum_hist_exp, (1 + FACT.shift_ks_contact - FACT.shift_ks_contact/10*m)*cum_hist)[1]
            #                    stats_shifted[m] = distribution_compare(cum_hist_exp, (1 + FACT.shift_ks_contact - FACT.shift_ks_contact/10*m)*cum_hist)
            #                
            #                dist_compare_contact[i,j,k] = min(stats_shifted)
                            dist_compare_contact[i,j,k] = distribution_compare(cum_hist_exp, cum_hist)

                    else: 
                        print("Given type must be 'free', 'barrier'.")

        if type == 'free':
            dist_compare_cat[np.where(self.dist_compare_cat == 0)] = np.max(dist_compare_cat)
            self.dist_compare_cat = dist_compare_cat
        elif type == 'barrier':
            dist_compare_contact[np.where(dist_compare_contact == 0)] = np.max(dist_compare_contact)
            self.dist_compare_contact = dist_compare_contact
    
    
    def plot_dist_compare(self, scan_x_label, scan_y_label, type = 'free'):
        """ Plot MSE
        
        Parameters
        ----------
        scan_x_label: str
            Label of x axis
        scan_y_label: str
            Label of y axis
        type: str
            Type of empirical data, freely growing (type='free') or stalled ('barrier'), or 'both'.
        """
        dim1 = len(self.params1)
        dim2 = len(self.params2)
        dim3 = len(self.params3)  
        
        
        if (type == 'barrier') or (type == 'both'):
            cols = min(dim3,2)
            rows = int(np.ceil(dim3/cols))
            
            fig, axs = plt.subplots(nrows=rows, ncols=cols, sharex='col', sharey='row', 
                                    figsize=(6 * cols, 6* dim1/dim2 *rows ))
            
            for i in range(rows):
                for j in range(cols):
                    count = j * rows + i
                    if count < dim3:
                        if cols == 1:
                            ax = axs
                        elif rows >1:
                            ax = axs[j,i]
                        else:
                            ax = axs[count]
                        
                        map = ax.imshow(self.dist_compare_contact[:,:,count]**0.5, aspect='auto', 
                                        vmin = 0, vmax= np.max(self.dist_compare_contact[:,:,count]**0.5),
                                        cmap='gist_heat_r', interpolation='nearest')
                        ax.set_ylabel(scan_y_label)
                        ax.set_xlabel(scan_x_label)
                        fig.colorbar(map, ax = ax, label = 'sqrt(MSE)')
               
                    ax.set_xticks(np.arange(dim1))
                    ax.set_yticks(np.arange(dim2))
                    ax.set_xticklabels(self.params1)
                    ax.set_yticklabels(self.params2)  
            
            fig.suptitle('MSE between cumulative contact time distributions', fontsize=16)  
            
#            if self.sim_parameters["record_data"]:      
    #            filename = foldername + "/" + DATE + '_screen_04'
    #            plt.savefig(filename+'.pdf', format='pdf', dpi=1000) #, transparent=True)
    #            plt.savefig(filename+'.png', format='png', dpi=200) #, transparent=True )      
            plt.show()                  

        if (type == 'free') or (type == 'both'):
            cols = min(dim3,2)
            rows = int(np.ceil(dim3/cols))
            
            fig, axs = plt.subplots(nrows=rows, ncols=cols, sharex='col', sharey='row', 
                                    figsize=(6 * cols, 6* dim1/dim2 *rows ))
            
            for i in range(rows):
                for j in range(cols):
                    count = j * rows + i
                    if count < dim3:
                        if cols == 1:
                            ax = axs
                        elif rows >1:
                            ax = axs[j,i]
                        else:
                            ax = axs[count]
                        
                        map = ax.imshow(self.dist_compare_cat[:,:,count]**0.5, aspect='auto', 
                                        vmin = 0, vmax= np.max(self.dist_compare_cat[:,:,count]**0.5),
                                        cmap='gist_heat_r', interpolation='nearest')

                        ax.set_ylabel(scan_y_label)
                        ax.set_xlabel(scan_x_label)
                        fig.colorbar(map, ax = ax, label = 'sqrt(MSE)')#'P (Kolmogorov-Smirnov)')
               
                    ax.set_xticks(np.arange(dim1))
                    ax.set_yticks(np.arange(dim2))
                    ax.set_xticklabels(self.params1)
                    ax.set_yticklabels(self.params2)  
            
            fig.suptitle('MSE between cumulative catastrophe time distributions', fontsize=16)  
            
#            if self.sim_parameters["record_data"]:        
#                filename = foldername + "/" + DATE + '_screen_05'
#                plt.savefig(filename+'.pdf', format='pdf', dpi=1000) #, transparent=True)
#                plt.savefig(filename+'.png', format='png', dpi=200) #, transparent=True ) 
            
            plt.show() 

        if type == 'both':
            cols = min(dim3,2)
            rows = int(np.ceil(dim3/cols))
            
            fig, axs = plt.subplots(nrows=rows, ncols=cols, sharex='col', sharey='row', 
                                    figsize=(6 * cols, 6* dim1/dim2 *rows ))
            
            for i in range(rows):
                for j in range(cols):
                    count = j * rows + i
                    if count < dim3:
                        if cols == 1:
                            ax = axs
                        elif rows >1:
                            ax = axs[j,i]
                        else:
                            ax = axs[count]
                        
                        map = ax.imshow(np.multiply(self.dist_compare_cat[:,:,count]**0.5,self.dist_compare_contact[:,:,count]**0.5), 
                                        vmin= 0, vmax= np.mean(self.dist_compare_cat[:,:,count]**0.5)*np.mean(self.dist_compare_contact[:,:,count]**0.5), 
                                        aspect='auto', cmap='gist_heat_r', interpolation='nearest')
                        ax.set_ylabel(scan_y_label)
                        ax.set_xlabel(scan_x_label)
                        fig.colorbar(map, ax = ax, label = 'sqrt(MSE*MSE)')#'P (Kolmogorov-Smirnov)')
               
                    ax.set_xticks(np.arange(dim1))
                    ax.set_yticks(np.arange(dim2))
                    ax.set_xticklabels(self.params1)
                    ax.set_yticklabels(self.params2)  
            
            fig.suptitle('Combined MSE between cumulative distributions (=MSE1 * MSE2)', fontsize=16)  
            
#            if self.sim_parameters["record_data"]:   
#                filename = foldername + "/" + DATE + '_screen_06'
#                plt.savefig(filename+'.pdf', format='pdf', dpi=1000) #, transparent=True)
#                plt.savefig(filename+'.png', format='png', dpi=200) #, transparent=True ) 
            
            plt.show() 
    
    
    
            
# ----------------------------------------------------------------------------
# Helper functions
# ----------------------------------------------------------------------------
from scipy.interpolate import interp1d

# Alternative version I tested and played around with
#def distribution_compare(Cum_hist1, Cum_hist2, num_interpol=1000):
#    """ Compare two distributions
#    
#    Parameters
#    ----------
#    Cum_hist1: numpy array
#        (Cumulative) distribution 1
#    Cum_hist2: numpy array
#        (Cumulative) distribution 1
#    num_interpol: int, optional
#        Number of interpolation points. 
#        In general: More is better. Less is faster. Much more than resolution of data will not add anything. 
#    """
#    y1 = 1/(len(Cum_hist1)-1) * np.arange(0, len(Cum_hist1) , 1)
#    y2 = 1/(len(Cum_hist2)-1) * np.arange(0, len(Cum_hist2) , 1)
#    fit1 = interp1d(y1, Cum_hist1, kind='cubic')
#    fit2 = interp1d(y2, Cum_hist2, kind='cubic')
#    ynew = np.linspace(0,0.95, num=num_interpol) # only look at first 95% (avoid overly strong influence of 'weird' tails)
#    return (np.mean((fit1(ynew) - fit2(ynew))**2) / np.mean([np.mean(Cum_hist1), np.mean(Cum_hist2)]))

def distribution_compare(Cum_hist1, Cum_hist2, num_interpol=10000):
    """ Compare two distributions
    
    Parameters
    ----------
    Cum_hist1: numpy array
        (Cumulative) distribution 1
    Cum_hist2: numpy array
        (Cumulative) distribution 1
    num_interpol: int, optional
        Number of interpolation points. 
        In general: More is better. Less is faster. Much more than resolution of data will not add anything. 
    """
    y1 = 1/(len(Cum_hist1)-1) * np.arange(0, len(Cum_hist1) , 1)
    y2 = 1/(len(Cum_hist2)-1) * np.arange(0, len(Cum_hist2) , 1)
    fit1 = interp1d(Cum_hist1, y1, kind='nearest')
    fit2 = interp1d(Cum_hist2, y2, kind='nearest')
    xnew = np.linspace(0,min(max(Cum_hist1),max(Cum_hist2)), num=num_interpol) # only compare the actual overlapping part
    return (np.mean((fit1(xnew) - fit2(xnew))**2))





        

