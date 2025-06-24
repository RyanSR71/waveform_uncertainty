import numpy as np
import bilby
import random
import time as tm
import sys
import scipy
import lal
from bilby.core import utils
from bilby.core.series import CoupledTimeAndFrequencySeries
from bilby.core.utils import PropertyAccessor
from bilby.gw.conversion import convert_to_lal_binary_neutron_star_parameters

def progressBar(count_value, total, suffix=''):
        bar_length = 100
        filled_up_Length = int(round(bar_length* count_value / float(total)))
        percentage = round(100.00 * count_value/float(total),1)
        bar = '=' * filled_up_Length + '-' * (bar_length - filled_up_Length)
        sys.stdout.write('[%s] %s%s %s\r' %(bar, percentage, '%', suffix))
        sys.stdout.flush()

def match(signal,data,PSDs,duration):
    
    signal_match = np.sqrt(bilby.gw.utils.matched_filter_snr(signal,signal,PSDs,4))
    data_match = np.sqrt(bilby.gw.utils.matched_filter_snr(data,data,PSDs,4))
    normalized_match = np.abs(bilby.gw.utils.matched_filter_snr(signal,data,PSDs,4)/(signal_match*data_match))
    
    return normalized_match

def td_waveform(fd_waveform,sampling_frequency):
    reversed_fd_waveform = fd_waveform[::-1]
    total_fd_strain = np.concatenate((fd_waveform,np.conjugate(reversed_fd_waveform[1:-1])))
    td_waveform = sampling_frequency*np.real(np.fft.ifft(total_fd_strain))
    return td_waveform

def smooth_interpolation(full_grid,nodes,parameters,gamma):
    spline = scipy.interpolate.interp1d(nodes,parameters)(nodes)
    temp_grid = np.geomspace(full_grid[1],full_grid[-1],200)
    data = np.interp(temp_grid,nodes,spline)
    new_data = data.copy()
    
    r = int(gamma*len(data))
    if r != 0:
        lower_index = int(gamma*len(data))
        upper_index = int((1-gamma)*len(data))
        for i in range(lower_index,upper_index):
            new_data[i] = (1/(2*r))*np.sum(data[i-r:i+r])
    
    output = np.interp(full_grid,temp_grid,new_data)
    if np.abs(output[0]) > 0:
        output -= output[0]

    return output

def variable_prior(uncertainty,k,xi_low,xi_high):
    frequency_grid = np.linspace(0.001,1,len(uncertainty))
    desired_frequency_nodes = np.geomspace(xi_low,xi_high,k+1)
    
    indexes = [list(frequency_grid).index(min(frequency_grid, key=lambda x:np.abs(x-node))) for node in desired_frequency_nodes]
    frequency_nodes = np.array(frequency_grid[indexes])
    
    coef = np.array([uncertainty[indexes[i]] for i in range(k+1)])

    return frequency_nodes, coef

def maxL(result):
    '''
    Calculates the set of parameters in a posterior that together yield the highest likelihood

    Parameters
    ==================
    result: bilby.core.result.Result
        bilby result object from a parameter estimation run

    Returns
    ==================
    maxL_dict: dictionary
        dictionary of the maximum likelihood values of each of the injected parameters
    '''
    maxL_index = np.argmax(result.log_likelihood_evaluations)
    
    maxL_dict = dict()
    for parameter in result.priors.keys():
        maxL_dict[parameter] = result.posterior[parameter][maxL_index]
        
    return maxL_dict
