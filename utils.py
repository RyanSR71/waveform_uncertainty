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



def injection(data,**kwargs):
    '''
    Pulls a sample out of a parameter dictionary
    
    Parameters
    ==================
    data: dictionary
        dictionary of parameter samples
    index: int, optional
        position within the dict to pull the sample
        default: random integer between zero and the length of the data
        
    Returns
    ==================
    injection: dictionary
        dictionary of injection parameters
    ''' 
    index = kwargs.get('index',random.randint(0,len(data)))
    
    injection = dict()
    for key in data.keys():
        injection[f'{key}']=data[f'{key}'][index]
        
    return injection



def progressBar(count_value, total, suffix=''): #To Do
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
