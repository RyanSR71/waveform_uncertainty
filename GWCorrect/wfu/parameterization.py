import numpy as np
import bilby
import random
import time
import sys
import scipy
import lal
import logging
import tqdm
from bilby.core import utils
from bilby.core.series import CoupledTimeAndFrequencySeries
from bilby.core.utils import PropertyAccessor
from bilby.gw.conversion import convert_to_lal_binary_neutron_star_parameters
from .utils import ProgressBar

def fd_model_difference(hf1,hf2,**kwargs):
    '''
    Generates frequency domain waveform differences between two models hf1 and hf2
    
    Parameters
    ===================
    hf1: bilby.gw.waveform_generator.WaveformGenerator
        frequency domain waveform generator object
    hf2: bilby.gw.waveform_generator.WaveformGenerator
        frequency domain waveform generator object
    injection: dict, optional
        injection parameters if waveform generators do not have parameters in them
        default: None
    npoints: int, optional
        length of the desired frequency grid
        default: 1000
    psd_data: numpy.ndarray, optional
        array containing the psd data and their corresponding frequencies
        default: None
    correction_parameter: float, optional
        fraction of maximum amplitude to cut off the amplitude at
        default: 0.0001
    ref_amplitude: numpy.ndarray, optional
        array of gravitational waveform amplitude
        default None
        
    Returns
    ==================
    frequency_grid: numpy.ndarray
        frequency grid that corresponds to the uncertainty arrays
    amplitude_difference: numpy.ndarray
        array of amplitude difference values
    phase_difference: numpy.ndarray
        array of phase difference values
    '''
    injection = kwargs.get('injection',None)
    npoints = kwargs.get('npoints',1000)
    psd_data = kwargs.get('psd_data',None)
    correction_parameter = kwargs.get('correction_parameter',0.0001)
    ref_amplitude = kwargs.get('ref_amplitude',None)

    f_low = hf1.waveform_arguments['minimum_frequency']
    f_high = hf1.waveform_arguments['maximum_frequency']
    f_ref = hf1.waveform_arguments['reference_frequency']
    
    bilby.core.utils.log.setup_logger(log_level=30)
    np.seterr(all='ignore')

    # adding injection parameters to waveform generators
    if injection is not None:
        hf1.frequency_domain_strain(parameters=injection)
        hf2.frequency_domain_strain(parameters=injection)
    
    # setting up frequency grid and frequency indexes
    start_index = np.argmin(np.abs(hf1.frequency_array - f_ref))+1
    frequency_grid = np.geomspace(f_low,f_high,npoints)
    wf_freqs = np.geomspace(start_index,len(hf1.frequency_array)-1,npoints).astype(int)

    # full waveforms from polarizations
    waveform_1 = hf1.frequency_domain_strain()['plus'][wf_freqs]-hf1.frequency_domain_strain()['cross'][wf_freqs]*1j
    waveform_2 = hf2.frequency_domain_strain()['plus'][wf_freqs]-hf2.frequency_domain_strain()['cross'][wf_freqs]*1j
    
    # waveform amplitudes
    amplitude_1 = np.abs(waveform_1)
    amplitude_2 = np.abs(waveform_2)

    # waveform phases
    phase_1 = np.angle(waveform_1)
    phase_2 = np.angle(waveform_2)
                     
    amplitude_difference = (amplitude_2/amplitude_1) - 1
    unaligned_phase_difference = phase_2-phase_1

    # removing phase shifts of 2pi
    while any(value > 6 for value in [np.abs(unaligned_phase_difference[i+1]-unaligned_phase_difference[i]) for i in range(len(unaligned_phase_difference)-2)]):
        unaligned_phase_difference = np.unwrap(unaligned_phase_difference)

    # finding the position of the cuttoff/discontinuity correction frequency, f_disc
    final_index_1 = list(amplitude_1).index(min(amplitude_1, key=lambda x:np.abs(x-correction_parameter*np.max(amplitude_1))))
    final_index_2 = list(amplitude_2).index(min(amplitude_2, key=lambda x:np.abs(x-correction_parameter*np.max(amplitude_2))))
    final_index = min([final_index_1,final_index_2])
    
    if ref_amplitude is None:
        ref_amplitude = np.copy(amplitude_2)
    ref_amplitude = np.interp(hf1.frequency_array[wf_freqs][0:final_index],f_high*np.linspace(0,1,len(ref_amplitude)),ref_amplitude)
    
    if psd_data is None:
        ref_sigma = np.max(ref_amplitude**2)
    else:
        ref_sigma = np.interp(hf1.frequency_array[wf_freqs][0:final_index], psd_data[:,0],psd_data[:,1])
        
    align_weights = (hf1.frequency_array[wf_freqs][0:final_index]*ref_amplitude**2)/(ref_sigma)
    fit = np.polyfit(hf1.frequency_array[wf_freqs][0:final_index],unaligned_phase_difference[0:final_index],1,w=align_weights)
    unaligned_phase_difference_no_shifts = unaligned_phase_difference[0:final_index]-np.poly1d(fit)(hf1.frequency_array[wf_freqs][0:final_index])
    aligned_phase_difference = np.copy(unaligned_phase_difference)
    aligned_phase_difference[0:final_index] = unaligned_phase_difference_no_shifts
        
    # making the discontinuity correction to amplitude_difference
    amplitude_difference[final_index:] = amplitude_difference[final_index-1]
    phase_difference = np.copy(aligned_phase_difference)
    phase_difference[final_index:] = aligned_phase_difference[final_index-1]
    
    return frequency_grid,amplitude_difference,phase_difference



def parameterization(hf1,hf2,prior,nsamples,**kwargs):
    '''
    Generates samples of waveform uncertainty between two approximants and parameterizes the data with Chebyshev polynomial functions.

    Parameters
    ==================
    hf1: bilby.gw.waveform_generator.WaveformGenerator
        frequency domain waveform generator object
    hf2: bilby.gw.waveform_generator.WaveformGenerator
        frequency domain waveform generator object
    prior: bilby.core.prior.dict.PriorDict
        bilby prior object
    nsamples: int
        number of draws of waveform uncertainty desired
    spline_resolution: int, 500
        number of spline nodes desired
        default: 500
    npoints: int, optional
        length of the desired frequency grid
        default: 1000
    psd_data: numpy.ndarray, optional
        array containing the psd data and their corresponding frequencies
        default: None
    correction_parameter: float, optional
        fraction of maximum amplitude to cut off the amplitude at
        default: 0.0001
    ref_amplitude: numpy.ndarray, optional
        reference amplitude for residual phase calculation
        default: None

    Returns
    ==================
    parameterized_data: numpy.ndarray
        table containing the index, frequency_grid, dA_fit_parameters, dphi_fit_parameters, final_index, dA_final_point, dphi_final_point,
        and injection_parameters for each draw of waveform uncertainty

        frequency_grid: numpy.ndarray
            frequencies corresponding to the frequency parameters specified
        frequency_nodes: numpy.ndarray
            frequency nodes for the splines
        dA_parameters: numpy.ndarray
            amplitude difference spline parameters
        dphi_parameters: numpy.ndarray
            phase difference spline parameters
        injection: dictionary
            source parameters injected into the waveform generators
    '''
    npoints = kwargs.get('npoints',1000)
    psd_data = kwargs.get('psd_data',None)
    correction_parameter = kwargs.get('correction_parameter',0.0001)
    ref_amplitude = kwargs.get('ref_amplitude',None)
    spline_resolution = kwargs.get('spline_resolution',500)

    parameterization_data = np.zeros([nsamples,5],dtype=object)

    # setting the reference amplitude
    if ref_amplitude is None:
        injection = prior.sample()
        ref_waveform = hf2.frequency_domain_strain(parameters=injection)['plus']-hf2.frequency_domain_strain(parameters=injection)['cross']*1j
        ref_amplitude = np.abs(ref_waveform)

    log = logging.getLogger(__name__)
    log.setLevel(logging.INFO)
    log.addHandler(ProgressBar())
    
    for index in tqdm.tqdm(range(nsamples), desc = 'Generating Waveform Differences'):
        
        injection = prior.sample()
        
        # calculating waveform model differences
        frequency_grid,amplitude_difference,phase_difference = fd_model_difference(hf1,hf2,injection=injection,npoints=npoints,psd_data=psd_data,correction_parameter=correction_parameter,ref_amplitude=ref_amplitude)

        spline_indexes = np.linspace(0,len(frequency_grid)-1,spline_resolution).astype(int)
        frequency_nodes = frequency_grid[spline_indexes]
        amplitude_parameters = amplitude_difference[spline_indexes]
        phase_parameters = phase_difference[spline_indexes]
        
        #constructing output matrix
        parameterization_data[index][0] = np.array(frequency_grid)
        parameterization_data[index][1] = np.array(frequency_nodes)
        parameterization_data[index][2] = np.array(amplitude_parameters)
        parameterization_data[index][3] = np.array(phase_parameters)
        parameterization_data[index][4] = injection
    
    return parameterization_data



def recovery_from_parameterization(parameterization_draw, **kwargs):
    '''
    Takes a row from a parameterization matrix and returns the frequency grid, amplitude difference, and phase difference

    Parameters
    ==================
    parameterization_draw: numpy.ndarray
        one row of a parameterization matrix
    dimensionless: bool, optional
        whether or not the output is returned in dimensionless frequency units
        Default: False
    xi_low: float, optional
        if dimensionless is True, this is the lower bound on the dimensionless frequency grid
        default: 0.001
    xi_high: float, optional
        if dimensionless is True, this is the upper bound on the dimensionless frequency grid
        default: 1
    resolution: int, optional
        if dimensionless is True, this is the number of points in the dimensionless frequency grid
        default: 1000

    Returns:
    ==================
    frequency_grid: numpy.ndarray
        array of frequency points
    amplitude_difference: numpy.ndarray
        array of amplitude differences
    phase_difference: numpy.ndarray
        array of phase differences
    '''
    dimensionless = kwargs.get('dimensionless',False)
    xi_low = kwargs.get('xi_low',0.001)
    xi_high = kwargs.get('xi_high',1)
    resolution = kwargs.get('resolution',1000)
    
    frequency_grid = parameterization_draw[0].copy()
    nodes = parameterization_draw[1].copy()
    dA_params = parameterization_draw[2].copy()
    dphi_params = parameterization_draw[3].copy()
    
    if dimensionless is True:
        M = bilby.gw.conversion.generate_mass_parameters(parameterization_draw[4])['total_mass']
        total_frequency_grid = np.geomspace(xi_low,xi_high,resolution)
        frequency_grid *= M/203025.4467280836
        nodes *= M/203025.4467280836
    else:
        total_frequency_grid = frequency_grid.copy()
        
    dA_spline = scipy.interpolate.CubicSpline(nodes,dA_params)
    amplitude_difference = np.interp(total_frequency_grid,frequency_grid,dA_spline(frequency_grid))
    
    dphi_spline = scipy.interpolate.CubicSpline(nodes,dphi_params)
    phase_difference = np.interp(total_frequency_grid,frequency_grid,dphi_spline(frequency_grid))
    
    return total_frequency_grid,amplitude_difference,phase_difference



def uncertainties_from_parameterization(parameterization, **kwargs):
    '''
    Takes in an entire parameterization matrix and returns the Gaussian parameters for the amplitude and phase differences

    Parameters:
    ==================
    parameterization: numpy.ndarry
        parameterization matrix
    dimensionless: bool, optional
        whether or not the output is returned in dimensionless frequency units
        Default: False
    xi_low: float, optional
        if dimensionless is True, this is the lower bound on the dimensionless frequency grid
        default: 0.001
    xi_high: float, optional
        if dimensionless is True, this is the upper bound on the dimensionless frequency grid
        default: 1
    resolution: int, optional
        if dimensionless is True, this is the number of points in the dimensionless frequency grid
        default: 1000

    Returns:
    ==================
    mean_amplitude_difference: numpy.ndarray
        array of the mean amplitude differences as a function of frequency
    amplitude_uncertainty: numpy.ndarray
        array of the standard deviation of the amplitude differences across frequency
    mean_phase_difference: numpy.ndarry
        array of the mean phase differences as a function of frequency
    phase_uncertainty: numpy.ndarray
        array of the stand deviation of the phase differences across frequency
    '''
    dimensionless = kwargs.get('dimensionless',False)
    xi_low = kwargs.get('xi_low',0.001)
    xi_high = kwargs.get('xi_high',1)
    resolution = kwargs.get('resolution',1000)
    
    total_dAs = []
    total_dphis = []
    
    for i in range(len(parameterization)):
        frequency_grid,amplitude_difference,phase_difference = recovery_from_parameterization(parameterization[i],
                                                                                              dimensionless=dimensionless,
                                                                                             xi_low=xi_low,xi_high=xi_high,
                                                                                             resolution=resolution)
        
        total_dAs.append(amplitude_difference)
        total_dphis.append(phase_difference)
    
    mean_amplitude_difference = np.mean(total_dAs,axis=0)
    amplitude_uncertainty = np.std(total_dAs,axis=0)
    
    mean_phase_difference = np.mean(total_dphis,axis=0)
    phase_uncertainty = np.std(total_dphis,axis=0)
    
    return mean_amplitude_difference,amplitude_uncertainty,mean_phase_difference,phase_uncertainty
