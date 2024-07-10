"WaveformUncertainty package"
__version__ = "0.1.0"

import numpy as np
import bilby
import random
import time as tm
import sys
import scipy
from bilby.core import utils
from bilby.core.series import CoupledTimeAndFrequencySeries
from bilby.core.utils import PropertyAccessor
from bilby.gw.conversion import convert_to_lal_binary_neutron_star_parameters



def fd_model_difference(hf1,hf2,**kwargs):
    '''
    Generates frequency domain waveform differences between two models hf1 and hf2
    
    Parameters
    ===================
    hf1 & hf2: bilby.gw.waveform_generator.WaveformGenerator
        frequency domain waveform generator object WITH injected parameters (for strain calculation)
    f_low: float, optional
        minimum frequency
        default: 20.0
    f_high: float, optional
        maximum frequency
        default: 2048.0
    f_ref: float, optional
        reference frequency
        default: 50.0
    npoints: int, optional
        length of the desired frequency grid
        default: 1000
    polarization: string, optional
        polarization of the strain data (plus or cross)
        default: 'plus'
    psd_data: numpy.ndarray, optional
        array containing the psd data and their corresponding frequencies
        default: None
    correction_parameter: float, optional
        value at which to cut the second derivative of amplitude difference (see WFU_equations #1)
        default: -10e-6 FIX THIS
        
    Returns
    ==================
    frequency_grid: numpy.ndarray
        frequency grid that corresponds to the uncertainty arrays
    amplitude_difference: numpy.ndarray
        array of amplitude difference values
    phase_difference: numpy.ndarray
        array of phase difference values; if psd data is None, raw_phase_difference will be returned, residual_phase_difference otherwise
    amplitude_difference_final_point: float
        value of the amplitude difference where the discontinuity correction starts
    phase_difference_final_point: float
        value of the phase difference where the discontinuity correction starts
    final_index: int
        position within the frequency grid where the discontinuity correction starts
    '''
    f_low = kwargs.get('f_low',20.0)
    f_high = kwargs.get('f_high',2048.0)
    f_ref = kwargs.get('f_ref',50.0)
    npoints = kwargs.get('npoints',1000)
    polarization = kwargs.get('polarization','plus')
    psd_data = kwargs.get('psd_data',None)
    correction_parameter = kwargs.get('correction_parameter',-10e-6)
    
    bilby.core.utils.log.setup_logger(log_level=30)
    np.seterr(all='ignore')
    
    start_index = np.argmin(np.abs(hf1.frequency_array - f_ref))+1
    frequency_grid = np.geomspace(f_low,f_high,npoints)
    wf_freqs = np.geomspace(start_index,len(hf1.frequency_array)-1,npoints).astype(int)
    
    amplitude_1 = np.abs(hf1.frequency_domain_strain()[f'{polarization}'][wf_freqs])
    amplitude_2 = np.abs(hf2.frequency_domain_strain()[f'{polarization}'][wf_freqs])
    
    phase_1 = np.unwrap(np.unwrap(np.unwrap(np.unwrap(np.unwrap(np.angle(hf1.frequency_domain_strain()[f'{polarization}'][wf_freqs]))))))
    phase_2 = np.unwrap(np.unwrap(np.unwrap(np.unwrap(np.unwrap(np.angle(hf2.frequency_domain_strain()[f'{polarization}'][wf_freqs]))))))
                     
    amplitude_difference = (amplitude_1-amplitude_2)/amplitude_1
    raw_phase_difference = phase_1-phase_2
    
    if psd_data is not None:
        ref_amplitude = np.abs(hf1.frequency_domain_strain()['plus'][wf_freqs])
        ref_sigma = np.interp(hf1.frequency_array[wf_freqs], psd_data[:,0],psd_data[:,1])
        align_weights = ref_amplitude*ref_amplitude / ref_sigma * hf1.frequency_array[wf_freqs]
        fit = np.polyfit(hf1.frequency_array[wf_freqs],raw_phase_difference,1,w=align_weights)
        residual_phase_difference = raw_phase_difference - np.poly1d(fit)(hf1.frequency_array[wf_freqs])
        
    def derivative(x,y):
        output = []
        for i in range(0,len(x)-1):
            slope = (y[i+1]-y[i])/(x[i+1]-x[i])
            output.append(slope)
        new_x = list(np.copy(x))
        new_x.remove(x[-1])
        return new_x, output  
    f_prime, amplitude_difference_first_derivative = derivative(hf1.frequency_array[wf_freqs],amplitude_difference)
    f_double_prime,amplitude_difference_second_derivative = derivative(f_prime,amplitude_difference_first_derivative)
    for i in range(len(amplitude_difference_second_derivative)):
        if amplitude_difference_second_derivative[i] < correction_parameter:
            final_index = i
            break
    else: 
        final_index = len(hf1.frequency_array[wf_freqs])-1
        
    amplitude_difference[final_index:] = amplitude_difference[final_index-1]
    
    amplitude_difference_final_point = amplitude_difference[final_index-1]
    
    
    if psd_data is not None:
        residual_phase_difference[final_index:] = residual_phase_difference[final_index-1]
        residual_phase_difference_final_point = residual_phase_difference[final_index-1]
        return frequency_grid,amplitude_difference,residual_phase_difference,amplitude_difference_final_point,residual_phase_difference_final_point,final_index
    else:
        raw_phase_difference[final_index:] = raw_phase_difference[final_index-1]
        raw_phase_difference_final_point = raw_phase_difference[final_index-1]
        return frequency_grid,amplitude_difference,raw_phase_difference,amplitude_difference_final_point,raw_phase_difference_final_point,final_index

    
    
def parameter_dict_from_prior(prior,nsamples):
    '''
    Takes a bilby prior object and constructs a dictionary of random samples
    
    Parameters
    ==================
    prior: bilby.core.prior.dict.PriorDict
        prior object
    nsamples: int
        number of samples desired
        
    Returns
    ==================
    parameter_dict: dictionary
        dictionary of nsamples random samples
    '''
    parameter_dict = dict()
    for key in prior.keys():
        parameter = np.zeros(nsamples)
        for i in range(nsamples):
            parameter[i] = prior[f'{key}'].sample()
        parameter_dict[f'{key}'] = parameter
    
    return parameter_dict



def injection(data,**kwargs):
    '''
    Pulls a sample out of a parameter dictionary
    
    Parameters
    ==================
    data: dictionary
        Dictionary of parameter samples (from paramater_dict_from_prior())
    index: int, optional
        position within the dict to pull the sample
        default: random integer between zero and the length of the data
    precession: bool, optional
        only True if waveform approximant supports precessing spins 
        if False, precessing spin parameters will be removed/replaced with non-precessing parameters
        default: False
    tides: bool, optional
        only True if waveform approximant supports tidal deformabilities
        if False, tidal parameters will be removed
        default: True
        
    Returns
    ==================
    injection: dictionary
        dictionary of injection parameters
    ''' 
    index = kwargs.get('index',random.randint(0,len(data)))
    precession = kwargs.get('precession',False)
    tides = kwargs.get('tides',True)
    
    injection = dict()
    for key in data.keys():
        injection[f'{key}']=data[f'{key}'][index]
    
    if tides==False:
        parameters = ['lambda_1','lambda_2','lambda_tilde','delta_lambda_tilde']
        for parameter in parameters:
            if parameter in data.keys():
                del injection[f'{parameter}']
                
    if 'a_1' in data.keys() and precession==False:
        parameters = ['chi_1','chi_2','chi_1_in_plane','chi_2_in_plane','cos_tilt_1','cos_tilt_2']
        for parameter in parameters:
            if parameter in data.keys():
                del injection[f'{parameter}']
        injection['chi_1'] = injection['a_1']*np.cos(injection['tilt_1'])
        injection['chi_2'] = injection['a_2']*np.cos(injection['tilt_2'])
        del injection['a_1']
        del injection['a_2']
        del injection['tilt_1']
        del injection['tilt_2']
        
    return injection



def parameterization(approximant1,approximant2,parameter_data,nsamples,**kwargs):
    '''
    Generates samples of waveform uncertainty between two approximants and parameterizes the data with Chebyshev polynomial functions.

    Parameters
    ==================
    approximant1: string
        name of the first waveform approximant
    approximant2: string
        name of the second waveform approximant
    parameter_data: dictionary or bilby.core.prior.dict.PriorDict
        dictionary containing neutron star parameter samples or a bilby prior object that will be converted into a dictionary
    nsamples: int
        number of draws of waveform uncertainty desired
    precession: bool, optional
        True if both waveform approximants support precessing spins 
        if False, precessing spin parameters will be removed/replaced with non-precessing parameters
        default: False
    tides: bool, optional
        only True if both waveform approximants support tidal deformabilities
        if False, tidal parameters will be removed
        default: True
    fit_parameters: int, optional
        number of terms to use in the parameterization
        default: 15
    npoints: int, optional
        length of the desired frequency grid
        default: 1000
    f_low: float, optional
        lower bound on the frequency grid
        default: 20.0
    f_high: float, optional
        upper bound on the frequency grid
        default: 2048.0
    f_ref: float, optional
        reference frequency
        default: 50.0
    sampling_frequency: float, optional
        sampling frequency or sampling rate
        default: 4096
    duration: float, optional
        duration of the signal
        default: 256
    max_dA_error: float [%], optional
        maximum allowed error between the amplitude uncertainty and its parameterization
        default: 2 
    max_dphi_error: float [degrees], optional
        maximum allowed error between the phase uncertainty and its parameterization
        default: 2
    psd_data: numpy.ndarray, optional
        array containing the psd data and their corresponding frequencies
        default: None
    correction_parameter: float, optional
        value at which to cut the second derivative of amplitude difference (see WFU_equations #1)
        default: -10e-6
    polarization: string, optional
        polarization of the strain data (plus or cross)
        default: 'plus'
    fit_threshold: float [%], optional
        minimum parameterization success rate
        default: 75

    Returns
    ==================
    parameterized_data: numpy.ndarray
        table containing the index, frequency_grid, dA_fit_parameters, dphi_fit_parameters, final_index, dA_final_point, dphi_final_point,
        and injection_parameters for each draw of waveform uncertainty

        index: int
            position within the list of indexes the waveform uncertainty draws were drawn from (only for debugging purposes)
        frequency_grid: numpy.ndarray
            frequencies corresponding to the frequency parameters specified
        dA_fit_parameters: numpy.ndarray
            parameters corresponding to the parameterization of the amplitude uncertainty
        dphi_fit_parameters: numpy.ndarray
            parameters corresponding to the parameterization of the phase uncertainty
        final_index: int
            position within the frequency grid where the tidal correction starts
        dA_final_point: float
            value of the amplitude uncertainty that the tidal correction levels off at
        dphi_final_point: float
            value of the phase uncertainty that the tidal correction levels off at
        injection parameters: dictionary
            neutron star parameters injected into the waveform generators

    '''
    f_low = kwargs.get('f_low',20.0)
    f_high = kwargs.get('f_high',2048.0)
    f_ref = kwargs.get('f_ref',50.0)
    npoints = kwargs.get('npoints',1000)
    polarization = kwargs.get('polarization','plus')
    psd_data = kwargs.get('psd_data',None)
    correction_parameter = kwargs.get('correction_parameter',-10e-6)
    precession = kwargs.get('precession',False)
    tides = kwargs.get('tides',True)
    fit_parameters = kwargs.get('fit_parameters',15)
    sampling_frequency = kwargs.get('sampling_frequency',4096)
    duration = kwargs.get('duration',256)
    max_amplitude_error = kwargs.get('max_amplitude_error',2)
    max_phase_error = kwargs.get('max_phase_error',2)
    fit_threshold = kwargs.get('fit_threshold',75)
        
    np.seterr(all='ignore')
    progress = 1
    bilby.core.utils.log.setup_logger(log_level=30)
    
    if type(parameter_data)==bilby.core.prior.dict.PriorDict:
        data = parameter_dict_from_prior(parameter_data,2*nsamples)
    else:
        data = parameter_data
    
    index_samples=list(range(len(data["a_1"])))
    indexes=[]
    for draws in range(len(index_samples)):
        index=random.choice(index_samples)
        indexes.append(index)
        index_samples.remove(index)
    
    wfargs1 = dict(waveform_approximant=approximant1, reference_frequency=f_ref, catch_waveform_errors=True)
    wfargs2 = dict(waveform_approximant=approximant2, reference_frequency=f_ref, catch_waveform_errors=True)
    
    def recovery(identity,data):
        if str(identity) == 'amplitude_difference':
            parameterized_curve = np.polynomial.chebyshev.chebval(data[1][0:data[4]],data[2])
            post_waveform_uncertainties = np.copy(data[1])
            post_waveform_uncertainties[data[4]-2:] = data[5]
        elif str(identity) == 'phase_difference':
            parameterized_curve = np.polynomial.chebyshev.chebval(data[1][0:data[4]],data[3])
            post_waveform_uncertainties = np.copy(data[1])
            post_waveform_uncertainties[data[4]-2:] = data[6]
        return np.concatenate((parameterized_curve[:-2],post_waveform_uncertainties[data[4]-2:]))
    
    def progressBar(count_value, total, suffix=''):
        bar_length = 100
        filled_up_Length = int(round(bar_length* count_value / float(total)))
        percentage = round(100.00 * count_value/float(total),1)
        bar = '=' * filled_up_Length + '-' * (bar_length - filled_up_Length)
        sys.stdout.write('[%s] %s%s %s\r' %(bar, percentage, '%', suffix))
        sys.stdout.flush()
    
    trial = 0
    final_indexes = []
    extraneous_indexes = []
    
    output_matrix = np.zeros([len(indexes),8],dtype=object)
    start = tm.time()
    
    print("Parameterizing...")
    
    for index in indexes:
    
        progressBar(progress,(nsamples))

        hf1 = bilby.gw.WaveformGenerator(
            parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_neutron_star_parameters,
            parameters=injection(data,index=index,precession=precession,tides=tides), 
            waveform_arguments=wfargs1,
            frequency_domain_source_model=bilby.gw.source.lal_binary_neutron_star, 
            sampling_frequency=sampling_frequency, 
            duration=duration
        )
        hf2 = bilby.gw.WaveformGenerator(
            parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_neutron_star_parameters,
            parameters=injection(data,index=index,precession=precession,tides=tides), 
            waveform_arguments=wfargs2,
            frequency_domain_source_model=bilby.gw.source.lal_binary_neutron_star, 
            sampling_frequency=sampling_frequency, 
            duration=duration
        )
        
        frequency_grid,amplitude_difference,phase_difference,amplitude_difference_final_point,phase_difference_final_point,final_index = fd_model_difference(hf1,hf2,f_low=f_low,f_high=f_high,f_ref=f_ref,npoints=npoints,polarization=polarization,psd_data=psd_data,correction_parameter=correction_parameter)
        
        amplitude_difference_fit = np.polynomial.chebyshev.Chebyshev.fit((frequency_grid[0:final_index]),amplitude_difference[0:final_index],fit_parameters-1)
        amplitude_difference_parameters = amplitude_difference_fit.convert().coef

        phase_difference_fit = np.polynomial.chebyshev.Chebyshev.fit((frequency_grid[0:final_index]),phase_difference[0:final_index],fit_parameters-1)
        phase_difference_parameters = phase_difference_fit.convert().coef

        output_matrix[index][0] = index
        output_matrix[index][1] = np.array(frequency_grid)
        output_matrix[index][2] = np.array(amplitude_difference_parameters)
        output_matrix[index][3] = np.array(phase_difference_parameters)
        output_matrix[index][4] = final_index
        output_matrix[index][5] = amplitude_difference_final_point
        output_matrix[index][6] = phase_difference_final_point
        output_matrix[index][7] = injection(data,index=index,precession=precession,tides=tides)

        final_indexes.append(index)
        
        relative_amplitude_error = np.max(np.abs(100*((amplitude_difference-recovery('amplitude_difference',output_matrix[index])))/([max_amplitude_error]*len(amplitude_difference))))
        relative_phase_error = np.max(np.abs((100*((phase_difference-recovery('phase_difference',output_matrix[index])))/([max_phase_error*(2*np.pi/360)]*len(phase_difference)))))

        if relative_amplitude_error>100 or relative_phase_error>100:
            final_indexes.remove(index)
            extraneous_indexes.append(index)
        else:
            progress += 1

        if 100*len(final_indexes)/(len(final_indexes)+len(extraneous_indexes)) < fit_threshold:
            trial += 1
        if trial == 20:
            parameterization_rate = np.round(100*len(final_indexes)/(len(final_indexes)+len(extraneous_indexes)),2)
            raise Exception(f"Parameterization Rate Too Low; It must be above {fit_threshold}%, but was at {parameterization_rate}%.")
        if len(final_indexes) == nsamples:
            break
    
    empty_data_matrix = np.zeros([len(final_indexes),8],dtype=object)
    count = 0
    for index in final_indexes:
        for k in range(8):
            empty_data_matrix[count][k] = output_matrix[index][k]
        count += 1
    parameterized_data = np.copy(empty_data_matrix)
    
    print("")
    print("Done!")
    print("")
    print(f"Time Elapsed: {round(tm.time()-start,3)} seconds")
    print("")
    final_parameterization_rate = round(100*(nsamples/(len(final_indexes)+len(extraneous_indexes))),4)
    print(f"Parameterized Model Difference data was created with {len(final_indexes)} sample sets at a parameterization rate of {final_parameterization_rate}%.") 
    print("")
    
    return parameterized_data



def recovery_from_parameterization(identity,data):
    '''
    Converts a parameterized set of waveform difference (output of WaveformUncertainty.parameterization()) back into waveform difference arrays
    
    Parameters
    ==================
    identity: string
        specifies which waveform difference should be returned (amplitude_difference or phase_difference)
    data: numpy.ndarray
        one index of the output matrix from WaveformUncertainty.parameterization(); input WaveformUncertainty.parameterization()[index]
    
    Returns
    ==================
    difference_array: numpy.ndarray
        array of the waveform difference converted from the parameterization; has the same shape as the frequency grid within the original matrix
    '''
    if str(identity) == 'amplitude_difference':
        parameterized_curve = np.polynomial.chebyshev.chebval(data[1][0:data[4]],data[2])
        post_waveform_uncertainties = np.copy(data[1])
        post_waveform_uncertainties[data[4]-2:] = data[5]
    
    elif str(identity) == 'phase_difference':
        parameterized_curve = np.polynomial.chebyshev.chebval(data[1][0:data[4]],data[3])
        post_waveform_uncertainties = np.copy(data[1])
        post_waveform_uncertainties[data[4]-2:] = data[6]
    
    else:
        raise Exception('Identity of the uncertainty must be "amplitude_difference" or "phase_difference".')
    
    return np.concatenate((parameterized_curve[:-2],post_waveform_uncertainties[data[4]-2:]))



def uncertainties_from_parameterization(data,**kwargs):
    '''
    Takes all of the sets in a parameterized waveform difference matrix and takes the mean and standard deviation of amplitude and phase difference
    
    Parameters
    ==================
    data: numpy.ndarray
        WaveformUncertainty.parameterization() output matrix
    linear: bool, optional
        if True, default geometric frequency grid will be replaced with a linear one; useful for waveform uncertainty sampling
        default: False
    resolution: float, optional
        distance between points in the linear frequency grid
        default: None
        
    Returns
    ==================
    mean_amplitude_difference: numpy.ndarray
        array of the mean value of the amplitude difference corresponding to the frequency grid
    amplitude_uncertainty: numpy.ndarray
        standard deviations of the amplitude difference
    mean_phase_difference: numpy.ndarray
        array of the mean value of the phase difference corresponding to the frequency grid
    phase_uncertainty: numpy.ndarray
        standard deviations of the phase difference
    linear_frequency_grid: numpy.ndarray
        only if linear=True, new linear frequency grid
    '''
    linear = kwargs.get('linear',False)
    resolution = kwargs.get('resolution',None)
    
    amplitude_difference_data = []
    for index in range(len(data)):
        amplitude_difference_data.append(recovery_from_parameterization('amplitude_difference',data[index]))

    phase_difference_data = []
    for index in range(len(data)):
        phase_difference_data.append(recovery_from_parameterization('phase_difference',data[index]))

    if linear==True and resolution is not None:
    
        linear_frequency_grid = np.arange(data[0][1][0],data[0][1][-1]+resolution,resolution)
        amplitude_difference = []
        for i in range(len(amplitude_difference_data)):
            amplitude_difference.append(np.interp(linear_frequency_grid,data[0][1],amplitude_difference_data[i]))

        phase_difference = []
        for i in range(len(phase_difference_data)):
            phase_difference.append(np.interp(linear_frequency_grid,data[0][1],phase_difference_data[i]))
        
        mean_amplitude_difference = np.mean(amplitude_difference,axis=0)
        amplitude_uncertainty = np.std(amplitude_difference,axis=0)

        mean_phase_difference = np.mean(phase_difference,axis=0)
        phase_uncertainty = np.std(phase_difference,axis=0)
        
        return mean_amplitude_difference,amplitude_uncertainty,mean_phase_difference,phase_uncertainty,linear_frequency_grid   
    
    else:
        mean_amplitude_difference = np.mean(amplitude_difference_data,axis=0)
        amplitude_uncertainty = np.std(amplitude_difference_data,axis=0)

        mean_phase_difference = np.mean(phase_difference_data,axis=0)
        phase_uncertainty = np.std(phase_difference_data,axis=0)
        
        return mean_amplitude_difference,amplitude_uncertainty,mean_phase_difference,phase_uncertainty
        


def WFU_prior(mean_amplitude_difference,amplitude_uncertainty,mean_phase_difference,phase_uncertainty,frequency_grid,nnodes,**kwargs):
    '''
    Automatically generates a bilby prior object containing Gaussian waveform uncertainty parameter priors (alphas and betas)
    If given a pre-existing prior object, the waveform uncertainty parameters will be added to it

    Parameters
    =================
    mean_amplitude_difference: numpy.ndarray
        array of mean amplitude difference values
    amplitude_uncertainty: numpy.ndarray
        array of the amplitude uncertainty; defined as the standard deviation of amplitude differences
    mean_phase_difference: numpy.ndarray
        array of mean phase difference values
    phase_uncertainty: numpy.ndarray
        array of the phase uncertainty; defined as the standard deviation of phase differences
    frequency_grid: numpy.ndarray
        frequency grid corresponding to the waveform uncertainties
    nnodes: int
        number of frequency nodes desired
    prior: bilby.core.prior.PriorDict, optional
        if given, the output prior will simply be added to this
        default: None
    spacing: string, optional
        dictates the type of progression the frequency nodes will have; 'linear' or 'geometric'
        default: 'linear'

    Returns
    ==================
    prior: bilby.core.prior.PriorDict
        prior containing the waveform uncertainty parameters (alphas and betas)
    frequency_nodes: numpy.ndarray
        frequency nodes used by __WaveformGeneratorWFU() to generate waveform difference splines
    '''
    prior = kwargs.get('prior',None)
    spacing = kwargs.get('spacing','linear')
    
    if prior is None:
        prior = bilby.core.prior.PriorDict()
    
    if spacing == 'linear':
        frequency_scale = np.linspace(0,len(frequency_grid)-1,nnodes).astype(int)
    elif spacing == 'geometric':
        start_index = int(len(frequency_grid)/16)
        frequency_scale = [0]
        for i in range(nnodes-1):
            frequency_scale.append(np.geomspace(start_index,len(frequency_grid)-1,nnodes-1).astype(int)[i])
    
    frequency_nodes = frequency_grid[frequency_scale]
    
    for i in range(len(frequency_scale)):
        prior[f'alpha_{i+1}'] = bilby.core.prior.Gaussian(name=f'alpha_{i+1}',latex_label=r'$\alpha_{n}$'.replace('n',str(i+1)),
                                                        mu=mean_amplitude_difference[frequency_scale[i]],
                                                          sigma=amplitude_uncertainty[frequency_scale[i]])
    for i in range(len(frequency_scale)):
        prior[f'beta_{i+1}'] = bilby.core.prior.Gaussian(name=f'beta_{i+1}',latex_label=r'$\beta_{n}$'.replace('n',str(i+1)),
                                                        mu=mean_phase_difference[frequency_scale[i]],
                                                         sigma=phase_uncertainty[frequency_scale[i]])
    
    return prior,frequency_nodes



class WaveformGeneratorWFU(object):
    '''
    Modified WaveformGenerator object from bilby.gw to include waveform uncertainty corrections in the strain calculation
    To sample waveform uncertainty, include all relevant "alpha" and "beta" parameters in the prior.
    Note: make sure the number of alphas, betas, and waveform_uncertainty_nodes are the same
    
    New Parameters
    ==================
    waveform_uncertainty_nodes: numpy.ndarray, optional
        array of frequency nodes to be used in generating the dA and dphi splines
        default: None
    dA_sampling: bool, optional
        if True, the waveform generator will attempt to pull alpha parameters from the parameter dictionary (either an injection or the prior)
        default: None
    dphi_sampling: bool, optional
        if True, the waveform generator will attempt to pull beta parameters from the parameter dictionary (either an injection or the prior)
        default: None
    '''
    duration = PropertyAccessor('_times_and_frequencies', 'duration')
    sampling_frequency = PropertyAccessor('_times_and_frequencies', 'sampling_frequency')
    start_time = PropertyAccessor('_times_and_frequencies', 'start_time')
    frequency_array = PropertyAccessor('_times_and_frequencies', 'frequency_array')
    time_array = PropertyAccessor('_times_and_frequencies', 'time_array')
    def __init__(self, duration=None, sampling_frequency=None, start_time=0, frequency_domain_source_model=None,
                 time_domain_source_model=None, parameters=None,
                 waveform_uncertainty_nodes=None,dA_sampling=False,dphi_sampling=False,
                 parameter_conversion=None,
                 waveform_arguments=None):
        self._times_and_frequencies = CoupledTimeAndFrequencySeries(duration=duration,
                                                                    sampling_frequency=sampling_frequency,
                                                                    start_time=start_time)
        self.frequency_domain_source_model = frequency_domain_source_model
        self.time_domain_source_model = time_domain_source_model
        self.source_parameter_keys = self.__parameters_from_source_model()
        self.dA_sampling = dA_sampling
        self.dphi_sampling=dphi_sampling
        
        if parameter_conversion is None:
            self.parameter_conversion = convert_to_lal_binary_black_hole_parameters
        else:
            self.parameter_conversion = parameter_conversion
        if waveform_arguments is not None:
            self.waveform_arguments = waveform_arguments
        else:
            self.waveform_arguments = dict()
        if waveform_uncertainty_nodes is not None:
            self.waveform_uncertainty_nodes = waveform_uncertainty_nodes
        else:
            self.waveform_uncertainty_nodes = None
        if isinstance(parameters, dict):
            self.parameters = parameters
        self._cache = dict(parameters=None, waveform=None, model=None)
        
        utils.logger.info(
            "Waveform generator initiated with\n"
            "  frequency_domain_source_model: {}\n"
            "  time_domain_source_model: {}\n"
            "  parameter_conversion: {}"
            .format(utils.get_function_path(self.frequency_domain_source_model),
                    utils.get_function_path(self.time_domain_source_model),
                    utils.get_function_path(self.parameter_conversion))
        )
        
    def __repr__(self):
        if self.frequency_domain_source_model is not None:
            fdsm_name = self.frequency_domain_source_model.__name__
        else:
            fdsm_name = None
        if self.time_domain_source_model is not None:
            tdsm_name = self.time_domain_source_model.__name__
        else:
            tdsm_name = None
        if self.parameter_conversion is None:
            param_conv_name = None
        else:
            param_conv_name = self.parameter_conversion.__name__
        return self.__class__.__name__ + '(duration={}, sampling_frequency={}, start_time={}, ' \
                                         'frequency_domain_source_model={}, time_domain_source_model={}, ' \
                                         'parameter_conversion={}, ' \
                                         'waveform_uncertainty_nodes={}, ' \
                                         'dA_sampling={}, ' \
                                         'dphi_sampling={}, ' \
                                         'waveform_arguments={})'\
            .format(self.duration, self.sampling_frequency, self.start_time, fdsm_name, tdsm_name,
                    param_conv_name,self.waveform_uncertainty_nodes, self.waveform_arguments)
    
    def frequency_domain_strain(self, parameters=None):
        return self._calculate_strain(model=self.frequency_domain_source_model,
                                      model_data_points=self.frequency_array,
                                      parameters=parameters,
                                      transformation_function=utils.nfft,
                                      transformed_model=self.time_domain_source_model,
                                      transformed_model_data_points=self.time_array,
                                      waveform_uncertainty_nodes=self.waveform_uncertainty_nodes)
    
    def time_domain_strain(self, parameters=None):
        return self._calculate_strain(model=self.time_domain_source_model,
                                      model_data_points=self.time_array,
                                      parameters=parameters,
                                      transformation_function=utils.infft,
                                      transformed_model=self.frequency_domain_source_model,
                                      transformed_model_data_points=self.frequency_array)

    def _calculate_strain(self, model, model_data_points, transformation_function, transformed_model,
                          transformed_model_data_points, parameters, waveform_uncertainty_nodes):
        if parameters is not None:
            self.parameters = parameters
        if self.parameters == self._cache['parameters'] and self._cache['model'] == model and \
                self._cache['transformed_model'] == transformed_model:
            return self._cache['waveform']
        if model is not None:
            model_strain = self._strain_from_model(model_data_points, model)
        elif transformed_model is not None:
            model_strain = self._strain_from_transformed_model(transformed_model_data_points, transformed_model,
                                                               transformation_function)
        else:
            raise RuntimeError("No source model given")
        self._cache['waveform'] = model_strain
        self._cache['parameters'] = self.parameters.copy()
        self._cache['model'] = model
        self._cache['transformed_model'] = transformed_model
        
        '''
        The following block performs the waveform uncertainty correction:
        '''
        
        if self.waveform_uncertainty_nodes is not None:
            alphas = np.zeros(len(self.waveform_uncertainty_nodes))
            betas = np.zeros(len(self.waveform_uncertainty_nodes))
            for i in range(len(self.waveform_uncertainty_nodes)):
                if self.dA_sampling == True:
                    alphas[i]=parameters['alpha_' + f'{i+1}']
                if self.dphi_sampling == True:
                    betas[i]=parameters['beta_' + f'{i+1}']
            dA = scipy.interpolate.CubicSpline(self.waveform_uncertainty_nodes,alphas)(self.frequency_array)
            dphi = scipy.interpolate.CubicSpline(self.waveform_uncertainty_nodes,betas)(self.frequency_array)
            model_strain['plus'] *= (1+dA)*(2+dphi*1j)/(2-dphi*1j)
            model_strain['cross'] *= (1+dA)*(2+dphi*1j)/(2-dphi*1j)
            
        return model_strain
    
    def _strain_from_model(self, model_data_points, model):
        return model(model_data_points, **self.parameters)
    
    def _strain_from_transformed_model(self, transformed_model_data_points, transformed_model, transformation_function):
        transformed_model_strain = self._strain_from_model(transformed_model_data_points, transformed_model)
        if isinstance(transformed_model_strain, np.ndarray):
            return transformation_function(transformed_model_strain, self.sampling_frequency)
        model_strain = dict()
        for key in transformed_model_strain:
            if transformation_function == utils.nfft:
                model_strain[key], _ = \
                    transformation_function(transformed_model_strain[key], self.sampling_frequency)
            else:
                model_strain[key] = transformation_function(transformed_model_strain[key], self.sampling_frequency)
        return model_strain
    
    @property
    def parameters(self):
        return self.__parameters
    
    @parameters.setter 
    def parameters(self, parameters):
        if not isinstance(parameters, dict):
            raise TypeError('"parameters" must be a dictionary.')
        new_parameters = parameters.copy()
        new_parameters, _ = self.parameter_conversion(new_parameters)
        self.__parameters = new_parameters
        self.__parameters.update(self.waveform_arguments)
        
    def __parameters_from_source_model(self):
        if self.frequency_domain_source_model is not None:
            model = self.frequency_domain_source_model
        elif self.time_domain_source_model is not None:
            model = self.time_domain_source_model
        else:
            raise AttributeError('Either time or frequency domain source '
                                 'model must be provided.')
        return set(utils.infer_parameters_from_function(model))


class WaveformGeneratorWFU_inj(object):
    
    '''
    Temporary waveform generator; takes in arrays of amplitude and phase difference to directly modify the strain data
    '''
    
    duration = PropertyAccessor('_times_and_frequencies', 'duration')
    sampling_frequency = PropertyAccessor('_times_and_frequencies', 'sampling_frequency')
    start_time = PropertyAccessor('_times_and_frequencies', 'start_time')
    frequency_array = PropertyAccessor('_times_and_frequencies', 'frequency_array')
    time_array = PropertyAccessor('_times_and_frequencies', 'time_array')
    def __init__(self, duration=None, sampling_frequency=None, start_time=0, frequency_domain_source_model=None,
                 time_domain_source_model=None, parameters=None,
                 amplitude_difference=None,phase_difference=None,
                 parameter_conversion=None,
                 waveform_arguments=None):
        self._times_and_frequencies = CoupledTimeAndFrequencySeries(duration=duration,
                                                                    sampling_frequency=sampling_frequency,
                                                                    start_time=start_time)
        self.frequency_domain_source_model = frequency_domain_source_model
        self.time_domain_source_model = time_domain_source_model
        self.source_parameter_keys = self.__parameters_from_source_model()
        self.amplitude_difference = amplitude_difference
        self.phase_difference = phase_difference
        
        if parameter_conversion is None:
            self.parameter_conversion = convert_to_lal_binary_black_hole_parameters
        else:
            self.parameter_conversion = parameter_conversion
        if waveform_arguments is not None:
            self.waveform_arguments = waveform_arguments
        else:
            self.waveform_arguments = dict()
        if isinstance(parameters, dict):
            self.parameters = parameters
        self._cache = dict(parameters=None, waveform=None, model=None)
        
        utils.logger.info(
            "Waveform generator initiated with\n"
            "  frequency_domain_source_model: {}\n"
            "  time_domain_source_model: {}\n"
            "  parameter_conversion: {}"
            .format(utils.get_function_path(self.frequency_domain_source_model),
                    utils.get_function_path(self.time_domain_source_model),
                    utils.get_function_path(self.parameter_conversion))
        )
        
    def __repr__(self):
        if self.frequency_domain_source_model is not None:
            fdsm_name = self.frequency_domain_source_model.__name__
        else:
            fdsm_name = None
        if self.time_domain_source_model is not None:
            tdsm_name = self.time_domain_source_model.__name__
        else:
            tdsm_name = None
        if self.parameter_conversion is None:
            param_conv_name = None
        else:
            param_conv_name = self.parameter_conversion.__name__
        return self.__class__.__name__ + '(duration={}, sampling_frequency={}, start_time={}, ' \
                                         'frequency_domain_source_model={}, time_domain_source_model={}, ' \
                                         'parameter_conversion={}, ' \
                                         'amplitude_difference={}, ' \
                                         'phase_difference={}, ' \
                                         'waveform_arguments={})'\
            .format(self.duration, self.sampling_frequency, self.start_time, fdsm_name, tdsm_name,
                    param_conv_name,self.waveform_uncertainty_nodes, self.waveform_arguments)
    
    def frequency_domain_strain(self, parameters=None):
        return self._calculate_strain(model=self.frequency_domain_source_model,
                                      model_data_points=self.frequency_array,
                                      parameters=parameters,
                                      transformation_function=utils.nfft,
                                      transformed_model=self.time_domain_source_model,
                                      transformed_model_data_points=self.time_array,
                                      amplitude_difference=self.amplitude_difference,
                                      phase_difference=self.phase_difference)
    
    def time_domain_strain(self, parameters=None):
        return self._calculate_strain(model=self.time_domain_source_model,
                                      model_data_points=self.time_array,
                                      parameters=parameters,
                                      transformation_function=utils.infft,
                                      transformed_model=self.frequency_domain_source_model,
                                      transformed_model_data_points=self.frequency_array)

    def _calculate_strain(self, model, model_data_points, transformation_function, transformed_model,
                          transformed_model_data_points, parameters, amplitude_difference, phase_difference):
        if parameters is not None:
            self.parameters = parameters
        if self.parameters == self._cache['parameters'] and self._cache['model'] == model and \
                self._cache['transformed_model'] == transformed_model:
            return self._cache['waveform']
        if model is not None:
            model_strain = self._strain_from_model(model_data_points, model)
        elif transformed_model is not None:
            model_strain = self._strain_from_transformed_model(transformed_model_data_points, transformed_model,
                                                               transformation_function)
        else:
            raise RuntimeError("No source model given")
        self._cache['waveform'] = model_strain
        self._cache['parameters'] = self.parameters.copy()
        self._cache['model'] = model
        self._cache['transformed_model'] = transformed_model
        
        '''
        The following block performs the waveform uncertainty correction:
        '''
        
        if self.amplitude_difference is not None:
            model_strain['plus'] *= (1+self.amplitude_difference)*(2+self.phase_difference*1j)/(2-self.phase_difference*1j)
            model_strain['cross'] *= (1+self.amplitude_difference)*(2+self.phase_difference*1j)/(2-self.phase_difference*1j)
            
        return model_strain
    
    def _strain_from_model(self, model_data_points, model):
        return model(model_data_points, **self.parameters)
    
    def _strain_from_transformed_model(self, transformed_model_data_points, transformed_model, transformation_function):
        transformed_model_strain = self._strain_from_model(transformed_model_data_points, transformed_model)
        if isinstance(transformed_model_strain, np.ndarray):
            return transformation_function(transformed_model_strain, self.sampling_frequency)
        model_strain = dict()
        for key in transformed_model_strain:
            if transformation_function == utils.nfft:
                model_strain[key], _ = \
                    transformation_function(transformed_model_strain[key], self.sampling_frequency)
            else:
                model_strain[key] = transformation_function(transformed_model_strain[key], self.sampling_frequency)
        return model_strain
    
    @property
    def parameters(self):
        return self.__parameters
    
    @parameters.setter 
    def parameters(self, parameters):
        if not isinstance(parameters, dict):
            raise TypeError('"parameters" must be a dictionary.')
        new_parameters = parameters.copy()
        new_parameters, _ = self.parameter_conversion(new_parameters)
        self.__parameters = new_parameters
        self.__parameters.update(self.waveform_arguments)
        
    def __parameters_from_source_model(self):
        if self.frequency_domain_source_model is not None:
            model = self.frequency_domain_source_model
        elif self.time_domain_source_model is not None:
            model = self.time_domain_source_model
        else:
            raise AttributeError('Either time or frequency domain source '
                                 'model must be provided.')
        return set(utils.infer_parameters_from_function(model))



class WaveformGeneratorWFU_interp(object):
    
    '''
        Temporary waveform generator; uses 1D interpolation instead of cubic splines
    '''
    
    duration = PropertyAccessor('_times_and_frequencies', 'duration')
    sampling_frequency = PropertyAccessor('_times_and_frequencies', 'sampling_frequency')
    start_time = PropertyAccessor('_times_and_frequencies', 'start_time')
    frequency_array = PropertyAccessor('_times_and_frequencies', 'frequency_array')
    time_array = PropertyAccessor('_times_and_frequencies', 'time_array')
    def __init__(self, duration=None, sampling_frequency=None, start_time=0, frequency_domain_source_model=None,
                 time_domain_source_model=None, parameters=None,
                 waveform_uncertainty_nodes=None,dA_sampling=False,dphi_sampling=False,
                 parameter_conversion=None,
                 waveform_arguments=None):
        self._times_and_frequencies = CoupledTimeAndFrequencySeries(duration=duration,
                                                                    sampling_frequency=sampling_frequency,
                                                                    start_time=start_time)
        self.frequency_domain_source_model = frequency_domain_source_model
        self.time_domain_source_model = time_domain_source_model
        self.source_parameter_keys = self.__parameters_from_source_model()
        self.dA_sampling = dA_sampling
        self.dphi_sampling=dphi_sampling
        
        if parameter_conversion is None:
            self.parameter_conversion = convert_to_lal_binary_black_hole_parameters
        else:
            self.parameter_conversion = parameter_conversion
        if waveform_arguments is not None:
            self.waveform_arguments = waveform_arguments
        else:
            self.waveform_arguments = dict()
        if waveform_uncertainty_nodes is not None:
            self.waveform_uncertainty_nodes = waveform_uncertainty_nodes
        else:
            self.waveform_uncertainty_nodes = None
        if isinstance(parameters, dict):
            self.parameters = parameters
        self._cache = dict(parameters=None, waveform=None, model=None)
        
        utils.logger.info(
            "Waveform generator initiated with\n"
            "  frequency_domain_source_model: {}\n"
            "  time_domain_source_model: {}\n"
            "  parameter_conversion: {}"
            .format(utils.get_function_path(self.frequency_domain_source_model),
                    utils.get_function_path(self.time_domain_source_model),
                    utils.get_function_path(self.parameter_conversion))
        )
        
    def __repr__(self):
        if self.frequency_domain_source_model is not None:
            fdsm_name = self.frequency_domain_source_model.__name__
        else:
            fdsm_name = None
        if self.time_domain_source_model is not None:
            tdsm_name = self.time_domain_source_model.__name__
        else:
            tdsm_name = None
        if self.parameter_conversion is None:
            param_conv_name = None
        else:
            param_conv_name = self.parameter_conversion.__name__
        return self.__class__.__name__ + '(duration={}, sampling_frequency={}, start_time={}, ' \
                                         'frequency_domain_source_model={}, time_domain_source_model={}, ' \
                                         'parameter_conversion={}, ' \
                                         'waveform_uncertainty_nodes={}, ' \
                                         'dA_sampling={}, ' \
                                         'dphi_sampling={}, ' \
                                         'waveform_arguments={})'\
            .format(self.duration, self.sampling_frequency, self.start_time, fdsm_name, tdsm_name,
                    param_conv_name,self.waveform_uncertainty_nodes, self.waveform_arguments)
    
    def frequency_domain_strain(self, parameters=None):
        return self._calculate_strain(model=self.frequency_domain_source_model,
                                      model_data_points=self.frequency_array,
                                      parameters=parameters,
                                      transformation_function=utils.nfft,
                                      transformed_model=self.time_domain_source_model,
                                      transformed_model_data_points=self.time_array,
                                      waveform_uncertainty_nodes=self.waveform_uncertainty_nodes)
    
    def time_domain_strain(self, parameters=None):
        return self._calculate_strain(model=self.time_domain_source_model,
                                      model_data_points=self.time_array,
                                      parameters=parameters,
                                      transformation_function=utils.infft,
                                      transformed_model=self.frequency_domain_source_model,
                                      transformed_model_data_points=self.frequency_array)

    def _calculate_strain(self, model, model_data_points, transformation_function, transformed_model,
                          transformed_model_data_points, parameters, waveform_uncertainty_nodes):
        if parameters is not None:
            self.parameters = parameters
        if self.parameters == self._cache['parameters'] and self._cache['model'] == model and \
                self._cache['transformed_model'] == transformed_model:
            return self._cache['waveform']
        if model is not None:
            model_strain = self._strain_from_model(model_data_points, model)
        elif transformed_model is not None:
            model_strain = self._strain_from_transformed_model(transformed_model_data_points, transformed_model,
                                                               transformation_function)
        else:
            raise RuntimeError("No source model given")
        self._cache['waveform'] = model_strain
        self._cache['parameters'] = self.parameters.copy()
        self._cache['model'] = model
        self._cache['transformed_model'] = transformed_model
        
        '''
        The following block performs the waveform uncertainty correction:
        '''
        
        if self.waveform_uncertainty_nodes is not None:
            alphas = np.zeros(len(self.waveform_uncertainty_nodes))
            betas = np.zeros(len(self.waveform_uncertainty_nodes))
            for i in range(len(self.waveform_uncertainty_nodes)):
                if self.dA_sampling == True:
                    alphas[i]=parameters['alpha_' + f'{i+1}']
                if self.dphi_sampling == True:
                    betas[i]=parameters['beta_' + f'{i+1}']
            dA = scipy.interpolate.interp1d(self.waveform_uncertainty_nodes,alphas,fill_value="extrapolate")(self.frequency_array)
            dphi = scipy.interpolate.interp1d(self.waveform_uncertainty_nodes,betas,fill_value="extrapolate")(self.frequency_array)
            model_strain['plus'] *= (1+dA)*(2+dphi*1j)/(2-dphi*1j)
            model_strain['cross'] *= (1+dA)*(2+dphi*1j)/(2-dphi*1j)
            
        return model_strain
    
    def _strain_from_model(self, model_data_points, model):
        return model(model_data_points, **self.parameters)
    
    def _strain_from_transformed_model(self, transformed_model_data_points, transformed_model, transformation_function):
        transformed_model_strain = self._strain_from_model(transformed_model_data_points, transformed_model)
        if isinstance(transformed_model_strain, np.ndarray):
            return transformation_function(transformed_model_strain, self.sampling_frequency)
        model_strain = dict()
        for key in transformed_model_strain:
            if transformation_function == utils.nfft:
                model_strain[key], _ = \
                    transformation_function(transformed_model_strain[key], self.sampling_frequency)
            else:
                model_strain[key] = transformation_function(transformed_model_strain[key], self.sampling_frequency)
        return model_strain
    
    @property
    def parameters(self):
        return self.__parameters
    
    @parameters.setter 
    def parameters(self, parameters):
        if not isinstance(parameters, dict):
            raise TypeError('"parameters" must be a dictionary.')
        new_parameters = parameters.copy()
        new_parameters, _ = self.parameter_conversion(new_parameters)
        self.__parameters = new_parameters
        self.__parameters.update(self.waveform_arguments)
        
    def __parameters_from_source_model(self):
        if self.frequency_domain_source_model is not None:
            model = self.frequency_domain_source_model
        elif self.time_domain_source_model is not None:
            model = self.time_domain_source_model
        else:
            raise AttributeError('Either time or frequency domain source '
                                 'model must be provided.')
        return set(utils.infer_parameters_from_function(model))
