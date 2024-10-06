"WaveformUncertainty package"
__version__ = "0.7.0.3"

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
    polarization: string, optional
        polarization of the strain data (plus or cross)
        default: 'plus'
    psd_data: numpy.ndarray, optional
        array containing the psd data and their corresponding frequencies
        default: None
    correction_parameter_A: float, optional
        value at which to cut the second derivative of amplitude difference; if None, correction will not occur
        default: 5e-5
    correction_parameter_B: int, optional
        index at which to start the search for any discontinuity
        default: 0
    correction_parameter_C: int, optional
        number of dA derivatives to take for discontinuity correction
        default: 2
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
        array of phase difference values; if psd data is None, raw_phase_difference will be returned, residual_phase_difference otherwise
    amplitude_difference_final_point: float
        value of the amplitude difference where the discontinuity correction starts
    phase_difference_final_point: float
        value of the phase difference where the discontinuity correction starts
    final_index: int
        position within the frequency grid where the discontinuity correction starts
    '''
    injection = kwargs.get('injection',None)
    npoints = kwargs.get('npoints',1000)
    polarization = kwargs.get('polarization','plus')
    psd_data = kwargs.get('psd_data',None)
    correction_parameter_A = kwargs.get('correction_parameter_A',5e-5)
    correction_parameter_B = kwargs.get('correction_parameter_B',0)
    correction_parameter_C = kwargs.get('correction_parameter_C',2)
    ref_amplitude = kwargs.get('ref_amplitude',None)

    f_low = hf1.waveform_arguments['f_low']
    f_high = hf1.waveform_arguments['f_high']
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

    # waveform amplitudes
    amplitude_1 = np.abs(hf1.frequency_domain_strain()[f'{polarization}'][wf_freqs])
    amplitude_2 = np.abs(hf2.frequency_domain_strain()[f'{polarization}'][wf_freqs])

    # waveform phases
    phase_1 = np.angle(hf1.frequency_domain_strain()[f'{polarization}'][wf_freqs])
    phase_2 = np.angle(hf2.frequency_domain_strain()[f'{polarization}'][wf_freqs])
                     
    amplitude_difference = (amplitude_2-amplitude_1)/amplitude_1
    raw_phase_difference = phase_2-phase_1

    # removing phase shifts of 2pi
    while any(value > 6 for value in [np.abs(raw_phase_difference[i+1]-raw_phase_difference[i]) for i in range(len(raw_phase_difference)-2)]):
        raw_phase_difference = np.unwrap(raw_phase_difference)

    # fitting a line to raw_phase_difference weighted by PSDs and subtracting off that line
    if psd_data is not None:
        if ref_amplitude is None:
            ref_amplitude = np.abs(hf1.frequency_domain_strain()[f'{polarization}'][wf_freqs])
        ref_amplitude = np.interp(hf1.frequency_array[wf_freqs],f_high*np.linspace(0,1,len(ref_amplitude)),ref_amplitude)
        ref_sigma = np.interp(hf1.frequency_array[wf_freqs], psd_data[:,0],psd_data[:,1])
        align_weights = ref_amplitude*ref_amplitude / ref_sigma * hf1.frequency_array[wf_freqs]
        fit = np.polyfit(hf1.frequency_array[wf_freqs],raw_phase_difference,1,w=align_weights)
        residual_phase_difference = raw_phase_difference - np.poly1d(fit)(hf1.frequency_array[wf_freqs])

    # taking two frequency derivatives of amplitude_difference and comparing it to the correction parameter for f_COR calculation
    amplitude_difference_derivative = np.copy(amplitude_difference)
    for i in range(correction_parameter_C):
        amplitude_difference_derivative = np.gradient(amplitude_difference_derivative)/np.gradient(frequency_grid)

    for i in range(correction_parameter_B,len(amplitude_difference_derivative)):
        if correction_parameter_A is None:
            final_index = len(hf1.frequency_array[wf_freqs])-1
            break
        elif np.abs(amplitude_difference_derivative[i]) > correction_parameter_A:
            final_index = i
            break
        else: 
            final_index = len(hf1.frequency_array[wf_freqs])-1

    # making the discontinuity correction to amplitude_difference
    amplitude_difference[final_index:] = amplitude_difference[final_index-1]
    amplitude_difference_final_point = amplitude_difference[final_index-1]

    # making the discontinuity correction to phase difference (raw or residual)
    if psd_data is not None:
        phase_difference = np.copy(residual_phase_difference)
        phase_difference[final_index:] = residual_phase_difference[final_index-1]
        phase_difference_final_point = residual_phase_difference[final_index-1]
    else:
        phase_difference = np.copy(raw_phase_difference)
        phase_difference[final_index:] = raw_phase_difference[final_index-1]
        phase_difference_final_point = raw_phase_difference[final_index-1]

    return frequency_grid,amplitude_difference,phase_difference,amplitude_difference_final_point,phase_difference_final_point,final_index

    
    
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
        # reconstructing amplitude difference up to the discontinuity
        parameterized_curve = np.polynomial.chebyshev.chebval(data[1][0:data[4]],data[2])
        # reconstructing amplitude difference after the discontinuity
        post_waveform_uncertainties = np.copy(data[1])
        post_waveform_uncertainties[data[4]-2:] = data[5]
    
    elif str(identity) == 'phase_difference':
        # reconstructing phase difference up to the discontinuity
        parameterized_curve = np.polynomial.chebyshev.chebval(data[1][0:data[4]],data[3])
        # reconstructing phase difference after the discontinuity
        post_waveform_uncertainties = np.copy(data[1])
        post_waveform_uncertainties[data[4]-2:] = data[6]
    
    else:
        raise Exception('Identity of the uncertainty must be "amplitude_difference" or "phase_difference".')

    # combining both sides of the discontinuity
    difference_array = np.concatenate((parameterized_curve[:-2],post_waveform_uncertainties[data[4]-2:]))
    
    return difference_array



def parameterization(hf1,hf2,parameter_data,nsamples,**kwargs):
    '''
    Generates samples of waveform uncertainty between two approximants and parameterizes the data with Chebyshev polynomial functions.

    Parameters
    ==================
    hf1: bilby.gw.waveform_generator.WaveformGenerator
        frequency domain waveform generator object
    hf2: bilby.gw.waveform_generator.WaveformGenerator
        frequency domain waveform generator object
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
    max_dA_error: float [%], optional
        maximum allowed error between the amplitude uncertainty and its parameterization
        default: 1 
    max_dphi_error: float [degrees], optional
        maximum allowed error between the phase uncertainty and its parameterization
        default: 5
    psd_data: numpy.ndarray, optional
        array containing the psd data and their corresponding frequencies
        default: None
    correction_parameter_A: float, optional
        value at which to cut the second derivative of amplitude difference; if None, correction will not occur
        default: 5e-5
    correction_parameter_B: int, optional
        index at which to start the search for any discontinuity
        default: 0
    correction_parameter_C: int, optional
        number of amplitude difference derivatives to take for the discontinuity correction
        default: 2
    ref_amplitude: numpy.ndarray, optional
        reference amplitude for residual phase calculation
        default: None
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
    npoints = kwargs.get('npoints',1000)
    polarization = kwargs.get('polarization','plus')
    psd_data = kwargs.get('psd_data',None)
    correction_parameter_A = kwargs.get('correction_parameter_A',5e-5)
    correction_parameter_B = kwargs.get('correction_parameter_B',0)
    correction_parameter_C = kwargs.get('correction_parameter_C',2)
    ref_amplitude = kwargs.get('ref_amplitude',None)
    precession = kwargs.get('precession',False)
    tides = kwargs.get('tides',True)
    fit_parameters = kwargs.get('fit_parameters',15)
    max_amplitude_error = kwargs.get('max_amplitude_error',1)
    max_phase_error = kwargs.get('max_phase_error',5)
    fit_threshold = kwargs.get('fit_threshold',75)
    
    np.seterr(all='ignore')
    progress = 1
    bilby.core.utils.log.setup_logger(log_level=30)
    
    if type(parameter_data)==bilby.core.prior.dict.PriorDict:
        # converting a prior object to a dictionary
        data = parameter_dict_from_prior(parameter_data,int((nsamples/fit_threshold)+1))
    else:
        data = parameter_data

    # generating random order of samples
    index_samples=list(range(len(data[list(data.keys())[0]])))
    indexes=[]
    for draws in range(len(index_samples)):
        index=random.choice(index_samples)
        indexes.append(index)
        index_samples.remove(index)
    
    def progressBar(count_value, total, suffix=''): #To Do
        bar_length = 100
        filled_up_Length = int(round(bar_length* count_value / float(total)))
        percentage = round(100.00 * count_value/float(total),1)
        bar = '=' * filled_up_Length + '-' * (bar_length - filled_up_Length)
        sys.stdout.write('[%s] %s%s %s\r' %(bar, percentage, '%', suffix))
        sys.stdout.flush()

    # setting initial conditions
    trial = 0
    final_indexes = []
    extraneous_indexes = []
    
    output_matrix = np.zeros([len(indexes),8],dtype=object)
    start = tm.time()
    
    print("Generating Waveform Differences and Parameterizing...")

    # setting the reference amplitude
    if ref_amplitude is None:
        ref_amplitude = np.abs(hf1.frequency_domain_strain(parameters=injection(parameter_data,precession=precession,tides=tides))[f'{polarization}'])
    
    for index in indexes:
    
        progressBar(progress,(nsamples))

        # calculating waveform model differences
        frequency_grid,amplitude_difference,phase_difference,amplitude_difference_final_point,phase_difference_final_point,final_index = fd_model_difference(hf1,hf2,injection=injection(parameter_data,index=index,precession=precession,tides=tides),npoints=npoints,polarization=polarization,psd_data=psd_data,correction_parameter_A=correction_parameter_A,correction_parameter_B=correction_parameter_B,correction_parameter_C=correction_parameter_C,ref_amplitude=ref_amplitude)

        # chebyshev polynomial fits and saving coefficients
        amplitude_difference_fit = np.polynomial.chebyshev.Chebyshev.fit((frequency_grid[0:final_index]),amplitude_difference[0:final_index],fit_parameters-1)
        amplitude_difference_parameters = amplitude_difference_fit.convert().coef

        phase_difference_fit = np.polynomial.chebyshev.Chebyshev.fit((frequency_grid[0:final_index]),phase_difference[0:final_index],fit_parameters-1)
        phase_difference_parameters = phase_difference_fit.convert().coef

        #constructing output matrix
        output_matrix[index][0] = index
        output_matrix[index][1] = np.array(frequency_grid)
        output_matrix[index][2] = np.array(amplitude_difference_parameters)
        output_matrix[index][3] = np.array(phase_difference_parameters)
        output_matrix[index][4] = final_index
        output_matrix[index][5] = amplitude_difference_final_point
        output_matrix[index][6] = phase_difference_final_point
        output_matrix[index][7] = injection(data,index=index,precession=precession,tides=tides)

        final_indexes.append(index)

        # finding errors relative to error margins
        relative_amplitude_error = np.max(np.abs(100*((amplitude_difference-recovery_from_parameterization('amplitude_difference',output_matrix[index])))/([max_amplitude_error]*len(amplitude_difference))))
        relative_phase_error = np.max(np.abs((100*((phase_difference-recovery_from_parameterization('phase_difference',output_matrix[index])))/([max_phase_error*(2*np.pi/360)]*len(phase_difference)))))

        # removing index if outside error margins
        if relative_amplitude_error>100 or relative_phase_error>100:
            final_indexes.remove(index)
            extraneous_indexes.append(index)
        else:
            progress += 1

        # checking if parameterization rate is above fit threshold
        if 100*len(final_indexes)/(len(final_indexes)+len(extraneous_indexes)) < fit_threshold:
            trial += 1

        # allows for 20 failures before killing the run
        if trial == 20:
            parameterization_rate = np.round(100*len(final_indexes)/(len(final_indexes)+len(extraneous_indexes)),2)
            raise Exception(f"Parameterization Rate Too Low; It must be above {fit_threshold}%, but was at {parameterization_rate}%. {len(final_indexes)+len(extraneous_indexes)}")
        
        if len(final_indexes) == nsamples:
            break

    # constructing the final parameterization table
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
    print(f"Parameterized Waveform Model Difference data was created with {len(final_indexes)} sample sets at a final parameterization rate of {final_parameterization_rate}%.") 
    print("")
    
    return parameterized_data



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

    # grabbing all of the waveform difference data
    amplitude_difference_data = [recovery_from_parameterization('amplitude_difference',data[index]) for index in range(len(data))]
    phase_difference_data = [recovery_from_parameterization('phase_difference',data[index]) for index in range(len(data))]

    if linear==True and resolution is not None:

        # constructing linear frequency grid
        linear_frequency_grid = np.arange(data[0][1][0],data[0][1][-1]+resolution,resolution)

        # interpolating waveform differences to the new linear frequency grid
        amplitude_difference = [np.interp(linear_frequency_grid,data[0][1],amplitude_difference_data[i]) for i in range(len(amplitude_difference_data))]
        phase_difference = [np.interp(linear_frequency_grid,data[0][1],phase_difference_data[i]) for i in range(len(phase_difference_data))]

        # mean and standard devitation along the vertical axis
        mean_amplitude_difference = np.mean(amplitude_difference,axis=0)
        amplitude_uncertainty = np.std(amplitude_difference,axis=0)

        mean_phase_difference = np.mean(phase_difference,axis=0)
        phase_uncertainty = np.std(phase_difference,axis=0)
    
    else:

        linear_frequency_grid = None

        # mean and standard devitation along the vertical axis
        mean_amplitude_difference = np.mean(amplitude_difference_data,axis=0)
        amplitude_uncertainty = np.std(amplitude_difference_data,axis=0)

        mean_phase_difference = np.mean(phase_difference_data,axis=0)
        phase_uncertainty = np.std(phase_difference_data,axis=0)
        
    return mean_amplitude_difference,amplitude_uncertainty,mean_phase_difference,phase_uncertainty,linear_frequency_grid



def match(signal,data,PSDs,duration):
    
    signal_match = np.sqrt(bilby.gw.utils.matched_filter_snr(signal,signal,PSDs,4))
    data_match = np.sqrt(bilby.gw.utils.matched_filter_snr(data,data,PSDs,4))
    normalized_match = np.abs(bilby.gw.utils.matched_filter_snr(signal,data,PSDs,4)/(signal_match*data_match))
    
    return normalized_match



def WFU_dphi_prior(phase_uncertainty,frequency_grid,injection,hf,PSDs,match_boundary,duration,nnodes,**kwargs):
    prior = kwargs.get('prior',None)
    zero_resolution = kwargs.get('zero_resolution',10)
    polarization = kwargs.get('polarization','plus')
    match_resolution = kwargs.get('match_resolution',100)

    if prior is None:
        prior = bilby.core.prior.PriorDict()

    M = bilby.gw.conversion.generate_mass_parameters(injection)['total_mass']*lal.MSUN_SI
    f_light = (lal.C_SI**3)/(lal.G_SI*M)
    f_IM = 0.018*f_light
    
    low_frequency_nodes = np.arange(frequency_grid[0],f_IM,(f_IM-frequency_grid[0])/zero_resolution).astype(int)
    frequency_nodes = list(np.geomspace(f_IM,f_light,nnodes).astype(int))
    frequency_nodes.append(frequency_grid[-1])
    
    total_frequency_nodes = np.concatenate((low_frequency_nodes,frequency_nodes))
    indexes = list(range(-zero_resolution+1,nnodes+2))
    hf.frequency_nodes = total_frequency_nodes
    hf.indexes = indexes
    
    position = []
    for node in total_frequency_nodes:
        position.append(list(np.round(frequency_grid,1)).index(float(node)))

    for i in indexes:
        injection[f'dphi_{i}'] = 0
    reference_waveform = hf.frequency_domain_strain(parameters=injection)[polarization]
    
    for i in range(nnodes):
        good_dphis = []
        for j in np.linspace(0,5*phase_uncertainty[position[zero_resolution+i]],match_resolution):
            new_injection = injection.copy()
            new_injection[f'dphi_{i+1}'] = j
            waveform = hf.frequency_domain_strain(parameters=new_injection)[polarization]
            match_percent = match(reference_waveform,waveform,PSDs,duration)*100
            if match_percent >= match_boundary:
                good_dphis.append(j)
        prior[f'dphi_{i+1}'] = bilby.core.prior.TruncatedGaussian(name=f'dphi_{i+1}',latex_label=r'$\varphi_{num}$'.replace('num',str(i+1)),
                                                                mu=0,sigma=phase_uncertainty[position[zero_resolution+i]],
                                                                minimum=-good_dphis[-1],maximum=good_dphis[-1])
        print(f'dphi_{i+1} Prior Complete')
    
    prior[f'dphi_{nnodes+1}'] = bilby.core.prior.DeltaFunction(name=f'dphi_{nnodes+1}',peak=0)
    
    for i in range(zero_resolution):
        prior[f'dphi_{-i}'] = bilby.core.prior.DeltaFunction(name=f'dphi_{-i}',peak=0)
    
    return prior,total_frequency_nodes



def WFU_amplitude_prior(amplitude_uncertainty,frequency_grid,nnodes,**kwargs):
    f_M = kwargs.get('f_M',100)
    prior = kwargs.get('prior',None)
    zero_resolution = kwargs.get('zero_resolution',10)
    
    if prior is None:
        prior = bilby.core.prior.PriorDict()
    
    low_frequency_nodes = np.arange(frequency_grid[0],f_M,(f_M-frequency_grid[0])/zero_resolution)
    frequency_nodes = np.geomspace(f_M,frequency_grid[-1],nnodes).astype(int)
    
    total_frequency_nodes = np.concatenate((low_frequency_nodes,frequency_nodes))
    
    position = []
    for node in frequency_nodes:
        position.append(list(np.round(frequency_grid,1)).index(float(node)))
    
    for i in range(nnodes):
        prior[f'dA_{i+1}'] = bilby.core.prior.Gaussian(name=f'dA_{i+1}',latex_label=r'$\alpha_{val}$'.replace('val',str(i+1)),mu=0,sigma=amplitude_uncertainty[position[i]])
        
    prior[f'dA_{nnodes}'] = bilby.core.prior.DeltaFunction(name=f'dA_{nnodes}',peak=0)
    
    for i in range(zero_resolution):
        prior[f'dA_{-i}'] = bilby.core.prior.DeltaFunction(name=f'dA_{-i}',peak=0)
    
    return prior,total_frequency_nodes



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



def Q_factor(WFU_result,NHP_result,injection):
    WFU_maxL = maxL(WFU_result)
    NHP_maxL = maxL(NHP_result)

    try:
        if WFU_maxL['phase'] > np.pi:
            WFU_maxL['phase'] -= np.pi
        if NHP_maxL['phase'] > np.pi:
            NHP_maxL['phase'] -= np.pi
    except:
        pass
    
    NHP_sum = 0
    WFU_sum = 0
    
    for p in injection.keys():
        NHP_sum += ((NHP_maxL[p]-injection[p])/injection[p])**2
        WFU_sum += ((WFU_maxL[p]-injection[p])/injection[p])**2
    
    epsilon_NHP = np.sqrt(NHP_sum)
    epsilon_WFU = np.sqrt(WFU_sum)
    
    Q_factor = ((epsilon_NHP-epsilon_WFU)/epsilon_NHP)*100
    
    return Q_factor



class WaveformGeneratorWFU(object):
    '''
    Modified WaveformGenerator object from bilby.gw to include waveform uncertainty corrections in the strain calculation
    To sample waveform uncertainty, include all relevant "alpha" and "phi" parameters in the prior.
    Note: make sure the number of alphas, phis, and waveform_uncertainty_nodes are the same
    
    New Parameters
    ==================
    waveform_uncertainty_nodes: numpy.ndarray, optional
        array of frequency nodes to be used in generating the dA and dphi splines
        default: None
    dA_sampling: bool, optional
        if True, the waveform generator will attempt to pull alpha parameters from the parameter dictionary (either an injection or the prior)
        default: None
    dphi_sampling: bool, optional
        if True, the waveform generator will attempt to pull phi parameters from the parameter dictionary (either an injection or the prior)
        default: None
    indexes: numpy.array, optional
        the numbers of waveform correction parameters; e.g. phi_1 = 1, alpha_-8 = -8, etc.
        default: None
    '''
    duration = PropertyAccessor('_times_and_frequencies', 'duration')
    sampling_frequency = PropertyAccessor('_times_and_frequencies', 'sampling_frequency')
    start_time = PropertyAccessor('_times_and_frequencies', 'start_time')
    frequency_array = PropertyAccessor('_times_and_frequencies', 'frequency_array')
    time_array = PropertyAccessor('_times_and_frequencies', 'time_array')
    def __init__(self, duration=None, sampling_frequency=None, start_time=0, frequency_domain_source_model=None,
                 time_domain_source_model=None, parameters=None,
                 frequency_nodes=None,indexes=None,
                 parameter_conversion=None,
                 waveform_arguments=None):
        self._times_and_frequencies = CoupledTimeAndFrequencySeries(duration=duration,
                                                                    sampling_frequency=sampling_frequency,
                                                                    start_time=start_time)
        self.frequency_domain_source_model = frequency_domain_source_model
        self.time_domain_source_model = time_domain_source_model
        self.source_parameter_keys = self.__parameters_from_source_model()
        self.indexes=indexes
        
        if parameter_conversion is None:
            self.parameter_conversion = convert_to_lal_binary_black_hole_parameters
        else:
            self.parameter_conversion = parameter_conversion
        if waveform_arguments is not None:
            self.waveform_arguments = waveform_arguments
        else:
            self.waveform_arguments = dict()
        if frequency_nodes is not None:
            self.frequency_nodes = frequency_nodes
        else:
            self.frequency_nodes = None
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
                                         'frequency_nodes={}, ' \
                                         'indexes={}, ' \
                                         'waveform_arguments={})'\
            .format(self.duration, self.sampling_frequency, self.start_time, fdsm_name, tdsm_name,
                    param_conv_name, self.frequency_nodes, self.indexes, self.waveform_arguments)
    
    def frequency_domain_strain(self, parameters=None):
        return self._calculate_strain(model=self.frequency_domain_source_model,
                                      model_data_points=self.frequency_array,
                                      parameters=parameters,
                                      transformation_function=utils.nfft,
                                      transformed_model=self.time_domain_source_model,
                                      transformed_model_data_points=self.time_array,
                                      frequency_nodes=self.frequency_nodes,
                                      indexes=self.indexes
                                      )

    def time_domain_strain(self,parameters=None):
        fd_model_strain = self._calculate_strain(model=self.frequency_domain_source_model,
                                      model_data_points=self.frequency_array,
                                      parameters=parameters,
                                      transformation_function=utils.nfft,
                                      transformed_model=self.time_domain_source_model,
                                      transformed_model_data_points=self.time_array,
                                      frequency_nodes=self.frequency_nodes,
                                      indexes=self.indexes
                                      )
        
        plus_td_waveform = self.sampling_frequency*np.fft.ifft(fd_model_strain['plus'])
        plus_td_model_strain = np.interp(self.time_array,np.linspace(self.time_array[0],self.time_array[-1],len(plus_td_waveform)),plus_td_waveform)

        cross_td_waveform = self.sampling_frequency*np.fft.ifft(fd_model_strain['cross'])
        cross_td_model_strain = np.interp(self.time_array,np.linspace(self.time_array[0],self.time_array[-1],len(cross_td_waveform)),cross_td_waveform)

        model_strain = dict()
        model_strain['plus'] = plus_td_model_strain
        model_strain['cross'] = cross_td_model_strain

        return model_strain
    
    def _calculate_strain(self, model, model_data_points, transformation_function, transformed_model,
                          transformed_model_data_points, parameters, frequency_nodes,indexes):
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
        
        if self.frequency_nodes is not None:
            if self.indexes is not None:
                indexes = self.indexes
            else:
                indexes = list(range(1,len(self.frequency_nodes)+1))
            
            try:
                alphas = [parameters[f'dA_{i}'] for i in indexes]
                dA = scipy.interpolate.CubicSpline(self.frequency_nodes,alphas)(self.frequency_array)
            except:
                dA = 0
           
            try:
                phis = [parameters[f'dphi_{i}'] for i in indexes]
                dphi = scipy.interpolate.CubicSpline(self.frequency_nodes,phis)(self.frequency_array)
            except:
                dphi = 0
                
            model_strain['plus'] = model_strain['plus']*(1+dA)*np.exp(dphi*1j)
            model_strain['cross'] = model_strain['cross']*(1+dA)*np.exp(dphi*1j)
            
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
        for key in self.source_parameter_keys.symmetric_difference(new_parameters):
            # preventing waveform uncertainty parameters from being removed
            if self.frequency_nodes is not None:
                if self.indexes is not None:
                    indexes = self.indexes
                else:
                    indexes = list(range(1,len(self.frequency_nodes)+1))
                if key not in [f'dA_{i}' for i in indexes]+[f'dphi_{i}' for i in indexes]:  
                    new_parameters.pop(key)
            else:
                new_parameters.pop(key)
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
