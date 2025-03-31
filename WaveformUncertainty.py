"WaveformUncertainty package"
__version__ = "0.9.6"

import numpy as np
import bilby
import random
import time as tm
import sys
import scipy
import lal
import matplotlib.pyplot as plt
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
    correction_parameter = kwargs.get('correction_parameter',0.0001)
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

    # finding the position of the cuttoff/discontinuity correction frequency, f_disc
    final_index_1 = list(amplitude_1).index(min(amplitude_1, key=lambda x:np.abs(x-correction_parameter*np.max(amplitude_1))))
    final_index_2 = list(amplitude_2).index(min(amplitude_2, key=lambda x:np.abs(x-correction_parameter*np.max(amplitude_2))))
    final_index = min([final_index_1,final_index_2])
    
    # fitting a line to raw_phase_difference weighted by PSDs and subtracting off that line
    if psd_data is not None:
        if ref_amplitude is None:
            ref_amplitude = np.abs(hf1.frequency_domain_strain()[f'{polarization}'][wf_freqs][0:final_index])
        ref_amplitude = np.interp(hf1.frequency_array[wf_freqs][0:final_index],f_high*np.linspace(0,1,len(ref_amplitude)),ref_amplitude)
        ref_sigma = np.interp(hf1.frequency_array[wf_freqs][0:final_index], psd_data[:,0],psd_data[:,1])
        align_weights = (ref_amplitude**2)/(ref_sigma*hf1.frequency_array[wf_freqs][0:final_index])
        fit = np.polyfit(hf1.frequency_array[wf_freqs][0:final_index],raw_phase_difference[0:final_index],1,w=align_weights)
        raw_phase_difference_no_shifts = raw_phase_difference[0:final_index]-np.poly1d(fit)(hf1.frequency_array[wf_freqs][0:final_index])
        residual_phase_difference = np.copy(raw_phase_difference)
        residual_phase_difference[0:final_index] = raw_phase_difference_no_shifts
        
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
    correction_parameter: float, optional
        fraction of maximum amplitude to cut off the amplitude at
        default: 0.0001
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
    correction_parameter = kwargs.get('correction_parameter',0.0001)
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

    data = prior.sample(int((100*nsamples/fit_threshold)+1))

    # generating random order of samples
    index_samples=list(range(len(data[list(data.keys())[0]])))
    indexes=[]
    for draws in range(len(index_samples)):
        index=random.choice(index_samples)
        indexes.append(index)
        index_samples.remove(index)
    
    # setting initial conditions
    trial = 0
    final_indexes = []
    extraneous_indexes = []
    
    output_matrix = np.zeros([len(indexes),8],dtype=object)
    start = tm.time()
    
    print("Generating Waveform Differences and Parameterizing...")

    # setting the reference amplitude
    if ref_amplitude is None:
        ref_amplitude = np.abs(hf1.frequency_domain_strain(parameters=injection(data))[f'{polarization}'])
    
    for index in indexes:
    
        progressBar(progress,(nsamples))

        # calculating waveform model differences
        frequency_grid,amplitude_difference,phase_difference,amplitude_difference_final_point,phase_difference_final_point,final_index = fd_model_difference(hf1,hf2,injection=injection(data,index=index),npoints=npoints,polarization=polarization,psd_data=psd_data,correction_parameter=correction_parameter,ref_amplitude=ref_amplitude)

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
        output_matrix[index][7] = injection(data,index=index)

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
    geometrized_frequency_grid: numpy.ndarray, optional
        if given, data will be returned in geometrized units with points corresponding to this array
        default: None
    resolution: float, optional
        size of output frequency grid if no geometrized_frequency_grid is given
        default: 5000
        
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
    new_frequency_grid: numpy.ndarray
        output frequency grid
    '''
    geometrized_frequency_grid = kwargs.get('geometrized_frequency_grid',None)
    resolution = kwargs.get('resolution',5000)

    # grabbing all of the waveform difference data
    amplitude_difference_data = [recovery_from_parameterization('amplitude_difference',data[i]) for i in range(len(data))]
    phase_difference_data = [recovery_from_parameterization('phase_difference',data[i]) for i in range(len(data))]

    # constructing linear frequency grid
    new_frequency_grid = np.linspace(data[0][1][0],data[0][1][-1],resolution)

    # interpolating waveform differences to the new frequency grid
    amplitude_difference_array = [np.interp(new_frequency_grid,data[0][1],amplitude_difference_data[i]) for i in range(len(amplitude_difference_data))]
    phase_difference_array = [np.interp(new_frequency_grid,data[0][1],phase_difference_data[i]) for i in range(len(phase_difference_data))]

    if geometrized_frequency_grid is not None:
        total_mass_array = np.array([bilby.gw.conversion.generate_mass_parameters(data[i][7])['total_mass']*lal.MSUN_SI for i in range(len(data))])
        geometrized_frequency_array = np.array([new_frequency_grid*total_mass_array[i]*lal.G_SI*(299792458**(-3)) for i in range(len(data))])
        geometrized_amplitude_difference_array = np.array([np.interp(geometrized_frequency_grid,geometrized_frequency_array[i],amplitude_difference_array[i]) for i in range(len(data))])
        geometrized_phase_difference_array = np.array([np.interp(geometrized_frequency_grid,geometrized_frequency_array[i],phase_difference_array[i]) for i in range(len(data))])
        
        amplitude_difference_array = np.copy(geometrized_amplitude_difference_array)
        phase_difference_array = np.copy(geometrized_phase_difference_array)
        new_frequency_grid = np.copy(geometrized_frequency_grid)
    
    # mean and standard devitation along the vertical axis
    mean_amplitude_difference = np.mean(amplitude_difference_array,axis=0)
    amplitude_uncertainty = np.std(amplitude_difference_array,axis=0)

    mean_phase_difference = np.mean(phase_difference_array,axis=0)
    phase_uncertainty = np.std(phase_difference_array,axis=0)
       
    return mean_amplitude_difference,amplitude_uncertainty,mean_phase_difference,phase_uncertainty,new_frequency_grid



def dphi_prior(phase_uncertainty,k, **kwargs):
    '''
    Generates a Gaussian prior for the phase correction parameters
    
    Parameters
    ===================
    phase_uncertainty: numpy.ndarray
        array of standard deviation of a set of phase differences; by default, this should be as a function of dimensionless frequency, xi
    k: int
        number of phase correction parameters desired
    mean_phase_difference: numpy.ndarray, optional
        array of the means of a set of phase differences, by default, this should be as a function of dimensionless frequency, xi
        default: None
    prior: bilby.core.prior.PriorDict, optional
        bilby prior object; if given, dphi priors will be added to this dictionary
        default: None
    geometrized: bool, optional
        if True, will return geometrized frequency nodes; if False, normal frequency nodes (Hz)
        default: True
    xi_low: float, optional
        if geometrized is True; lower bound on the geometrized frequency band
        default: 0.018
    xi_high: float, optional
        if geometrized is True; upper bound on the geometrized frequency band
        default: 1/pi (0.318...)
    f_low: float, optional
        if geometrized is False; lower bound on the normal frequency band
        default: 20.0 Hz
    f_high: float, optional
        if geometrized is False; upper bound on the normal frequency band
        default: 1024.0 Hz
        
    Returns
    ==================
    frequency_nodes: numpy.ndarray
        array of frequency nodes
    prior: bilby.core.prior.PriorDict
        bilby prior object containing the phase correction priors
    '''
    f_low = kwargs.get('f_low',20)
    f_high = kwargs.get('f_high',1024)
    xi_low = kwargs.get('xi_low',0.018)
    xi_high = kwargs.get('xi_high',1/np.pi)
    prior = kwargs.get('prior',None)
    geometrized = kwargs.get('geometrized',True)
    mean_phase_difference = kwargs.get('mean_phase_difference',None)
    
    if prior is None:
        prior = bilby.core.prior.PriorDict()
    
    if mean_phase_difference is None:
        mean_phase_difference = np.array([0]*len(phase_uncertainty))
    
    if geometrized is True:
        frequency_grid = np.linspace(0.001,1,len(phase_uncertainty))
        desired_frequency_nodes = np.geomspace(xi_low,xi_high,k+1)
    else:
        frequency_grid = np.linspace(f_low,f_high,len(phase_uncertainty))
        desired_frequency_nodes = np.geomspace(f_low,f_high,k+1)
        
    indexes = [list(frequency_grid).index(min(frequency_grid, key=lambda x:np.abs(x-node))) for node in desired_frequency_nodes]
    frequency_nodes = np.array(frequency_grid[indexes])

    prior['dphi_0'] = bilby.core.prior.DeltaFunction(name='dphi_0',latex_label=r'$\varphi_0$',peak=0)
    for i in list(range(len(frequency_nodes)))[1:]:
        prior[f'dphi_{i}'] = bilby.core.prior.Gaussian(name=f'dphi_{i}',latex_label=r'$\varphi_num$'.replace('num',str(i)),
                                                     mu=mean_phase_difference[indexes[i]],sigma=phase_uncertainty[indexes[i]])
    
    return frequency_nodes, prior



def dA_prior(amplitude_uncertainty,k, **kwargs):
    '''
    Generates a Gaussian prior for the amplitude correction parameters
    
    Parameters
    ===================
    amplitude_uncertainty: numpy.ndarray
        array of standard deviation of a set of amplitude differences; by default, this should be as a function of dimensionless frequency, xi
    k: int
        number of amplitude correction parameters desired
    mean_amplitude_difference: numpy.ndarray, optional
        array of the means of a set of amplitude differences, by default, this should be as a function of dimensionless frequency, xi
        default: None
    prior: bilby.core.prior.PriorDict, optional
        bilby prior object; if given, dA priors will be added to this dictionary
        default: None
    geometrized: bool, optional
        if True, will return geometrized frequency nodes; if False, normal frequency nodes (Hz)
        default: True
    xi_low: float, optional
        if geometrized is True; lower bound on the geometrized frequency band
        default: 0.018
    xi_high: float, optional
        if geometrized is True; upper bound on the geometrized frequency band
        default: 1/pi (0.318...)
    f_low: float, optional
        if geometrized is False; lower bound on the normal frequency band
        default: 20.0 Hz
    f_high: float, optional
        if geometrized is False; upper bound on the normal frequency band
        default: 1024.0 Hz
        
    Returns
    ==================
    frequency_nodes: numpy.ndarray
        array of frequency nodes
    prior: bilby.core.prior.PriorDict
        bilby prior object containing the amplitude correction priors
    '''
    f_low = kwargs.get('f_low',20)
    f_high = kwargs.get('f_high',1024)
    xi_low = kwargs.get('xi_low',0.018)
    xi_high = kwargs.get('xi_high',1/np.pi)
    prior = kwargs.get('prior',None)
    geometrized = kwargs.get('geometrized',True)
    mean_amplitude_difference = kwargs.get('mean_amplitude_difference',None)
    
    if prior is None:
        prior = bilby.core.prior.PriorDict()
    
    if mean_amplitude_difference is None:
        mean_amplitude_difference = np.array([0]*len(amplitude_uncertainty))
    
    if geometrized is True:
        frequency_grid = np.linspace(0.001,1,len(amplitude_uncertainty))
        desired_frequency_nodes = np.geomspace(xi_low,xi_high,k+1)
    else:
        frequency_grid = np.linspace(f_low,f_high,len(amplitude_uncertainty))
        desired_frequency_nodes = np.geomspace(f_low,f_high,k+1)
        
    indexes = [list(frequency_grid).index(min(frequency_grid, key=lambda x:np.abs(x-node))) for node in desired_frequency_nodes]
    frequency_nodes = np.array(frequency_grid[indexes])

    prior['dA_0'] = bilby.core.prior.DeltaFunction(name='dA_0',latex_label=r'$\alpha_0$',peak=0)
    for i in list(range(len(frequency_nodes)))[1:]:
        prior[f'dA_{i}'] = bilby.core.prior.Gaussian(name=f'dA_{i}',latex_label=r'$\alpha_num$'.replace('num',str(i)),
                                                     mu=mean_amplitude_difference[indexes[i]],sigma=amplitude_uncertainty[indexes[i]])
    
    return frequency_nodes, prior



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



def TotalMassConstraint(*,name,f_low,f_high,**kwargs):
    '''
    Generates a bilby prior to constrain the total mass
    
    Parameters
    ===================
    name: string
        name of the prior
    f_low: float
        lower bound on the frequency band
    f_high: float
        upper bound of the frequency band
    unit: string, optional
        unit of the parameter
        default: r'$\mathrm{M}_\odot$' (solar mass)
    latex_label: string, optional
        label for the parameter in LaTeX
        default: r'$M$'
    xi_low: float, optional
        lower bound on the waveform uncertainty correction in dimensionless frequency
        default: 0.018
    xi_high: float, optional
        upper bound on the waveform uncertainty correction in dimensionless frequency
        default: 1/pi

    Returns
    ==================
    total_mass_prior: bilby.core.prior.base.Constraint
        bilby constraint prior object for the total mass
    '''
    unit = kwargs.get('unit',r'$\mathrm{M}_\odot$')
    latex_label = kwargs.get('latex_label',r'$M$')
    boundary = kwargs.get('boundary',None)
    xi_low = kwargs.get('xi_low',0.018)
    xi_high = kwargs.get('xi_high',1/np.pi)
    
    total_mass_prior = bilby.core.prior.Constraint(name=name,latex_label=latex_label,minimum=xi_high*203025.4467280836/f_high, maximum=xi_low*203025.4467280836/f_low, unit=unit)
    
    return total_mass_prior



def total_mass_conversion(parameters):
    '''
    Conversion function to generate the total mass from a set of parameters; to be used alongside the total mass prior
    
    Parameters
    ==================
    parameters: dict
        dictionary of binary black hole parameters
    
    Returns
    ==================
    parameters: dict
        input parameters, but with the total mass added
    '''
    parameters['total_mass'] = bilby.gw.conversion.generate_mass_parameters(parameters)['total_mass']
    
    return parameters



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
                 frequency_nodes=None,correct_amplitude=False,correct_phase=True,
                 geometrized=True,parameter_conversion=None,
                 waveform_arguments=None):
        self._times_and_frequencies = CoupledTimeAndFrequencySeries(duration=duration,
                                                                    sampling_frequency=sampling_frequency,
                                                                    start_time=start_time)
        self.frequency_domain_source_model = frequency_domain_source_model
        self.time_domain_source_model = time_domain_source_model
        self.source_parameter_keys = self.__parameters_from_source_model()
        self.geometrized=geometrized
        self.correct_amplitude=correct_amplitude
        self.correct_phase=correct_phase
        
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
                                         'geometrized={}, ' \
                                         'correct_amplitude={}, ' \
                                         'correct_phase={}, ' \
                                         'waveform_arguments={})'\
            .format(self.duration, self.sampling_frequency, self.start_time, fdsm_name, tdsm_name,
                    param_conv_name, self.frequency_nodes, self.geometrized, self.correct_amplitude, self.correct_phase, self.waveform_arguments)
    
    def frequency_domain_strain(self, parameters=None):
        return self._calculate_strain(model=self.frequency_domain_source_model,
                                      model_data_points=self.frequency_array,
                                      parameters=parameters,
                                      transformation_function=utils.nfft,
                                      transformed_model=self.time_domain_source_model,
                                      transformed_model_data_points=self.time_array,
                                      frequency_nodes=self.frequency_nodes,
                                      geometrized=self.geometrized,
                                      correct_amplitude=self.correct_amplitude,
                                      correct_phase=self.correct_phase,
                                      )

    def time_domain_strain(self,parameters=None):
        fd_model_strain = self._calculate_strain(model=self.frequency_domain_source_model,
                                      model_data_points=self.frequency_array,
                                      parameters=parameters,
                                      transformation_function=utils.nfft,
                                      transformed_model=self.time_domain_source_model,
                                      transformed_model_data_points=self.time_array,
                                      frequency_nodes=self.frequency_nodes,
                                      geometrized=self.geometrized,
                                      correct_amplitude=self.correct_amplitude,
                                      correct_phase=self.correct_phase,
                                      )

        model_strain = dict()
        model_strain['plus'] = td_waveform(fd_model_strain['plus'],self.sampling_frequency)
        model_strain['cross'] = td_waveform(fd_model_strain['cross'],self.sampling_frequency)

        return model_strain
    
    def _calculate_strain(self, model, model_data_points, transformation_function, transformed_model,
                          transformed_model_data_points, parameters, frequency_nodes, geometrized, correct_amplitude, correct_phase):
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
        
        '''
        The following block performs the waveform uncertainty correction:
        '''
        
        if self.frequency_nodes is not None:
            
            indexes = np.arange(0,len(self.frequency_nodes),1)
            
            if self.geometrized is True:
                M = bilby.gw.conversion.generate_mass_parameters(parameters)['total_mass']
                frequency_nodes = np.array(self.frequency_nodes)*float(203025.4467280836/M)
            else:
                frequency_nodes = self.frequency_nodes

            temp_frequency_grid = np.linspace(frequency_nodes[0],frequency_nodes[-1],1000)

            if self.correct_amplitude is True:
                try:
                    alphas = [parameters[f'dA_{i}'] for i in indexes]
                    dA_spline = scipy.interpolate.CubicSpline(frequency_nodes,alphas)(temp_frequency_grid)
                    dA = np.interp(self.frequency_array,temp_frequency_grid,dA_spline)
                except:
                    raise Exception('Amplitude Correction Failed!')
            else:
                dA = 0

            if self.correct_phase is True:
                try:
                    phis = [parameters[f'dphi_{i}'] for i in indexes]
                    dphi_spline = scipy.interpolate.CubicSpline(frequency_nodes,phis)(temp_frequency_grid)
                    dphi = np.interp(self.frequency_array,temp_frequency_grid,dphi_spline)
                except:
                    raise Exception('Phase Correction Failed!')
            else:
                dphi = 0
                    
        else:
            dA = 0
            dphi = 0
            
        model_strain['plus'] = model_strain['plus']*(1+dA)*np.exp(dphi*1j)
        model_strain['cross'] = model_strain['cross']*(1+dA)*np.exp(dphi*1j)
        
        self._cache['waveform'] = model_strain
        self._cache['parameters'] = self.parameters.copy()
        self._cache['model'] = model
        self._cache['transformed_model'] = transformed_model
        self._cache['amplitude_difference'] = dA
        self._cache['phase_difference'] = dphi
                              
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
                indexes = np.arange(0,len(self.frequency_nodes),1)
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



#######################################
##Temporary Utility Function Location##
#######################################

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

def td_waveform(fd_waveform,sampling_frequency):
    reversed_fd_waveform = fd_waveform[::-1]
    total_fd_strain = np.concatenate((fd_waveform,np.conjugate(reversed_fd_waveform[1:-1])))
    td_waveform = sampling_frequency*np.real(np.fft.ifft(total_fd_strain))
    return td_waveform
