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
from .utils import progressBar

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
