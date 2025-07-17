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
from .utils import smooth_interpolation, td_waveform, variable_prior, GC_waveform_correction



def GeneralCorrectionModelBBH(
        frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1,
        phi_12, a_2, tilt_2, phi_jl, theta_jn, phase,
        xi_0, delta_xi_tilde, dAs, dphis, **kwargs):
    """ A Binary Black Hole waveform model using lalsimulation

    Parameters
    ==========
    frequency_array: array_like
        The frequencies at which we want to calculate the strain
    mass_1: float
        The mass of the heavier object in solar masses
    mass_2: float
        The mass of the lighter object in solar masses
    luminosity_distance: float
        The luminosity distance in megaparsec
    a_1: float
        Dimensionless primary spin magnitude
    tilt_1: float
        Primary tilt angle
    phi_12: float
        Azimuthal angle between the two component spins
    a_2: float
        Dimensionless secondary spin magnitude
    tilt_2: float
        Secondary tilt angle
    phi_jl: float
        Azimuthal angle between the total binary angular momentum and the
        orbital angular momentum
    theta_jn: float
        Angle between the total binary angular momentum and the line of sight
    phase: float
        The phase at reference frequency or peak amplitude (depends on waveform)
    kwargs: dict
        Optional keyword arguments
        Supported arguments:

        - waveform_approximant
        - reference_frequency
        - minimum_frequency
        - maximum_frequency
        - catch_waveform_errors
        - pn_spin_order
        - pn_tidal_order
        - pn_phase_order
        - pn_amplitude_order
        - mode_array:
          Activate a specific mode array and evaluate the model using those
          modes only.  e.g. waveform_arguments =
          dict(waveform_approximant='IMRPhenomHM', mode_array=[[2,2],[2,-2])
          returns the 22 and 2-2 modes only of IMRPhenomHM.  You can only
          specify modes that are included in that particular model.  e.g.
          waveform_arguments = dict(waveform_approximant='IMRPhenomHM',
          mode_array=[[2,2],[2,-2],[5,5],[5,-5]]) is not allowed because the
          55 modes are not included in this model.  Be aware that some models
          only take positive modes and return the positive and the negative
          mode together, while others need to call both.  e.g.
          waveform_arguments = dict(waveform_approximant='IMRPhenomHM',
          mode_array=[[2,2],[4,-4]]) returns the 22 and 2-2 of IMRPhenomHM.
          However, waveform_arguments =
          dict(waveform_approximant='IMRPhenomXHM', mode_array=[[2,2],[4,-4]])
          returns the 22 and 4-4 of IMRPhenomXHM.
        - lal_waveform_dictionary:
          A dictionary (lal.Dict) of arguments passed to the lalsimulation
          waveform generator. The arguments are specific to the waveform used.

    Returns
    =======
    dict: A dictionary with the plus and cross polarisation strain modes
    """
    waveform_approximant = kwargs.get('waveform_approximant','IMRPhenomPv2')
    reference_frequency = kwargs.get('reference_frequency',50.0)
    minimum_frequency = kwargs.get('minimum_frequency',20.0)
    maximum_frequency = kwargs.get('maximum_frequency',frequency_array[-1])
    catch_waveform_errors = kwargs.get('catch_waveform_errors',False)
    pn_spin_order = kwargs.get('pn_spin_order',-1)
    pn_tidal_order = kwargs.get('pn_tidal_order',-1)
    pn_phase_order = kwargs.get('pn_phase_order',-1)
    pn_amplitude_order = kwargs.get('pn_amplitude_order',0)
    sigma_dA_spline = kwargs.get('sigma_dA_spline',None)
    sigma_dphi_spline = kwargs.get('sigma_dphi_spline',None)
    xi_high = kwargs.get('xi_high',1/np.pi)
    gamma = kwargs.get('gamma',0.025)
    
    waveform_arguments = dict(
        waveform_approximant=waveform_approximant, reference_frequency=reference_frequency,
        minimum_frequency=minimum_frequency, maximum_frequency=maximum_frequency,
        catch_waveform_errors=catch_waveform_errors, pn_spin_order=pn_spin_order, pn_tidal_order=pn_tidal_order,
        pn_phase_order=pn_phase_order, pn_amplitude_order=pn_amplitude_order,
    )
    
    model_strain = bilby.gw.source._base_lal_cbc_fd_waveform(
        frequency_array=frequency_array, mass_1=mass_1, mass_2=mass_2,
        luminosity_distance=luminosity_distance, theta_jn=theta_jn, phase=phase,
        a_1=a_1, a_2=a_2, tilt_1=tilt_1, tilt_2=tilt_2, phi_12=phi_12,
        phi_jl=phi_jl, **waveform_arguments)
    
    waveform_correction = GC_waveform_correction(frequency_array,xi_0,delta_xi_tilde,dAs,dphis,
                                                 sigma_dA_spline,sigma_dphi_spline,
                                                 mass_1,mass_2,xi_high,gamma)
    
    model_strain['plus'] *= waveform_correction
    model_strain['cross'] *= waveform_correction
    
    return model_strain



def binary_black_hole_correction_conversion(parameters):
    """
    Convert parameters we have into parameters we need.

    This is defined by the parameters of bilby.source.lal_binary_black_hole()


    Mass: mass_1, mass_2
    Spin: a_1, a_2, tilt_1, tilt_2, phi_12, phi_jl
    Extrinsic: luminosity_distance, theta_jn, phase, ra, dec, geocent_time, psi

    This involves popping a lot of things from parameters.
    The keys in added_keys should be popped after evaluating the waveform.

    Parameters
    ==========
    parameters: dict
        dictionary of parameter values to convert into the required parameters

    Returns
    =======
    converted_parameters: dict
        dict of the required parameters
    added_keys: list
        keys which are added to parameters during function call
    """

    converted_parameters = parameters.copy()
    original_keys = list(converted_parameters.keys())
    if 'luminosity_distance' not in original_keys:
        if 'redshift' in converted_parameters.keys():
            converted_parameters['luminosity_distance'] = \
                redshift_to_luminosity_distance(parameters['redshift'])
        elif 'comoving_distance' in converted_parameters.keys():
            converted_parameters['luminosity_distance'] = \
                comoving_distance_to_luminosity_distance(
                    parameters['comoving_distance'])

    for key in original_keys:
        if key[-7:] == '_source':
            if 'redshift' not in converted_parameters.keys():
                converted_parameters['redshift'] =\
                    luminosity_distance_to_redshift(
                        parameters['luminosity_distance'])
            converted_parameters[key[:-7]] = converted_parameters[key] * (
                1 + converted_parameters['redshift'])

    # we do not require the component masses be added if no mass parameters are present
    converted_parameters = bilby.gw.conversion.generate_component_masses(converted_parameters, require_add=False)

    for idx in ['1', '2']:
        key = 'chi_{}'.format(idx)
        if key in original_keys:
            if "chi_{}_in_plane".format(idx) in original_keys:
                converted_parameters["a_{}".format(idx)] = (
                    converted_parameters[f"chi_{idx}"] ** 2
                    + converted_parameters[f"chi_{idx}_in_plane"] ** 2
                ) ** 0.5
                converted_parameters[f"cos_tilt_{idx}"] = (
                    converted_parameters[f"chi_{idx}"]
                    / converted_parameters[f"a_{idx}"]
                )
            elif "a_{}".format(idx) not in original_keys:
                converted_parameters['a_{}'.format(idx)] = abs(
                    converted_parameters[key])
                converted_parameters['cos_tilt_{}'.format(idx)] = \
                    np.sign(converted_parameters[key])
            else:
                with np.errstate(invalid="raise"):
                    try:
                        converted_parameters[f"cos_tilt_{idx}"] = (
                            converted_parameters[key] / converted_parameters[f"a_{idx}"]
                        )
                    except (FloatingPointError, ZeroDivisionError):
                        logger.debug(
                            "Error in conversion to spherical spin tilt. "
                            "This is often due to the spin parameters being zero. "
                            f"Setting cos_tilt_{idx} = 1."
                        )
                        converted_parameters[f"cos_tilt_{idx}"] = 1.0

    for key in ["phi_jl", "phi_12"]:
        if key not in converted_parameters:
            converted_parameters[key] = 0.0

    for angle in ['tilt_1', 'tilt_2', 'theta_jn']:
        cos_angle = str('cos_' + angle)
        if cos_angle in converted_parameters.keys():
            with np.errstate(invalid="ignore"):
                converted_parameters[angle] = np.arccos(converted_parameters[cos_angle])

    if "delta_phase" in original_keys:
        with np.errstate(invalid="ignore"):
            converted_parameters["phase"] = np.mod(
                converted_parameters["delta_phase"]
                - np.sign(np.cos(converted_parameters["theta_jn"]))
                * converted_parameters["psi"],
                2 * np.pi)
    
    dA_keys = [key for key in parameters if 'dA_' in key]
    dphi_keys = [key for key in parameters if 'dphi_' in key]
    
    if len(dA_keys) != 0 and len(dphi_keys) != 0:
        raise Exception('Amplitude and Phase corrections do not have the same number of parameters!')
    elif len(dA_keys) == 0 and len(dphi_keys) == 0:
        raise Exception('No waveform correction detected!')
    
    if len(dA_keys) != 0:
        dAs = [parameters[key] for key in dA_keys]
        dphis = list(np.zeros(len(dA_keys)))
    if len(dphi_keys) != 0:
        dAs = list(np.zeros(len(dphi_keys)))
        dphis = [parameters[key] for key in dphi_keys]
    
    converted_parameters['dAs'] = dAs
    converted_parameters['dphis'] = dphis
    
    added_keys = [key for key in converted_parameters.keys()
                  if key not in original_keys]

    return converted_parameters, added_keys



class BasicCorrectionModel(object):
    '''
    Modified WaveformGenerator object from bilby.gw to include waveform uncertainty corrections in the strain calculation
    To sample waveform uncertainty, include all relevant "alpha" and "phi" parameters in the prior.
    Note: make sure the number of alphas, phis, and waveform_uncertainty_nodes are the same
    
    New Parameters
    ==================
   frequency_nodes: numpy.ndarray, optional
        array of frequency nodes to be used in generating the dA and dphi splines
        default: None
    correct_amplitude: bool, optional
        if True, the waveform generator will attempt to pull dA parameters from the parameter dictionary (either an injection or the prior)
        default: None
    correct_phase: bool, optional
        if True, the waveform generator will attempt to pull dphi parameters from the parameter dictionary (either an injection or the prior)
        default: None
    dimensionless: bool, optional
        whether or not to perfom the waveform correction in dimensionless frequency or not
        default: True
    '''
    duration = PropertyAccessor('_times_and_frequencies', 'duration')
    sampling_frequency = PropertyAccessor('_times_and_frequencies', 'sampling_frequency')
    start_time = PropertyAccessor('_times_and_frequencies', 'start_time')
    frequency_array = PropertyAccessor('_times_and_frequencies', 'frequency_array')
    time_array = PropertyAccessor('_times_and_frequencies', 'time_array')
    def __init__(self, duration=None, sampling_frequency=None, start_time=0, frequency_domain_source_model=None,
                 time_domain_source_model=None, parameters=None,
                 frequency_nodes=None,correct_amplitude=False,correct_phase=False,
                 dimensionless=True,parameter_conversion=None,
                 waveform_arguments=None):
        self._times_and_frequencies = CoupledTimeAndFrequencySeries(duration=duration,
                                                                    sampling_frequency=sampling_frequency,
                                                                    start_time=start_time)
        self.frequency_domain_source_model = frequency_domain_source_model
        self.time_domain_source_model = time_domain_source_model
        self.source_parameter_keys = self.__parameters_from_source_model()
        self.dimensionless=dimensionless
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
                                         'dimensionless={}, ' \
                                         'correct_amplitude={}, ' \
                                         'correct_phase={}, ' \
                                         'waveform_arguments={})'\
            .format(self.duration, self.sampling_frequency, self.start_time, fdsm_name, tdsm_name,
                    param_conv_name, self.frequency_nodes, self.dimensionless, self.correct_amplitude, self.correct_phase, self.waveform_arguments)
    
    def frequency_domain_strain(self, parameters=None):
        return self._calculate_strain(model=self.frequency_domain_source_model,
                                      model_data_points=self.frequency_array,
                                      parameters=parameters,
                                      transformation_function=utils.nfft,
                                      transformed_model=self.time_domain_source_model,
                                      transformed_model_data_points=self.time_array,
                                      frequency_nodes=self.frequency_nodes,
                                      dimensionless=self.dimensionless,
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
                                      dimensionless=self.dimensionless,
                                      correct_amplitude=self.correct_amplitude,
                                      correct_phase=self.correct_phase,
                                      )

        model_strain = dict()
        model_strain['plus'] = td_waveform(fd_model_strain['plus'],self.sampling_frequency)
        model_strain['cross'] = td_waveform(fd_model_strain['cross'],self.sampling_frequency)

        return model_strain
    
    def _calculate_strain(self, model, model_data_points, transformation_function, transformed_model,
                          transformed_model_data_points, parameters, frequency_nodes, dimensionless, correct_amplitude, correct_phase):
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
            
            if self.dimensionless is True:
                M = bilby.gw.conversion.generate_mass_parameters(parameters)['total_mass']
                frequency_nodes = np.array(self.frequency_nodes)*float(203025.4467280836/M)
            else:
                frequency_nodes = self.frequency_nodes

            try:
                gamma = parameters['gamma']
            except:
                gamma = 0.025
            
            if self.correct_amplitude is True:
                try:
                    alphas = [parameters[f'dA_{i}'] for i in indexes]
                    dA = smooth_interpolation(self.frequency_array,frequency_nodes,alphas,gamma)
                except:
                    raise Exception('Amplitude Correction Failed!')
            else:
                dA = 0

            if self.correct_phase is True:
                try:
                    phis = [parameters[f'dphi_{i}'] for i in indexes]
                    dphi = smooth_interpolation(self.frequency_array,frequency_nodes,phis,gamma)
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
            if transformation_function == bilby.utils.nfft:
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
                if key not in [f'dA_{i}' for i in indexes]+[f'dphi_{i}' for i in indexes]+['gamma_dphi','gamma_dA']:  
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



class GeneralCorrectionModel(object):
    '''
    Modified WaveformGenerator object from bilby.gw to include waveform uncertainty corrections in the strain calculation
    
    New Parameters
    ==================
    correction_arguments: dict, optional
        dictionary containing arguments for the waveform correction
        default: None

        contents:
            correct_amplitude: bool
                whether or not to attempt an amplitude correction
                include all dA parameters in the prior if True
            sigma_dA: numpy.ndarry
                if correct_amplitude is True, this is the array of standard deviations for the dA priors
            correct_phase: bool
                whether or not to attempt a phase correction
                include all dphi parameters in the prior if True
            sigma_dphi: numpy.ndarry
                if correct_phase is True, this is the array of standard deviations for the dphi priors
            nodes: int
                number of frequency nodes desired
            xi_high: float, optional
                absolute lower bound on dimensionless frequency, xi
                default: 1/pi
    '''
    duration = PropertyAccessor('_times_and_frequencies', 'duration')
    sampling_frequency = PropertyAccessor('_times_and_frequencies', 'sampling_frequency')
    start_time = PropertyAccessor('_times_and_frequencies', 'start_time')
    frequency_array = PropertyAccessor('_times_and_frequencies', 'frequency_array')
    time_array = PropertyAccessor('_times_and_frequencies', 'time_array')
    def __init__(self, duration=None, sampling_frequency=None, start_time=0, frequency_domain_source_model=None,
                 time_domain_source_model=None, parameters=None,
                 correction_arguments=None,parameter_conversion=None,
                 waveform_arguments=None):
        self._times_and_frequencies = CoupledTimeAndFrequencySeries(duration=duration,
                                                                    sampling_frequency=sampling_frequency,
                                                                    start_time=start_time)
        self.frequency_domain_source_model = frequency_domain_source_model
        self.time_domain_source_model = time_domain_source_model
        self.source_parameter_keys = self.__parameters_from_source_model()
        self.correction_arguments=correction_arguments
        
        if parameter_conversion is None:
            self.parameter_conversion = convert_to_lal_binary_black_hole_parameters
        else:
            self.parameter_conversion = parameter_conversion
        if waveform_arguments is not None:
            self.waveform_arguments = waveform_arguments
        else:
            self.waveform_arguments = dict()
        if self.correction_arguments is not None:
            self.correction_arguments = correction_arguments
        else:
            self.correction_arguments = None
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
                                         'correction_arguments={}, ' \
                                         'waveform_arguments={})'\
            .format(self.duration, self.sampling_frequency, self.start_time, fdsm_name, tdsm_name,
                    param_conv_name, self.correction_arguments, self.waveform_arguments)
    
    def frequency_domain_strain(self, parameters=None):
        return self._calculate_strain(model=self.frequency_domain_source_model,
                                      model_data_points=self.frequency_array,
                                      parameters=parameters,
                                      transformation_function=utils.nfft,
                                      transformed_model=self.time_domain_source_model,
                                      transformed_model_data_points=self.time_array,
                                      correction_arguments=self.correction_arguments,
                                      )

    def time_domain_strain(self,parameters=None):
        fd_model_strain = self._calculate_strain(model=self.frequency_domain_source_model,
                                      model_data_points=self.frequency_array,
                                      parameters=parameters,
                                      transformation_function=utils.nfft,
                                      transformed_model=self.time_domain_source_model,
                                      transformed_model_data_points=self.time_array,
                                      correction_arguments=self.correction_arguments,
                                      )

        model_strain = dict()
        model_strain['plus'] = td_waveform(fd_model_strain['plus'],self.sampling_frequency)
        model_strain['cross'] = td_waveform(fd_model_strain['cross'],self.sampling_frequency)

        return model_strain
    
    def _calculate_strain(self, model, model_data_points, transformation_function, transformed_model,
                          transformed_model_data_points, parameters, correction_arguments):
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
                              
        try:
            gamma = parameters['gamma']
        except:
            gamma = 0.025

        try:
            xi_high = correction_arguments['xi_high']
        except:
            xi_high = 1/np.pi

        M = bilby.gw.conversion.generate_mass_parameters(parameters)['total_mass']
        indexes = np.arange(0,correction_arguments['nodes']+1,1)
                              
        if correction_arguments['correct_amplitude'] is True:
            sigma_dA = correction_arguments['sigma_dA']
            xi_up = parameters['xi_0']*(1+((xi_high-parameters['xi_0'])/(parameters['xi_0']))*parameters['delta_xi_tilde'])
            dA_frequency_nodes,dA_coeffs = variable_prior(sigma_dA,correction_arguments['nodes'],
                                                      parameters['xi_0'], xi_up)
            dA_frequency_nodes *= float(203025.4467280836/M)
            try:
                prior_alphas = np.array([parameters[f'dA_{i}'] for i in indexes])
                alphas = prior_alphas*dA_coeffs
                dA = smooth_interpolation(self.frequency_array,dA_frequency_nodes,alphas,gamma)
            except:
                raise Exception('Amplitude Correction Failed!')
        else:
            dA = 0
        if correction_arguments['correct_phase'] is True:
            sigma_dphi = correction_arguments['sigma_dphi']
            xi_up = parameters['xi_0']*(1+((xi_high-parameters['xi_0'])/(parameters['xi_0']))*parameters['delta_xi_tilde'])
            dphi_frequency_nodes,dphi_coeffs = variable_prior(sigma_dphi,correction_arguments['nodes'],
                                                            parameters['xi_0'], xi_up)
            dphi_frequency_nodes *= float(203025.4467280836/M)
            try:
                prior_phis = np.array([parameters[f'dphi_{i}'] for i in indexes])
                phis = prior_phis*dphi_coeffs
                dphi = smooth_interpolation(self.frequency_array,dphi_frequency_nodes,phis,gamma)
            except:
                raise Exception('Phase Correction Failed!')
        else:
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
            indexes = np.arange(0,self.correction_arguments['nodes']+1,1)
            if key not in [f'dA_{i}' for i in indexes]+[f'dphi_{i}' for i in indexes]+['gamma','xi_0','delta_xi_tilde']:  
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
