from pycbc.waveform import td_approximants, fd_approximants
import numpy as np
import matplotlib.pyplot as plt
import pesummary
from pesummary.gw.file.strain import StrainData
from pesummary.io import read
import subprocess
import lal
import random
from pesummary.utils.array import Array
from pesummary.utils.samples_dict import MCMCSamplesDict
from pesummary.utils.samples_dict import MultiAnalysisSamplesDict
import time as tm
import sys
import bilby
from bilby.core.prior import Uniform, PowerLaw, Cosine, Sine, Constraint, Gaussian, DeltaFunction
from gwpy.timeseries import TimeSeries
import pandas as pd
import IPython as ipy
import scipy
from bilby.core import utils
from bilby.core.series import CoupledTimeAndFrequencySeries
from bilby.core.utils import PropertyAccessor
from bilby.gw.conversion import convert_to_lal_binary_neutron_star_parameters


class WaveformGeneratorWFU(object):
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
    
    
trigger_time = 1187008882.4

H1 = bilby.gw.detector.get_empty_interferometer("H1")
L1 = bilby.gw.detector.get_empty_interferometer("L1")

# Definite times in relation to the trigger time (time_of_event), duration and post_trigger_duration
post_trigger_duration = 2
duration = 100
analysis_start = trigger_time + post_trigger_duration - duration

sampling_frequency = 4096

# Use gwpy to fetch the open data
H1_analysis_data = TimeSeries.fetch_open_data(
    "H1", analysis_start, analysis_start + duration).resample(sampling_frequency)

L1_analysis_data = TimeSeries.fetch_open_data(
    "L1", analysis_start, analysis_start + duration).resample(sampling_frequency)
H1.set_strain_data_from_gwpy_timeseries(H1_analysis_data)
L1.set_strain_data_from_gwpy_timeseries(L1_analysis_data)
psd_duration = duration * 32 # 32 noise realizations, effectively
psd_start_time = analysis_start - psd_duration

H1_psd_data = TimeSeries.fetch_open_data(
    "H1", psd_start_time, psd_start_time + psd_duration).resample(sampling_frequency)

L1_psd_data = TimeSeries.fetch_open_data(
    "L1", psd_start_time, psd_start_time + psd_duration).resample(sampling_frequency)
psd_alpha = 2 * H1.strain_data.roll_off / duration # controlling windowing functions for the psd
H1_psd = H1_psd_data.psd(fftlength=duration, overlap=0, window=("tukey", psd_alpha), method="median")
L1_psd = L1_psd_data.psd(fftlength=duration, overlap=0, window=("tukey", psd_alpha), method="median")
H1.power_spectral_density = bilby.gw.detector.PowerSpectralDensity(
    frequency_array=H1_psd.frequencies.value, psd_array=H1_psd.value)
L1.power_spectral_density = bilby.gw.detector.PowerSpectralDensity(
    frequency_array=H1_psd.frequencies.value, psd_array=L1_psd.value)
f_max = 2048
H1.maximum_frequency = f_max
L1.maximum_frequency = f_max

interferometers = [H1, L1]

prior = bilby.core.prior.PriorDict()
trigger_time = 1187008882.4

prior['mass_1'] = Uniform(name=r'$m_{1}$',minimum=1.37563293,maximum=1.80367393)
prior['mass_2'] = Uniform(name=r'$m_{2}$',minimum=1.06371259,maximum=1.3758219)
prior['a_1'] = Uniform(name=r'$a_{1}$', minimum=0, maximum=0.05)
prior['a_2'] = Uniform(name=r'$a_{2}$', minimum=0, maximum=0.05)
prior['luminosity_distance'] = bilby.gw.prior.UniformSourceFrame(name='luminosity_distance', minimum=1, maximum=500.0, unit='Mpc')
prior['phi_12'] = Uniform(name=r'$\Phi_{12}$', minimum=0, maximum=6.2831853071795865,boundary='periodic')
prior['phi_jl'] = Uniform(name=r'$\Phi_{JL}$', minimum=0, maximum=6.2831853071795865,boundary='periodic')
prior['tilt_1'] = Uniform(name=r'$\theta_{1}$', minimum=0, maximum=6.2831853071795865, boundary='periodic')
prior['tilt_2'] = Uniform(name=r'$\theta_{2}$', minimum=0, maximum=6.2831853071795865, boundary='periodic')
prior['theta_jn'] = Sine(name=r'$\theta_{JN}$')
prior['phase'] = Uniform(name=r'$\Phi$', minimum=0, maximum=6.2831853071795865, boundary='periodic')
prior['lambda_1'] = Uniform(name=r"$\Lambda_{1}$", minimum=0.00147326, maximum=3154.41685213)
prior['lambda_2'] = Uniform(name=r"$\Lambda_{2}$", minimum=0.02966776, maximum=4598.76616739)
prior['geocent_time'] = Uniform(name=r"$t_{c}$", minimum=trigger_time-0.1, maximum=trigger_time+0.1)

prior['dec'] = DeltaFunction(name=r'$\delta$', peak=-0.408084)
prior['ra'] = DeltaFunction(name=r'\alpha_{r}$', peak=3.44616)
prior['psi'] = DeltaFunction(name=r'$\Psi$', peak=1.56379256)

'''
prior['alpha_1'] = Gaussian(name=r'$\alpha_{1}$', mu=-0.00018402606398450926, sigma=0.00041816788658705296)
prior['alpha_2'] = Gaussian(name=r'$\alpha_{2}$', mu=-0.00017364270287341895, sigma=0.0004034522758877465)
prior['alpha_3'] = Gaussian(name=r'$\alpha_{3}$', mu=-0.00014348552399630194, sigma=0.0003606729292871294)
prior['alpha_4'] = Gaussian(name=r'$\alpha_{4}$', mu=-0.00006157481730980426, sigma=0.0002442709550179268)
prior['alpha_5'] = Gaussian(name=r'$\alpha_{5}$', mu=0.00011969536115612494, sigma=0.00004047199756191817)
prior['alpha_6'] = Gaussian(name=r'$\alpha_{6}$', mu=0.00026954428003384435, sigma=0.0002408182715336409)
prior['alpha_7'] = Gaussian(name=r'$\alpha_{7}$', mu=0.0001931860309102457, sigma=0.0002762592144232504)
prior['alpha_8'] = Gaussian(name=r'$\alpha_{8}$', mu=0.0000018591905683778544, sigma=0.00047407638899255445)
prior['alpha_9'] = Gaussian(name=r'$\alpha_{9}$', mu=0.02609684856527983, sigma=0.024743687703602003)
prior['alpha_10'] = Gaussian(name=r'$\alpha_{10}$', mu=0.14220478498845276, sigma=0.14400944909014002)
'''

prior['beta_1'] = Gaussian(name=r'$\beta_{1}$', mu=0.02323465309328362, sigma=0.001719726076207944)
prior['beta_2'] = Gaussian(name=r'$\beta_{2}$', mu=0.022953500183714576, sigma=0.001687564544274498)
prior['beta_3'] = Gaussian(name=r'$\beta_{3}$', mu=0.02212099262556966, sigma=0.0015955662769974318)
prior['beta_4'] = Gaussian(name=r'$\beta_{4}$', mu=0.019720233646154717, sigma=0.0013586610184014127)
prior['beta_5'] = Gaussian(name=r'$\beta_{5}$', mu=0.01322755078369319, sigma=0.0009456734431152912)
prior['beta_6'] = Gaussian(name=r'$\beta_{6}$', mu=0.000021649200314916407, sigma=0.0006112668713038633)
prior['beta_7'] = Gaussian(name=r'$\beta_{7}$', mu=-0.0075365189575646195, sigma=0.001013474721440995)
prior['beta_8'] = Gaussian(name=r'$\beta_{8}$', mu=0.036520070930041254, sigma=0.0025976631174632887)
prior['beta_9'] = Gaussian(name=r'$\beta_{9}$', mu=0.19783526758528153, sigma=0.0966458034438057)
prior['beta_10'] = Gaussian(name=r'$\beta_{10}$', mu=0.987196509934645, sigma=0.6498396192690529)

frequency_nodes = np.arange(20,2048.1,0.1)[np.geomspace(1,20280,10).astype(int)]

ht1 = WaveformGeneratorWFU(
                    parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_neutron_star_parameters,
                    waveform_arguments=dict(waveform_approximant='IMRPhenomPv2_NRTidalv2',reference_frequency=50,catch_waveform_errors=True,f_low=20.0, f_high=2048.0),
                    frequency_domain_source_model=bilby.gw.source.lal_binary_neutron_star, 
                    sampling_frequency=4096, 
                    duration=256,
                    waveform_uncertainty_nodes=frequency_nodes,dA_sampling=False,dphi_sampling=True,
                )

likelihood = bilby.gw.likelihood.GravitationalWaveTransient(
    interferometers,
    ht1,
    priors=prior,
    time_marginalization=False, 
    phase_marginalization=False, 
    distance_marginalization=False,
)

result_short = bilby.run_sampler(
    likelihood, prior, sampler='nestle', outdir='/home/ryanmatthew.johnson/Waveform_Uncertainty/', label="GW170817-WFUtest",
    conversion_function=bilby.gw.conversion.generate_all_bns_parameters,
    nlive=5, 
    dlogz=1.,
    clean=True,
    maxiter=50,
)