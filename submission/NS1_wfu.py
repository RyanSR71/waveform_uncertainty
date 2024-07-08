import os
os.environ['LAL_DATA_PATH'] = '/home/carl-johan.haster/WaveForm_data/lalsuite-extra/data/lalsimulation'
import numpy as np
import bilby
import matplotlib.pyplot as plt
import random
import time as tm
import sys
import scipy
import lal
import pesummary
from pesummary.gw.file.strain import StrainData
from bilby.core.prior import Uniform, PowerLaw, Cosine, Sine, Constraint, Gaussian, DeltaFunction
from gwpy.timeseries import TimeSeries
from pesummary.io import read
import WaveformUncertainty as wfu


gw170817_posterior = read("/home/jocelyn.read/sample_reference/GW170817_bilby_pesummary.dat").samples_dict
maxL = gw170817_posterior.maxL
trigger_time = 1187008882.4

prior = bilby.core.prior.PriorDict()

prior['mass_1'] = Uniform(name='mass_1',latex_label=r'$m_{1}$',
                          minimum=float(np.min(gw170817_posterior['mass_1'])),maximum=float(np.max(gw170817_posterior['mass_1'])),unit=r'$\mathrm{M}_{\odot}$')
prior['a_1'] = Uniform(name='a_1',latex_label=r'$a_{1}$', minimum=float(np.min(gw170817_posterior['a_1'])),maximum=float(np.max(gw170817_posterior['a_1'])))
prior['tilt_1'] = Uniform(name='tilt_1',latex_label=r'$\theta_{1}$', 
                          minimum=float(np.min(gw170817_posterior['tilt_1'])),maximum=float(np.max(gw170817_posterior['tilt_1'])))
prior['lambda_1'] = Uniform(name='lambda_1',latex_label=r"$\Lambda_{1}$",minimum=float(np.min(gw170817_posterior['lambda_1'])),
                            maximum=float(np.max(gw170817_posterior['lambda_1'])))

prior['phase'] = Uniform(name='phase',latex_label=r'$\Phi$',minimum=float(np.min(gw170817_posterior['phase'])),
                         maximum=float(np.max(gw170817_posterior['phase'])))
prior['geocent_time'] = Uniform(name='geocent_time',latex_label=r"$t_{c}$",minimum=float(np.min(gw170817_posterior['geocent_time'])),
                                maximum=float(np.max(gw170817_posterior['geocent_time'])),unit='s')

prior['mass_2'] = DeltaFunction(name='mass_2',latex_label=r'$m_{2}$',peak=float(maxL['mass_2']),unit=r'$\mathrm{M}_{\odot}$')
prior['a_2'] = DeltaFunction(name='a_2',latex_label=r'$a_{2}$', peak=float(maxL['a_2']))
prior['luminosity_distance'] = DeltaFunction(name='luminosity_distance',latex_label=r'$d_{L}$',peak=float(maxL['luminosity_distance']), unit='Mpc')
prior['phi_12'] = DeltaFunction(name='phi_12',latex_label=r'$\Phi_{12}$', peak=float(maxL['phi_12']))
prior['phi_jl'] = DeltaFunction(name='phi_jl',latex_label=r'$\Phi_{JL}$', peak=float(maxL['phi_jl']))
prior['tilt_2'] = DeltaFunction(name='tilt_2',latex_label=r'$\theta_{2}$', peak=float(maxL['tilt_2']))
prior['theta_jn'] = DeltaFunction(name='theta_jn',latex_label=r'$\theta_{JN}$',peak=float(maxL['theta_jn']))
prior['lambda_2'] = DeltaFunction(name='lambda_2',latex_label=r"$\Lambda_{2}$",peak=float(maxL['lambda_2']))
prior['dec'] = DeltaFunction(name='dec',latex_label=r'$\delta$', peak=float(maxL['dec']))
prior['ra'] = DeltaFunction(name='ra',latex_label=r'\alpha_{r}$', peak=float(maxL['ra']))
prior['psi'] = DeltaFunction(name='psi',latex_label=r'$\Psi$', peak=float(maxL['psi']))

parameterization = np.load("/home/ryanmatthew.johnson/Waveform_Uncertainty/files/parameterization_nsamples_1000.npy",allow_pickle=True)

mean_amplitude_difference,amplitude_uncertainty,mean_phase_difference,phase_uncertainty,linear_frequency_grid = wfu.uncertainties_from_parameterization(parameterization,linear=True,resolution=0.1)

prior,frequency_nodes = wfu.WFU_prior(prior,mean_amplitude_difference=[0]*len(linear_frequency_grid),amplitude_uncertainty=amplitude_uncertainty,
                                      mean_phase_difference=[0]*len(linear_frequency_grid),phase_uncertainty=phase_uncertainty,spacing='geometric',
                                      frequency_grid=linear_frequency_grid,nnodes=6)


hf1 = wfu.WaveformGeneratorWFU(
                    parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_neutron_star_parameters,
                    waveform_arguments=dict(waveform_approximant='IMRPhenomPv2_NRTidalv2', reference_frequency=50, catch_waveform_errors=True, 
                                            f_low = 20.0, f_high=2048.0),
                    frequency_domain_source_model=bilby.gw.source.lal_binary_neutron_star, 
                    sampling_frequency=4096, 
                    duration=256,
                    waveform_uncertainty_nodes=frequency_nodes,dA_sampling=True,dphi_sampling=True,
                )


fixed_injection = dict()
for parameter in ['mass_1', 'mass_2', 'a_1', 'a_2', 'luminosity_distance', 'phi_12', 'phi_jl', 'tilt_1', 'tilt_2', 'theta_jn', 'phase', 'lambda_1', 'lambda_2', 'geocent_time', 'dec', 'ra', 'psi']:
    fixed_injection[parameter] = float(gw170817_posterior.maxL[parameter])
    

injected_amplitude_difference = np.interp(hf1.frequency_array,linear_frequency_grid,mean_amplitude_difference)
injected_phase_difference = np.interp(hf1.frequency_array,linear_frequency_grid,mean_phase_difference)
injected_waveform = wfu.WaveformGeneratorWFU_inj(
                    parameter_conversion=bilby.gw.conversion.convert_to_lal_binary_neutron_star_parameters,
                    waveform_arguments=dict(waveform_approximant='IMRPhenomPv2_NRTidalv2', reference_frequency=50, catch_waveform_errors=True, f_low = 20.0, f_high=2048.0),
                    frequency_domain_source_model=bilby.gw.source.lal_binary_neutron_star, 
                    sampling_frequency=4096, 
                    duration=256,
                    amplitude_difference=0.25*injected_amplitude_difference,phase_difference=0.25*injected_phase_difference,
                )


ifos = bilby.gw.detector.InterferometerList(['H1', 'L1'])
ifos.set_strain_data_from_power_spectral_densities(
    sampling_frequency=4096, duration=256,
    start_time=maxL['geocent_time'] - 3)
ifo_injection = ifos.inject_signal(
    waveform_generator=injected_waveform,
    parameters=fixed_injection)


likelihood = bilby.gw.likelihood.GravitationalWaveTransient(
    ifos,
    hf1,
    priors=prior,
    time_marginalization=True, 
    phase_marginalization=True, 
    distance_marginalization=False,
)


result = bilby.run_sampler(
    likelihood, prior, sampler='dynesty', outdir='/home/ryanmatthew.johnson/Waveform_Uncertainty/output/NS1/', 
    label="NS1_wfu_nlive_1000",
    conversion_function=bilby.gw.conversion.generate_all_bns_parameters,
    resume=True,
    nlive=1000, 
    dlogz=1.,
    clean=True,
    maxiter=None,
)