import numpy as np
import bilby
import random
import time
import sys
import scipy
import lal
from .utils import A_ASD_solutions, TFDG, EHG, xi_0_upper_bound, delta_xi_tilde_lower_bound
from bilby.core import utils
from bilby.core.series import CoupledTimeAndFrequencySeries
from bilby.core.utils import PropertyAccessor
from bilby.gw.conversion import convert_to_lal_binary_neutron_star_parameters



def dphi_prior(phase_uncertainty,k, **kwargs):
    '''
    Generates a Gaussian prior for the phase correction parameters for the BasicCorrectionModel
    
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
    Generates a Gaussian prior for the amplitude correction parameters for the BasicCorrectionModel
     
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



def xi_priors(waveform_generator,prior,psd_data,n,f_low,**kwargs):
    '''
    Generates xi_0 and delta_xi_tilde priors from a BBH/BNS/NSBH prior and adds them to the original prior.

    Parameters
    ==================
    waveform_generator: bilby.gw.WaveformGenerator
        bilby waveform generator object
    prior: bilby.core.prior.PriorDict
        bilby prior dictionary
    psd_data: numpy.ndarray
        array of power spectral density data; first column needs to be the frequency points and the second column needs to be the data
    n: int
        number of frequency nodes
    f_low: float
        lower bound on the frequency band (Hz)
    xi_0_latex_label: string, optional
        latex label for xi_0
        default: r'$xi_0$'
    delta_xi_tilde_latex_label: string, optional
        latex_label for delta_xi_tilde
        default: r'$delta tilde xi$'
    xi_low: float, optional
        lower bound on the dimensionless frequency band
        default: 0.018
    xi_high: float, optional
        upper bound on the dimensionless frequency band
        default: 1/pi
    samples: int, optional
        number of draws of amplitude to take to generate the priors
        default: 1000
    '''
    xi_0_latex_label = kwargs.get('xi_0_latex_label',r'$\xi_0$')
    delta_xi_tilde_latex_label = kwargs.get('delta_xi_tilde_latex_label',r'$\delta\tilde\xi$')
    xi_low = kwargs.get('xi_low',0.018)
    xi_high = kwargs.get('xi_high',1/np.pi)
    samples = kwargs.get('samples',1000)
    
    lower_xis, upper_xis = A_ASD_solutions(waveform_generator,psd_data,prior,samples,xi_low,xi_high)
    
    mu_1,sigma_1 = scipy.stats.norm.fit(lower_xis)
    mu_2,sigma_2 = scipy.stats.norm.fit(upper_xis)
    
    delta_xi_tildes = (np.array(upper_xis)-np.array(lower_xis))/(xi_high-np.array(lower_xis))
    
    mu_3,sigma_3 = scipy.stats.norm.fit(delta_xi_tildes)
    
    prior['xi_0'] = TFDG(name='xi_0',latex_label=xi_0_latex_label,
                        mu_1=mu_1,mu_2=mu_2,sigma_1=sigma_1,sigma_2=sigma_2,
                        minimum=xi_low,maximum=xi_0_upper_bound(n,xi_low=xi_low,xi_high=xi_high))
    
    prior['delta_xi_tilde'] = EHG(name='delta_xi_tilde',latex_label=delta_xi_tilde_latex_label,
                                  mu=mu_3,sigma=sigma_3,maximum=1,
                                  minimum=delta_xi_tilde_lower_bound(n,f_low,waveform_generator.duration,xi_low=xi_low,xi_high=xi_high))
    
    return prior



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
        default: r'$mathrm{M}_odot$' (solar mass)
    latex_label: string, optional
        label for the parameter in LaTeX
        default: r'$M$'
    boundary: string, optional
        boundary condition type for the prior
        default: None
    xi_low: float, optional
        lower bound on the waveform uncertainty correction in dimensionless frequency
        default: 0.018
    xi_high: float, optional
        upper bound on the waveform uncertainty correction in dimensionless frequency
        default: 1/pi, 0.318...

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
