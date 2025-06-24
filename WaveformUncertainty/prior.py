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



def xi_0_upper_bound(n, **kwargs):
    xi_low = kwargs.get('xi_low',0.018)
    xi_high = kwargs.get('xi_high',1/np.pi)
    def f(x):
        return x**(1-n) + (xi_low**(1-n)/(xi_high-xi_low))*x - ((xi_high*xi_low**(1-n))/(xi_high-xi_low))
    root = scipy.optimize.brentq(f, xi_low, xi_high)
    return root



class TFDG(bilby.core.prior.Prior):
    def __init__(self,mu_1,mu_2,sigma_1,sigma_2,minimum,maximum,name=None, latex_label=None):
        super(TFDG, self).__init__(
            name=name,latex_label=latex_label,minimum=minimum,maximum=maximum
        )
        self.mu_1 = float(mu_1)
        self.mu_2 = float(mu_2)
        self.sigma_1 = float(sigma_1)
        self.sigma_2 = float(sigma_2)
        
        
    def prob(self, val):
        in_region_1 = (val >= self.minimum) & (val <= self.mu_1)
        in_region_2 = (val > self.mu_1) & (val < self.mu_2)
        in_region_3 = (val >= self.mu_2) & (val <= self.maximum)
        N = (-np.sqrt(np.pi/2)*self.sigma_1*scipy.special.erf((self.minimum-self.mu_1)/(np.sqrt(2)*self.sigma_1))+np.sqrt(np.pi/2)*self.sigma_2*scipy.special.erf((self.maximum-self.mu_2)/(np.sqrt(2)*self.sigma_2))+self.mu_2-self.mu_1)**-1
        draw = N*np.exp(-0.5*((val-self.mu_1)/self.sigma_1)**2)*in_region_1+N*in_region_2+N*np.exp(-0.5*((val-self.mu_2)/self.sigma_2)**2)*in_region_3
        return draw
    
    
    def rescale(self, val):
        N = (-np.sqrt(np.pi/2)*self.sigma_1*scipy.special.erf((self.minimum-self.mu_1)/(np.sqrt(2)*self.sigma_1))+np.sqrt(np.pi/2)*self.sigma_2*scipy.special.erf((self.maximum-self.mu_2)/(np.sqrt(2)*self.sigma_2))+self.mu_2-self.mu_1)**-1
        A_1 = np.sqrt(np.pi/2)*N*self.sigma_1*scipy.special.erf((self.mu_1-self.minimum)/(np.sqrt(2)*self.sigma_1))
        A_2 = N*(self.mu_2-self.mu_1)
        
        if hasattr(val, "__len__"):
            draw = []
            for v in val:
        
                in_region_1 = (v >= 0) & (v <= A_1)
                in_region_2 = (v > A_1) & (v <= A_1+A_2)
                in_region_3 = (v >= A_1+A_2) & (v <= 1)

                if in_region_1:
                    draw.append(np.sqrt(2)*self.sigma_1*scipy.special.erfinv((np.sqrt(2*np.pi)*v-np.pi*np.sqrt(2/np.pi)*A_1)/(np.pi*N*self.sigma_1))+self.mu_1)
                elif in_region_2:
                    draw.append((N*self.mu_1+v-A_1)/N)
                elif in_region_3:
                    draw.append(self.mu_2-np.sqrt(2)*self.sigma_2*scipy.special.erfinv((np.sqrt(2*np.pi)*A_1+np.sqrt(2*np.pi)*A_2-np.sqrt(2*np.pi)*v)/(np.pi*N*self.sigma_2)))
                else:
                    raise Exception('Draw Failed!')
            return np.array(draw)
        
        else:
            in_region_1 = (val >= 0) & (val <= A_1)
            in_region_2 = (val > A_1) & (val <= A_1+A_2)
            in_region_3 = (val >= A_1+A_2) & (val <= 1)

            if in_region_1:
                draw = np.sqrt(2)*self.sigma_1*scipy.special.erfinv((np.sqrt(2*np.pi)*val-np.pi*np.sqrt(2/np.pi)*A_1)/(np.pi*N*self.sigma_1))+self.mu_1
            elif in_region_2:
                draw = (N*self.mu_1+val-A_1)/N
            elif in_region_3:
                draw = self.mu_2-np.sqrt(2)*self.sigma_2*scipy.special.erfinv((np.sqrt(2*np.pi)*A_1+np.sqrt(2*np.pi)*A_2-np.sqrt(2*np.pi)*val)/(np.pi*N*self.sigma_2))
            else:
                raise Exception('Draw Failed!')
            return draw



class EHG(bilby.core.prior.Prior):
    def __init__(self,mu_1,mu_2,sigma_1,sigma_2,minimum,maximum,name=None, latex_label=None):
        super(TFDG, self).__init__(
            name=name,latex_label=latex_label,minimum=minimum,maximum=maximum
        )
        self.mu = float(mu)
        self.sigma = float(sigma)        
        
    def prob(self, val):
        in_region_1 = (val >= self.minimum) & (val <= self.mu)
        in_region_2 = (val >= self.mu) & (val <= self.maximum)
        N = (self.mu-self.minimum-np.sqrt(np.pi/2)*scipy.special.erf((self.mu-self.maximum)/(np.sqrt(2)*self.sigma)))**(-1)
        draw = N*in_region_1+N*np.exp(-0.5*((val-self.mu)/self.sigma)**2)*in_region_2
        return draw
    
    
    def rescale(self, val):
        N = (self.mu-self.minimum-np.sqrt(np.pi/2)*scipy.special.erf((self.mu-self.maximum)/(np.sqrt(2)*self.sigma)))**(-1)
        A = N*(self.mu-self.minimum)
        
        if hasattr(val, "__len__"):
            draw = []
            for v in val:
        
                in_region_1 = (v >= 0) & (v <= A)
                in_region_2 = (v > A) & (v <= 1)

                if in_region_1:
                    draw.append((val+N*self.minimum)/N)
                elif in_region_2:
                    draw.append(self.mu+np.sqrt(2)*self.sigma*scipy.special.erfinv(np.sqrt(2/np.pi)*(val-A)/(N*self.sigma)))
                else:
                    raise Exception('Draw Failed!')
            return np.array(draw)
        
        else:
            in_region_1 = (val >= 0) & (val <= A)
            in_region_2 = (val > A) & (val <= 1)

            if in_region_1:
                draw = (val+N*self.minimum)/N
            elif in_region_2:
                draw = self.mu+np.sqrt(2)*self.sigma*scipy.special.erfinv(np.sqrt(2/np.pi)*(val-A)/(N*self.sigma))
            else:
                raise Exception('Draw Failed!')
            return draw
