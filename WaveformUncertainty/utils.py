import numpy as np
import bilby
import random
import time
import sys
import scipy
import lal
import tqdm
import logging
from bilby.core import utils
from bilby.core.series import CoupledTimeAndFrequencySeries
from bilby.core.utils import PropertyAccessor
from bilby.gw.conversion import convert_to_lal_binary_neutron_star_parameters



def progressBar(count_value, total, suffix=''):
        bar_length = 100
        filled_up_Length = int(round(bar_length* count_value / float(total)))
        percentage = round(100.00 * count_value/float(total),1)
        bar = '=' * filled_up_Length + '-' * (bar_length - filled_up_Length)
        sys.stdout.write('[%s] %s%s %s\r' %(bar, percentage, '%', suffix))
        sys.stdout.flush()



class ProgressBar(logging.Handler):
    def __init__(self, level=logging.NOTSET):
        super().__init__(level)

    def emit(self, record):
        try:
            msg = self.format(record)
            tqdm.tqdm.write(msg)
            self.flush()
        except Exception:
            self.handleError(record)



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



def smooth_interpolation(full_grid,nodes,parameters,gamma):
    spline = scipy.interpolate.interp1d(nodes,parameters)(nodes)
    temp_grid = np.geomspace(full_grid[1],full_grid[-1],200)
    data = np.interp(temp_grid,nodes,spline)
    new_data = data.copy()
    
    r = int(gamma*len(data))
    if r != 0:
        lower_index = int(gamma*len(data))
        upper_index = int((1-gamma)*len(data))
        for i in range(lower_index,upper_index):
            new_data[i] = (1/(2*r))*np.sum(data[i-r:i+r])
    
    output = np.interp(full_grid,temp_grid,new_data)
    if np.abs(output[0]) > 0:
        output -= output[0]

    return output



def variable_prior(uncertainty,k,xi_low,xi_high):
    frequency_grid = np.linspace(0.001,1,len(uncertainty))
    desired_frequency_nodes = np.geomspace(xi_low,xi_high,k+1)
    
    indexes = [list(frequency_grid).index(min(frequency_grid, key=lambda x:np.abs(x-node))) for node in desired_frequency_nodes]
    frequency_nodes = np.array(frequency_grid[indexes])
    
    coef = np.array([uncertainty[indexes[i]] for i in range(k+1)])

    return frequency_nodes, coef



def GC_waveform_correction(frequency_array,xi_0,delta_xi_tilde,dAs,dphis,sigma_dA_spline,sigma_dphi_spline,mass_1,mass_2,xi_high,gamma):
    n = len(dAs)-1
    total_mass = mass_1+mass_2
    dimensionless_frequency_nodes = np.array([xi_0*(1+((xi_high-xi_0)/(xi_0))*delta_xi_tilde)**(k/n) for k in range(n+1)])
    frequency_nodes = dimensionless_frequency_nodes*203025.4467280836/total_mass
    if sigma_dA_spline is not None:
        sigma_dA = sigma_dA_spline(dimensionless_frequency_nodes)
    else: 
        sigma_dA = np.ones(n+1)
    if sigma_dphi_spline is not None:
        sigma_dphi = sigma_dphi_spline(dimensionless_frequency_nodes)
    else:
        sigma_dphi = np.ones(n+1)
    
    amplitude_correction = smooth_interpolation(frequency_array,frequency_nodes,dAs*sigma_dA,gamma)
    phase_correction = smooth_interpolation(frequency_array,frequency_nodes,dphis*sigma_dphi,gamma)
    
    waveform_correction = (1+amplitude_correction)*np.exp(phase_correction*1j)
    return waveform_correction



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



def A_ASD_solutions(waveform_generator,psd_data,prior,samples,xi_low,xi_high,desc):
    lower_xis = []
    upper_xis = []
    log = logging.getLogger(__name__)
    log.setLevel(logging.INFO)
    log.addHandler(ProgressBar())

    for trial in tqdm.tqdm(range(samples), desc = desc):
        roots = []
        injection = prior.sample()
        geometrized_frequency_grid = np.geomspace(xi_low,xi_high,1000)
        amplitude = np.abs(waveform_generator.frequency_domain_strain(parameters=injection)['plus'])
        M = bilby.gw.conversion.generate_mass_parameters(injection)['total_mass']
        freqs = waveform_generator.frequency_array/float(203025.4467280836/M)

        effective_amplitude = np.interp(geometrized_frequency_grid,freqs,2*amplitude*np.sqrt(freqs)*np.sqrt(float(203025.4467280836/M)))
        ASD = np.interp(geometrized_frequency_grid,psd_data[:,0]/float(203025.4467280836/M),np.sqrt(psd_data[:,1]))
        nodes = np.linspace(0,len(geometrized_frequency_grid)-1,100).astype(int)
        parameters = (effective_amplitude-ASD)[nodes]
        spline = scipy.interpolate.CubicSpline(geometrized_frequency_grid[nodes],parameters)
        roots = spline.roots()
        try:
            if roots[0] >= xi_low and roots[-1] <= xi_high:
                lower_xis.append(roots[0])
                upper_xis.append(roots[-1])
        except:
            pass
        progress += 1
    return lower_xis, upper_xis



def xi_0_upper_bound(n, **kwargs):
    xi_low = kwargs.get('xi_low',0.018)
    xi_high = kwargs.get('xi_high',1/np.pi)
    def f(x):
        return x**(1-n) + (xi_low**(1-n)/(xi_high-xi_low))*x - ((xi_high*xi_low**(1-n))/(xi_high-xi_low))
    # we need to have the lower bound slightly larger than xi_low, so we multiply by 1.000001 arbitrarily
    root = scipy.optimize.brentq(f, 1.000001*xi_low, xi_high)
    return root



def delta_xi_tilde_lower_bound(n,f_low,duration,**kwargs):
    xi_low = kwargs.get('xi_low',0.018)
    xi_high = kwargs.get('xi_high',1/np.pi)
    
    minimum = (xi_low/(xi_high-xi_low))*(4/(duration*f_low))**n
    
    return minimum



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
    def __init__(self,mu,sigma,minimum,maximum,name=None, latex_label=None):
        super(EHG, self).__init__(
            name=name,latex_label=latex_label,minimum=minimum,maximum=maximum
        )
        self.mu = float(mu)
        self.sigma = float(sigma)        
        
    def prob(self, val):
        in_region_1 = (val >= self.minimum) & (val <= self.mu)
        in_region_2 = (val > self.mu) & (val <= self.maximum)
        N = (np.sqrt(np.pi/2)*self.sigma*scipy.special.erf((self.maximum-self.mu)/(np.sqrt(2)*self.sigma))+self.mu-self.minimum)**-1
        draw = N*in_region_1+N*np.exp(-0.5*((val-self.mu)/self.sigma)**2)*in_region_2
        return draw
            
    
    def rescale(self, val):
        N = (np.sqrt(np.pi/2)*self.sigma*scipy.special.erf((self.maximum-self.mu)/(np.sqrt(2)*self.sigma))+self.mu-self.minimum)**-1
        A = N*(self.mu-self.minimum)
        
        if hasattr(val, "__len__"):
            draw = []
            for v in val:
                in_region_1 = (v >= 0) & (v < A)
                in_region_2 = (v >= A) & (v <= 1)

                if in_region_1:
                    draw.append((N*self.minimum+v)/N)
                elif in_region_2:
                    draw.append(self.mu-np.sqrt(2)*self.sigma*scipy.special.erfinv((np.sqrt(2*np.pi)*A-np.sqrt(2*np.pi)*v)/(np.pi*N*self.sigma)))
                else:
                    raise Exception('Draw Failed!')
            return np.array(draw)
        
        else:
            in_region_1 = (val >= 0) & (val < A)
            in_region_2 = (val >= A) & (val <= 1)

            if in_region_1:
                draw = ((N*self.minimum+val)/N)
            elif in_region_2:
                draw = (self.mu-np.sqrt(2)*self.sigma*scipy.special.erfinv((np.sqrt(2*np.pi)*A-np.sqrt(2*np.pi)*val)/(np.pi*N*self.sigma)))
            else:
                raise Exception('Draw Failed!')
            return draw
