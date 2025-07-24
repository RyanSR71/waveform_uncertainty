import numpy as np
from numpy.linalg import inv
import bilby
import random
import time as tm
import sys
import scipy
import lal
import math
from .utils import inversion_function, beta_from_beta_tilde_wrapped, apply_ppe_correction

def ppECorrectionModel(
        frequency_array, mass_1, mass_2, luminosity_distance, a_1, tilt_1,
        phi_12, a_2, tilt_2, phi_jl, theta_jn, phase,
        beta_tilde, delta_epsilon_tilde, b, **kwargs):
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
    Mf_IM = kwargs.get('Mf_IM',0.018)
    ratio_f_MR_to_f_RD = kwargs.get('ratio_f_MR_to_f_RD',0.75)
    aligned = kwargs.get('aligned',True)
    
    waveform_arguments = dict(
        waveform_approximant=waveform_approximant, reference_frequency=reference_frequency,
        minimum_frequency=minimum_frequency, maximum_frequency=maximum_frequency,
        catch_waveform_errors=catch_waveform_errors, pn_spin_order=pn_spin_order, pn_tidal_order=pn_tidal_order,
        pn_phase_order=pn_phase_order, pn_amplitude_order=pn_amplitude_order,
    )
    
    unmodified_strain = bilby.gw.source._base_lal_cbc_fd_waveform(
        frequency_array=frequency_array, mass_1=mass_1, mass_2=mass_2,
        luminosity_distance=luminosity_distance, theta_jn=theta_jn, phase=phase,
        a_1=a_1, a_2=a_2, tilt_1=tilt_1, tilt_2=tilt_2, phi_12=phi_12,
        phi_jl=phi_jl, **waveform_arguments)
    
    beta = inversion_function(0,3,b)*beta_from_beta_tilde_wrapped(beta_tilde,minimum_frequency,1/np.pi,b,0.018,mass_1+mass_2)
    delta_epsilon = inversion_function(0,3,b)*beta_from_beta_tilde_wrapped(delta_epsilon_tilde,minimum_frequency,1/np.pi,b,0.018,mass_1+mass_2)
    
    model_strain = dict()
    model_strain['plus'] = apply_ppe_correction(unmodified_strain['plus'],frequency_array,mass_1+mass_2,beta,b,delta_epsilon,
                                                Mfreq_IM=Mfreq_IM,ratio_f_MR_to_f_RD=ratio_f_MR_to_f_RD,aligned=aligned,window_size=41)
    model_strain['cross'] = apply_ppe_correction(unmodified_strain['cross'],frequency_array,mass_1+mass_2,beta,b,delta_epsilon,
                                                 Mfreq_IM=Mfreq_IM,ratio_f_MR_to_f_RD=ratio_f_MR_to_f_RD,aligned=aligned,window_size=41)
    
    return model_strain
