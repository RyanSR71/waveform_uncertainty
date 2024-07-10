WaveformUncertainty.parameterization
====================================

.. code-block:: python

   WaveformUncertainty.parameterization(approximant1,approximant2,parameter_data,nsamples,
                                        precession=False,tides=True,fit_parameters=15,
                                        npoints=1000,f_low=20.0,f_high=2048.0,f_ref=50.0,
                                        sampling_frequency=4096,duration=256,
                                        max_amplitude_error=2,max_phase_error=2,psd_data=None,
                                        correction_parameter=-10e-6,polarization='plus',
                                        fit_threshold=75)

Generates samples of waveform differences between two approximants and parameterizes the data.

Parameters:
-----------
approximant1: string
    name of the first waveform approximant
approximant2: string
    name of the second waveform approximant
parameter_data: dictionary or bilby.core.prior.dict.PriorDict
    dictionary containing neutron star parameter samples or a bilby prior object that will be converted into a dictionary
nsamples: int
    number of draws of waveform uncertainty desired
precession: bool, optional, (False)
    True if both waveform approximants support precessing spins; 
    if False, precessing spin parameters will be removed/replaced with non-precessing parameters
tides: bool, optional, (True)
    only True if both waveform approximants support tidal deformabilities;
    if False, tidal parameters will be removed
fit_parameters: int, optional, (15)
    number of terms to use in the parameterization
npoints: int, optional, (1000)
    length of the desired frequency grid
f_low: float, optional, (20.0)
    lower bound on the frequency grid
f_high: float, optional, (2048.0)
    upper bound on the frequency grid
f_ref: float, optional, (50.0)
    reference frequency
sampling_frequency: float, optional, (4096)
    sampling frequency or sampling rate
duration: float, optional, (256)
    duration of the signal
max_dA_error: float [%], optional, (2)
    maximum allowed error between the amplitude uncertainty and its parameterization
max_dphi_error: float [degrees], optional, (2)
    maximum allowed error between the phase uncertainty and its parameterization
psd_data: numpy.ndarray, optional, (None)
    array containing the psd data and their corresponding frequencies
correction_parameter: float, optional, (-10e-6)
    value at which to cut the second derivative of amplitude difference (see WFU Equations #1)
polarization: string, optional, ('plus')
    polarization of the strain data {'plus','cross'}
fit_threshold: float [%], optional, (75)
    minimum parameterization success rate
  
Returns:
--------
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
