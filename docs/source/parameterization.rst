WaveformUncertainty.parameterization
====================================

.. code-block:: python

   WaveformUncertainty.parameterization(hf1,hf2,parameter_data,nsamples,
                                        precession=False,tides=True,fit_parameters=15,
                                        npoints=1000,max_amplitude_error=2,max_phase_error=2,
                                        psd_data=None,correction_parameter_A=5e-5,
                                        correction_parameter_B=0,correction_parameter_C=2,
                                        ref_amplitude=None,polarization='plus',
                                        fit_threshold=75)

Generates samples of waveform differences between two approximants and parameterizes the data (See `Equations and Notation <https://waveformuncertainty.readthedocs.io/en/latest/WFU_Equations.html#parameterization>`_)

Parameters:
-----------
hf1: bilby.gw.waveform_generator.WaveformGenerator
    frequency domain waveform generator object
hf2: bilby.gw.waveform_generator.WaveformGenerator
    frequency domain waveform generator object
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
max_dA_error: float [%], optional, (2)
    maximum allowed error between the amplitude uncertainty and its parameterization
max_dphi_error: float [degrees], optional, (2)
    maximum allowed error between the phase uncertainty and its parameterization
psd_data: numpy.ndarray, optional, (None)
    array containing the psd data and their corresponding frequencies
correction_parameter_A: float, optional, (5e-5)
    value at which to cut the second derivative of amplitude difference; if None, correction will not occur
correction_parameter_B: int, optional, (0)
    index at which to start the search for any discontinuity
correction_parameter_C: int, optional, (2)
    number of amplitude difference derivatives for the discontinuity correction
ref_amplitude: numpy.ndarray, optional, (None)
   reference amplitude for residual phase calculation; will be generated automatically if not given
polarization: string, optional, ('plus')
    polarization of the strain data {'plus','cross'}
fit_threshold: float [%], optional, (75)
    minimum parameterization success rate
  
Returns:
--------
parameterized_data: numpy.ndarray
    table containing the index, frequency_grid, dA_fit_parameters, dphi_fit_parameters, final_index, dA_final_point, dphi_final_point,
    and injection_parameters for each draw of waveform difference
      
    index: int
        position within the list of indexes the waveform difference draws were drawn from (only for debugging purposes)
    frequency_grid: numpy.ndarray
        frequencies corresponding to the frequency parameters specified
    dA_fit_parameters: numpy.ndarray
        parameters corresponding to the parameterization of the amplitude difference
    dphi_fit_parameters: numpy.ndarray
        parameters corresponding to the parameterization of the phase difference
    final_index: int
        position within the frequency grid where the discontinuity correction occurs
    dA_final_point: float
        value of the amplitude difference at the discontinuity correction
    dphi_final_point: float
        value of the phase difference at the discontinuity correction
    injection parameters: dictionary
        neutron star parameters injected into the waveform generators
