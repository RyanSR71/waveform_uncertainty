WaveformUncertainty.parameterization
====================================

.. code-block:: python

   WaveformUncertainty.parameterization(hf1,hf2,prior,nsamples,fit_parameters=15,
                                        npoints=1000,max_amplitude_error=1,max_phase_error=5,
                                        psd_data=None,correction_parameter_A=0.01,
                                        correction_parameter_B=256,correction_parameter_C=2,
                                        ref_amplitude=None,polarization='plus',
                                        fit_threshold=75)

Generates samples of waveform differences between two approximants and parameterizes the data (See `Equations and Notation <https://waveformuncertainty.readthedocs.io/en/latest/WFU_Equations.html#parameterization>`_)

.. math::

   \Delta\mathcal{A}_{\mu}(f;\theta)\approx\Delta\mathcal{A}_{T}(f;a,f_{\mathrm{disc}},\Delta\mathcal{A}_{\mu}(f_{\mathrm{disc}};\theta))= \begin{cases} 
      \sum_{i=0}^{N-1}a_{i}T_{i}(f) & f \leq f_{\mathrm{disc}} \\
      \Delta\mathcal{A}_{\mu}(f_{\mathrm{disc}};\theta) & f > f_{\mathrm{disc}} 
   \end{cases}

.. math::

   \Delta\phi_{\mu}(f;\theta)\approx\Delta\phi_{T}(f;p,f_{\mathrm{disc}},\Delta\phi_{\mu}(f_{\mathrm{disc}};\theta))= \begin{cases} 
      \sum_{i=0}^{N-1}p_{i}T_{i}(f) & f \leq f_{\mathrm{disc}} \\
      \Delta\phi_{\mu}(f_{\mathrm{disc}};\theta) & f > f_{\mathrm{disc}} 
   \end{cases}

Parameters:
-----------
hf1: bilby.gw.waveform_generator.WaveformGenerator
    frequency domain waveform generator object
hf2: bilby.gw.waveform_generator.WaveformGenerator
    frequency domain waveform generator object
prior: bilby.core.prior.dict.PriorDict
    bilby prior object
nsamples: int
    number of draws of waveform uncertainty desired
fit_parameters: int, optional, (15)
    number of terms to use in the parameterization
npoints: int, optional, (1000)
    length of the desired frequency grid
max_amplitude_error: float [%], optional, (1)
    maximum allowed error between the amplitude uncertainty and its parameterization
max_phase_error: float [degrees], optional, (5)
    maximum allowed error between the phase uncertainty and its parameterization
psd_data: numpy.ndarray, optional, (None)
    array containing the psd data and their corresponding frequencies
correction_parameter_A: float, optional, (0.01)
    value at which to cut the second derivative of amplitude difference; if None, correction will not occur
correction_parameter_B: int, optional, (256)
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
