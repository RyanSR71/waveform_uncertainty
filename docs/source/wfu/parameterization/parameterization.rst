parameterization
================

.. code-block:: python

   GWCorrect.parameterization.parameterization(hf1,hf2,prior,nsamples,fit_parameters=15,
                                        npoints=1000,max_amplitude_error=1,max_phase_error=5,
                                        psd_data=None,correction_parameter=0.0001,
                                        ref_amplitude=None,polarization='plus',
                                        fit_threshold=75)

Generates samples of waveform differences between two approximants and parameterizes the data.

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
spline_resolution: int, 500
   number of spline nodes desired
   default: 500
npoints: int, optional
   length of the desired frequency grid
   default: 1000
psd_data: numpy.ndarray, optional
   array containing the psd data and their corresponding frequencies
   default: None
correction_parameter: float, optional
   fraction of maximum amplitude to cut off the amplitude at
   default: 0.0001
ref_amplitude: numpy.ndarray, optional
   reference amplitude for residual phase calculation
   default: None
polarization: string, optional
   polarization of the strain data (plus or cross)
   default: 'plus'
  
Returns:
--------
parameterized_data: numpy.ndarray
   table containing the index, frequency_grid, dA_fit_parameters, dphi_fit_parameters, 
   final_index, dA_final_point, dphi_final_point, and injection_parameters for each draw
   
   frequency_grid: numpy.ndarray
      frequencies corresponding to the frequency parameters specified
   frequency_nodes: numpy.ndarray
      frequency nodes for the splines
   dA_parameters: numpy.ndarray
      amplitude difference spline parameters
   dphi_parameters: numpy.ndarray
      phase difference spline parameters
   injection: dictionary
      source parameters injected into the waveform generators
