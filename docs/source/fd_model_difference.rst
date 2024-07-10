WaveformUncertainty.fd_model_difference
=======================================

.. code-block:: python

   WaveformUncertainty.fd_model_difference(hf1,hf2,f_low=20.0,f_high=2048.0,f_ref=50.0,
                                             npoints=1000,polarization='plus',psd_data=None,
                                             correction_parameter=-10e-6)

Generates frequency domain waveform differences between two models hf1 and hf2

.. _parameters:

Parameters:
-----------
hf1: bilby.gw.waveform_generator.WaveformGenerator
   frequency domain waveform generator object WITH injected parameters (for strain calculation)
hf2: bilby.gw.waveform_generator.WaveformGenerator
   frequency domain waveform generator object WITH injected parameters (for strain calculation)
f_low: float, optional, (20.0)
   minimum frequency
f_high: float, optional, (2048.0)
   maximum frequency
f_ref: float, optional, (50.0)
   reference frequency
npoints: int, optional, (1000)
   length of the desired frequency grid
polarization: string, optional, ('plus')
   polarization of the strain data (plus or cross)
psd_data: numpy.ndarray, optional, (None)
   array containing the psd data and their corresponding frequencies
correction_parameter: float, optional, (-10e-6)
   value at which to cut the second derivative of amplitude difference (see WFU_equations.pdf #1)

.. returns:

Returns:
--------
frequency_grid: numpy.ndarray
   frequency grid that corresponds to the uncertainty arrays
amplitude_difference: numpy.ndarray
   array of amplitude difference values
phase_difference: numpy.ndarray
   array of phase difference values; if psd data is None, raw_phase_difference will be returned, residual_phase_difference otherwise
amplitude_difference_final_point: float
   value of the amplitude difference where the discontinuity correction starts
phase_difference_final_point: float
   value of the phase difference where the discontinuity correction starts
final_index: int
   position within the frequency grid where the discontinuity correction starts
