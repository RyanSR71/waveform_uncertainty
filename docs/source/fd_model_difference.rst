WaveformUncertainty.fd_model_difference
=======================================

.. code-block:: python

   WaveformUncertainty.fd_model_difference(hf1,hf2,injection=None,npoints=1000,
                                           polarization='plus',psd_data=None,
                                           correction_parameter_A=0.01,correction_parameter_B=256,
                                           correction_parameter_C=2,ref_amplitude=None)

Generates frequency domain waveform differences between two models hf1 and hf2.

.. math::

   \Delta\mathcal{A}_{\mu}(f;\theta)= \begin{cases} 
      \frac{|\mu_2(f;\theta)|-|\mu_1(f;\theta)|}{|\mu_1(f)|} & f \leq f_{\mathrm{disc}} \\
      \Delta\mathcal{A}_\mu(f_{\mathrm{disc}};\theta) & f > f_{\mathrm{disc}} 
   \end{cases}

.. math::

   \Delta\phi_{\mu}(f;\theta)= \begin{cases} 
      \arctan\left(\frac{\mathrm{Im}[\mu_2(f)]}{\mathrm{Re}[\mu_2(f)]}\right)-\arctan\left(\frac{\mathrm{Im}[\mu_1(f)]}{\mathrm{Re}[\mu_1(f)]}\right)-2\pi ft_c-\phi_c & f \leq f_{\mathrm{disc}} \\
      \Delta\phi_\mu(f_{\mathrm{disc}};\theta) & f > f_{\mathrm{disc}} 
   \end{cases}

Parameters:
-----------
hf1: bilby.gw.waveform_generator.WaveformGenerator
   frequency domain waveform generator object
hf2: bilby.gw.waveform_generator.WaveformGenerator
   frequency domain waveform generator object
injection: dictionary, optional, None
   dictionary of injection parameters if waveform generators do not have parameters; if they do not, this argument is not optional 
npoints: int, optional, (1000)
   length of the desired frequency grid
polarization: string, optional, ('plus')
   polarization of the strain data {'plus','cross'}
psd_data: numpy.ndarray, optional, (None)
   array containing the psd data and their corresponding frequencies
correction_parameter_A: float, optional, (0.01)
   value at which to cut the second derivative of amplitude difference; if None, correction will not occur
correction_parameter_B: int, optional, (256)
   index at which to start the search for any discontinuity
correction_parameter_C: int, optional, (2)
   number of amplitude difference derivatives to take for the discontinuity correction
ref_amplitude: numpy.ndarray, optional, (None)
   reference amplitude for residual phase calculation; will be generated automatically if not given

Returns:
--------
frequency_grid: numpy.ndarray
   array of frequencies that corresponds to the waveform difference arrays
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
