WaveformUncertainty.fd_model_difference
=======================================

.. code-block:: python

   WaveformUncertainty.fd_model_difference(hf1,hf2,injection=None,npoints=1000,
                                           polarization='plus',psd_data=None,
                                           correction_parameter=0.0001,ref_amplitude=None)

Generates frequency domain waveform differences between two models hf1 and hf2.

.. math::

   \Delta\mathcal{A}_{\mu}(f;\vartheta)= \begin{cases} 
      \frac{\mathcal{A}_2(f;\vartheta)}{\mathcal{A}_1(f;\vartheta)} & f \leq f_{\mathrm{disc}} \\
      \Delta\mathcal{A}_\mu(f_{\mathrm{disc}};\vartheta) & f > f_{\mathrm{disc}} 
   \end{cases}

.. math::

.. math::

   \Delta\phi_{\mu}(f;\vartheta)= \begin{cases} 
      \phi_2(f;\vartheta)-\phi_1(f;\vartheta)-2\pi f\Delta t_c-\Delta\phi_c & f \leq f_{\mathrm{disc}} \\
      \Delta\phi_\mu(f_{\mathrm{disc}};\vartheta) & f > f_{\mathrm{disc}} 
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
correction_parameter: float, optional, (0.0001)
   fraction of the peak frequency domain amplitudes at which to cut off the amplitudes; part of the f_disc calculation
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
