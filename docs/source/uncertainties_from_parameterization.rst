WaveformUncertainty.uncertainties_from_parameterization
=======================================================

.. code-block:: python

   WaveformUncertainty.uncertainties_from_parameterization(data,linear=False,resolution=None)

Takes all of the sets in a parameterized waveform difference matrix and takes the mean and standard deviation of amplitude and phase difference

Parameters
==================
data: numpy.ndarray
    WaveformUncertainty.parameterization() output table
linear: bool, optional, (False)
    if True, default geometric frequency grid will be replaced with a linear one; useful for waveform uncertainty sampling
resolution: float, (None)
    distance between points in the linear frequency grid
      
Returns
==================
mean_amplitude_difference: numpy.ndarray
    array of the mean value of the amplitude difference corresponding to the frequency grid
amplitude_uncertainty: numpy.ndarray
    standard deviations of the amplitude difference
mean_phase_difference: numpy.ndarray
    array of the mean value of the phase difference corresponding to the frequency grid
phase_uncertainty: numpy.ndarray
    standard deviations of the phase difference
linear_frequency_grid: numpy.ndarray
    only if linear=True, new linear frequency grid
