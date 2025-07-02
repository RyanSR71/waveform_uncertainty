uncertainties_from_parameterization
===================================

.. code-block:: python

   WaveformUncertainty.parameterization.uncertainties_from_parameterization(parameterization,dimensionless=False,
                                                                            xi_low=0.001,xi_high=1,resolution=1000)

Takes all of the sets in a parameterized waveform difference matrix and takes the mean and standard deviation of amplitude and phase difference.

Parameters:
-----------
parameterization: numpy.ndarry
   parameterization matrix
dimensionless: bool, optional
   whether or not the output is returned in dimensionless frequency units
   Default: False
xi_low: float, optional
   if dimensionless is True, this is the lower bound on the dimensionless frequency grid
   default: 0.001
xi_high: float, optional
   if dimensionless is True, this is the upper bound on the dimensionless frequency grid
   default: 1
resolution: int, optional
   if dimensionless is True, this is the number of points in the dimensionless frequency grid
   default: 1000
      
Returns:
--------
mean_amplitude_difference: numpy.ndarray
   array of the mean amplitude differences as a function of frequency
amplitude_uncertainty: numpy.ndarray
   array of the standard deviation of the amplitude differences across frequency
mean_phase_difference: numpy.ndarry
   array of the mean phase differences as a function of frequency
phase_uncertainty: numpy.ndarray
   array of the stand deviation of the phase differences across frequency
