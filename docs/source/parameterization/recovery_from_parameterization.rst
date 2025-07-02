recovery_from_parameterization
==============================

.. code-block:: python

   WaveformUncertainty.parameterization.recovery_from_parameterization(parameterization_draw,dimensionless=False,
                                                                       xi_low=0.01,xi_high=1,resolution=1000)

Converts a draw of parameterized waveform differences back into waveform difference arrays.

Parameters:
-----------
parameterization_draw: numpy.ndarray
   one row of a parameterization matrix
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
frequency_grid: numpy.ndarray
   array of frequency points
amplitude_difference: numpy.ndarray
   array of amplitude differences
phase_difference: numpy.ndarray
   array of phase differences
