recovery_from_parameterization
==============================

.. code-block:: python

   GWCorrect.parameterization.recovery_from_parameterization(parameterization_draw,dimensionless=False,
                                                                       xi_low=0.01,xi_high=1,resolution=1000)

Converts a draw of parameterized waveform differences back into waveform difference arrays.

Parameters:
-----------
parameterization_draw: numpy.ndarray
   one row of a parameterization matrix
dimensionless: bool, optional, (False)
   whether or not the output is returned in dimensionless frequency units
xi_low: float, optional, (0.001)
   if dimensionless is True, this is the lower bound on the dimensionless frequency grid
xi_high: float, optional, (1)
   if dimensionless is True, this is the upper bound on the dimensionless frequency grid
resolution: int, optional, (1000)
   if dimensionless is True, this is the number of points in the dimensionless frequency grid

Returns:
--------
frequency_grid: numpy.ndarray
   array of frequency points
amplitude_difference: numpy.ndarray
   array of amplitude differences
phase_difference: numpy.ndarray
   array of phase differences
