match
=====

.. code-block:: python

  GWCorrect.ppE.prior.match(signal,data,duration,PSDs=None)

Computes the normalized match between two waveforms.

.. math::

  \mathrm{match}(\tilde{h}_1,\tilde{h}_2)=\underset{\Delta t_c,\Delta\phi_c}{\mathrm{max}}\left|\frac{\int\hat{h}_1^*(f)\hat{h}_2(f)\exp(-2\pi if\Delta t_c+i\Delta\phi_c)df}{\sqrt{\int \hat{h}_1^*(f)\hat{h}_1(f)df}\sqrt{\int \hat{h}_2^*(f)\hat{h}_2(f)df}}\right|

Parameters:
-----------
signal: numpy.ndarray
    strain array of the first gravitational wave
data: numpy.ndarray
    strain array of the second gravitational wave
duration: float
    time duration of the gravitational waves 
PSDs: numpy.ndarray, optional, (None)
    array of power spectral densities to weight the waveforms

Returns:
--------
match: float
    the match between the two waveforms; ranges from 0 to 1
