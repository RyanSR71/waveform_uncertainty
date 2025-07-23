match
=====

.. code-block::

  GWCorrect.ppE.match(signal,data,PSDs,duration)

Computes the normalized match between two waveforms.

.. math::

  \mathrm{match}(\tilde{h}_1,\tilde{h}_2)=\underset{\Delta t_c,\Delta\phi_c}{\mathrm{max}}\left|\frac{\int\hat{h}_1^*(f)\hat{h}_2(f)\exp(-2\pi if\Delta t_c+i\Delta\phi_c)df}{\sqrt{\int \hat{h}_1^*(f)\hat{h}_1(f)df}\sqrt{\int \hat{h}_2^*(f)\hat{h}_2(f)df}}\right|

Parameters:
-----------
signal: numpy.ndarray
    strain array of the first gravitational wave
data: numpy.ndarray
    strain array of the second gravitational wave
PSDs: numpy.ndarray
    array of power spectral densities to weight the waveforms
duration: float
    time duration of the gravitational waves                                                                                                                                         
