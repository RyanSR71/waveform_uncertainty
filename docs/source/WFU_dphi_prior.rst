WaveformUncertainty.WFU_dphi_prior
=============================

.. code-block:: python

   WaveformUncertainty.WFU_dphi_prior(phase_uncertainty,frequency_grid,injection,hf,PSDs,match_boundary,
                                      duration,nnodes,prior=None,polarization='plus',match_resolution=100)

Automatically generates a bilby prior object containing truncated Gaussian priors for each dphi parameter.

.. math::

   \mathcal{TN}(0,\delta\phi_\mu(f),\mp\varphi^{\gamma\%})=\frac{2(2\pi)^{-\frac{1}{2}}}{\delta\phi_{\mu}(f_k)}\left[\frac{\mathrm{exp}\left[-\frac{1}{2}\left(\frac{\Delta\phi}{\delta\phi_{\mu}(f_k)}\right)^2\right]}{\mathrm{erf}\left(\frac{\varphi_k^{\gamma\%}}{\sqrt{2}\delta\phi_{\mu}(f_k)}\right)+\mathrm{erf}\left(\frac{\varphi_k^{\gamma\%}}{\sqrt{2}\delta\phi_{\mu}(f_k)}\right)}\right]

.. math::

   f_\mathrm{IM}=\frac{0.018c^3}{GM}\quad\quad f_\mathrm{light}=\frac{c^3}{\pi GM}

.. math::

   \pi(\varphi_k)=\begin{cases}
        \delta(\Delta\phi) & f_k\leq f_\mathrm{IM} \\
        \mathcal{TN}(0,\delta\phi_{\mu}(f_k),\mp\varphi_k^{\gamma\%}) & f_\mathrm{IM}<f_k<f_\mathrm{light} \\
        \delta(\Delta\phi) & f_k\geq f_\mathrm{light}
    \end{cases}

Parameters:
-----------
phase_uncertainty: numpy.ndarray
   array of the phase uncertainty; defined as the standard deviation of phase differences
frequency_grid: numpy.ndarray
   frequency array corresponding to the waveform uncertainties
injection: dictionary
   dictionary of injection parameters
hf: bilby.gw.WaveformGenerator
   bilby waveform generator object
PSDs: numpy.ndarray
   power spectral density data
match_boundary: float
   required faithfullness for each dphi parameter
duration: int [s]
   duration of the signal
nnodes: int
   number of frequency nodes desired
prior: bilby.core.prior.PriorDict, optional, (None)
   if given, the output prior will simply be added to this
polarization: string, optional, ('plus')
   polarization of the signal
match_resolution: int, optional, (100)
   number of draws of the match function taken to find the faithfullness boundaries
      
Returns:
--------
prior: bilby.core.prior.PriorDict
   prior containing the dphi parameters
frequency_nodes: numpy.ndarray
   frequency nodes used by __WaveformGeneratorWFU() to generate waveform difference splines
indexes: numpy.ndarray
   list of dphi parameters present in the prior; i.e. dphi_1-dphi_5 gives indexes=[1,2,3,4,5]
