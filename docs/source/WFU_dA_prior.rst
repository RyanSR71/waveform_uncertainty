WaveformUncertainty.WFU_dA_prior
=============================

.. code-block:: python

   WaveformUncertainty.WFU_dA_prior(amplitude_uncertainty,frequency_grid,injection,hf,PSDs,match_boundary,
                                      duration,nnodes,prior=None,polarization='plus',match_resolution=100)

Automatically generates a bilby prior object containing truncated Gaussian priors for each dA parameter.

.. math::

   \mathfrak{M}(h_1,h_2)=\underset{t_c,\phi_c}{\max}\hspace{0.1cm}\left|\frac{\int df\hat{h}_1^*(f)\hat{h}_2(f)\mathrm{e}^{-2\pi ift_c+i\phi_c}}{\sqrt{\int df\hat{h}_1^*(f)\hat{h}_1(f)}\sqrt{\int df\hat{h}_2^*(f)\hat{h}_2(f)}}\right|\times 100\%

.. math::

   \mathfrak{M}(\mu(f;\theta),\mu_{\mathcal{C}}(f;\theta,\alpha_k^{\gamma\%}))=\gamma\%

.. math::

   \mathcal{TN}(0,\delta\mathcal{A}_\mu(f),\mp\alpha^{\gamma\%})=\frac{(2\pi)^{-\frac{1}{2}}}{\delta\mathcal{A}_{\mu}(f_k)}\left[\frac{\mathrm{exp}\left[-\frac{1}{2}\left(\frac{\Delta\mathcal{A}}{\delta\mathcal{A}_{\mu}(f_k)}\right)^2\right]}{\mathrm{erf}\left(\frac{\alpha_k^{\gamma\%}}{\sqrt{2}\delta\mathcal{A}_{\mu}(f_k)}\right)}\right]

.. math::

   f_\mathrm{IM}=\frac{0.018c^3}{GM}\quad\quad f_\mathrm{light}=\frac{c^3}{\pi GM}

.. math::

   \pi(\alpha_k)=\begin{cases}
        \delta(\Delta\mathcal{A}) & f_k\leq f_\mathrm{IM} \\
        \mathcal{TN}(0,\delta\mathcal{A}_{\mu}(f_k),\mp\alpha_k^{\gamma\%}) & f_\mathrm{IM}<f_k<f_\mathrm{light} \\
        \delta(\Delta\mathcal{A}) & f_k\geq f_\mathrm{light}
    \end{cases}

Parameters:
-----------
amplitude_uncertainty: numpy.ndarray
   array of the amplitude uncertainty; defined as the standard deviation of amplitude differences
frequency_grid: numpy.ndarray
   frequency array corresponding to the waveform uncertainties
injection: dictionary
   dictionary of injection parameters
hf: bilby.gw.WaveformGenerator
   bilby waveform generator object
PSDs: numpy.ndarray
   power spectral density data
match_boundary: float
   required faithfullness for each dA parameter
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
   prior containing the dA parameters
frequency_nodes: numpy.ndarray
   frequency nodes used by __WaveformGeneratorWFU() to generate waveform difference splines
indexes: numpy.ndarray
   list of dA parameters present in the prior; i.e. dA_1-dA_5 gives indexes=[1,2,3,4,5]
