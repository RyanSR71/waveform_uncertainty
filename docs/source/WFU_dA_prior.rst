WaveformUncertainty.WFU_dA_prior
=============================

.. code-block:: python

   WaveformUncertainty.WFU_dA_prior(amplitude_uncertainty,k,mean_amplitude_difference=None,
                                    prior=None,geometrized=True,xi_low=0.018,xi_high=0.318,
                                    f_low=20.0,f_high=1024.0)

Automatically generates a bilby prior object containing truncated Gaussian priors for each dA parameter.

.. math::

   \Pi(\alpha_k)=\mathcal{N}(0,\delta\mathcal{A}_\mu(f_k))

Parameters:
-----------
amplitude_uncertainty: numpy.ndarray
   array of standard deviation of a set of amplitude differences; by default, this should be as a function of dimensionless frequency, xi
k: int
   number of amplitude correction parameters desired
mean_amplitude_difference: numpy.ndarray, optional
   array of the means of a set of amplitude differences, by default, this should be as a function of dimensionless frequency, xi
   if not given, the means of the alpha distributions are set to zero
   default: None
prior: bilby.core.prior.PriorDict, optional
   bilby prior object; if given, dA priors will be added to this dictionary
   default: None
geometrized: bool, optional
   if True, will return dimensionless frequency nodes; if False, standard frequency nodes (Hz)
   default: True
xi_low: float, optional
   if geometrized is True; lower bound on the dimensionless frequency band
   default: 0.018
xi_high: float, optional
   if geometrized is True; upper bound on the dimensionless frequency band
   default: 1/pi (0.318...)
f_low: float, optional
   if geometrized is False; lower bound on the standard frequency band
   default: 20.0 Hz
f_high: float, optional
   if geometrized is False; upper bound on the standard frequency band
   default: 1024.0 Hz
      
Returns:
--------
frequency_nodes: numpy.ndarray
   array of frequency nodes
prior: bilby.core.prior.PriorDict
   bilby prior object containing the amplitude correction priors
