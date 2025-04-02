WaveformUncertainty.dA_prior
=============================

.. code-block:: python

   WaveformUncertainty.dA_prior(amplitude_uncertainty,k,mean_amplitude_difference=None,prior=None,
                                geometrized=True,xi_low=0.018,xi_high=0.318,f_low=20.0,f_high=1024.0)

Automatically generates a bilby prior object containing Gaussian priors for each dA parameter.

.. math::

   \Pi(\alpha_k)=\mathcal{N}(0,\delta\mathcal{A}_\mu(f_k))

Parameters:
-----------
amplitude_uncertainty: numpy.ndarray
   array of standard deviation of a set of amplitude differences; by default, this should be as a function of dimensionless frequency, xi
k: int
   number of amplitude correction parameters desired
mean_amplitude_difference: numpy.ndarray, optional, (None)
   array of the means of a set of amplitude differences, by default, this should be as a function of dimensionless frequency, xi
   if not given, the means of the dA distributions are set to zero
prior: bilby.core.prior.PriorDict, optional, (None)
   bilby prior object; if given, dA priors will be added to this dictionary
geometrized: bool, optional, (True)
   if True, will return dimensionless frequency nodes; if False, standard frequency nodes (Hz)
xi_low: float, optional, (0.018)
   if geometrized is True; lower bound on the dimensionless frequency band
xi_high: float, optional, (1/pi, 0.318...)
   if geometrized is True; upper bound on the dimensionless frequency band
f_low: float, optional, (20.0 Hz)
   if geometrized is False; lower bound on the standard frequency band
f_high: float, optional, (1024.0 Hz)
   if geometrized is False; upper bound on the standard frequency band
      
Returns:
--------
frequency_nodes: numpy.ndarray
   array of frequency nodes
prior: bilby.core.prior.PriorDict
   bilby prior object containing the amplitude correction priors
