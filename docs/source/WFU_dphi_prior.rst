WaveformUncertainty.dphi_prior
=============================

.. code-block:: python

   WaveformUncertainty.dphi_prior(phase_uncertainty,k,mean_phase_difference=None,prior=None,
                                  geometrized=True,xi_low=0.018,xi_high=0.318,f_low=20.0,f_high=1024.0)

Automatically generates a bilby prior object containing Gaussian priors for each dphi parameter.

.. math::

   \Pi(\varphi_k)=\mathcal{N}(0,\delta\phi_\mu(f_k))

Parameters:
-----------
phase_uncertainty: numpy.ndarray
   array of standard deviation of a set of phase differences; by default, this should be as a function of dimensionless frequency, xi
k: int
   number of phase correction parameters desired
mean_phase_difference: numpy.ndarray, optional, (None)
   array of the means of a set of phase differences, by default, this should be as a function of dimensionless frequency, xi
   if not given, the means of the dphi distributions are set to zero
prior: bilby.core.prior.PriorDict, optional, (None)
   bilby prior object; if given, dphi priors will be added to this dictionary
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
   bilby prior object containing the phase correction priors
