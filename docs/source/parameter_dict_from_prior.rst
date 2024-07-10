WaveformUncertainty.parameter_dict_from_prior
=============================================

.. code-block:: python

   WaveformUncertainty.parameter_dict_from_prior(prior,nsamples)

Takes a bilby prior object and constructs a dictionary of random samples.

Parameters:
----------
prior: bilby.core.prior.dict.PriorDict
    prior object
nsamples: int
    number of samples desired
      
Returns:
-------
parameter_dict: dictionary
    dictionary of nsamples random samples
