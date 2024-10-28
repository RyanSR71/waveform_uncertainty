WaveformUncertainty.maxL
========================

.. code-block:: python

  WaveformUncertainty.maxL(result)

Finds the set of parameters in a parameter estimation posterior that together yield the highest likelihood.

.. math::
  \mathcal{L}(h|\lambda^\mathrm{ml})=\mathrm{max}[\mathcal{L}(h|\lambda)]

Parameters:
-----------
result: bilby.core.result.Result
  bilby result object; output from a parameter estimation run

Returns:
--------
maxL_dict: dictionary
  dictionary of the maximum likelihood values of each of the injected parameters
