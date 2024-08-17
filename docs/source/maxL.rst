WaveformUncertainty.maxL
========================

.. code-block:: python

  WaveformUncertainty.maxL(result)

Finds the set of parameters in a parameter estimation posterior that together yield the highest likelihood.

Parameters:
-----------
result: bilby.core.result.Result
  bilby result object; output from a parameter estimation run

Returns:
--------
maxL_dict: dictionary
  dictionary containing all of the parameters in the posterior and their most likely values
