WaveformUncertainty.recovery_from_parameterization
==================================================

.. code-block:: python

   WaveformUncertainty.recovery_from_parameterization(identity,data)

Converts a parameterized set of waveform difference back into waveform difference arrays (See `WaveformUncertainty.parameterization <https://waveformuncertainty.readthedocs.io/en/latest/parameterization.html>`_)

Parameters:
-----------
identity: string
    specifies which waveform difference should be returned; {'amplitude_difference', 'phase_difference'}
data: numpy.ndarray
    one index of the output matrix from WaveformUncertainty.parameterization(); input WaveformUncertainty.parameterization()[index]
  
Returns:
--------
difference_array: numpy.ndarray
    array of the waveform difference converted from the parameterization; has the same shape as the frequency grid within the original matrix
