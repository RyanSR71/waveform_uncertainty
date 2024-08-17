WaveformUncertainty.Q_factor
============================

.. code-block:: python

  WaveformUncertainty.Q_factor(WFU_result,NHP_result,injection)

Calculates the quality factor of a correction by comparing the posteriors of the waveform uncertainty corrected parameter estimation run and the null-hypothesis (uncorrection) parameter estimation run.

Parameters:
-----------
WFU_result: bilby.core.result.Result
  bilby result object; output of the waveform uncertainty corrected parameter estimation run
NHP_result: bilby.core.result.Result
  bilby result object; output of the null-hypothesis parameter estimation run
injection: dictionary
  dictionary of the injection parameters given to the samplers before each parameter estimation run
