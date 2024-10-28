WaveformUncertainty.Q_factor
============================

.. code-block:: python

  WaveformUncertainty.Q_factor(WFU_result,NHP_result,injection)

Calculates the quality factor of a correction by comparing the posteriors of the waveform uncertainty corrected parameter estimation run and the null-hypothesis (uncorrection) parameter estimation run.

.. math::
  \epsilon_i=\frac{1}{\sqrt{n}}\sqrt{\sum_p\left(\frac{\theta_i^{\mathrm{ml}}[p]-\theta^{\mathrm{inj}}[p]}{\theta^{\mathrm{inj}}[p]}\right)^2}

.. math::
  Q=\frac{\epsilon_{\varnothing}-\epsilon_{\mathcal{C}}}{\epsilon_{\varnothing}}\times 100\%

Parameters:
-----------
WFU_result: bilby.core.result.Result
  bilby result object; output of the waveform uncertainty corrected parameter estimation run
NHP_result: bilby.core.result.Result
  bilby result object; output of the null-hypothesis parameter estimation run
injection: dictionary
  dictionary of the injection parameters given to the samplers before each parameter estimation run

Returns:
--------
Q_factor: float [%]
  quality factor
