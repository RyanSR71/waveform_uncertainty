WaveformUncertainty.WaveformGeneratorWFU
========================================

.. code-block:: python

   class WaveformUncertainty.WaveformGeneratorWFU(duration=None,sampling_frequency=None,start_time=0,
                                                  frequency_domain_source_model=None,
                                                  time_domain_source_model=None,parameters=None,
                                                  parameter_conversion=None,waveform_arguments=None,
                                                  frequency_nodes=None,correct_amplitude=False,
                                                  correct_phase=False,geometrized=True)

Bases: ``object``

Modified WaveformGenerator object from `bilby.gw.WaveformGenerator <https://lscsoft.docs.ligo.org/bilby/api/bilby.gw.waveform_generator.WaveformGenerator.html#bilby.gw.waveform_generator.WaveformGenerator>`_ to include waveform uncertainty corrections in the strain calculation.

.. math::

   \mu_\mathcal{C}(f;\vartheta,\alpha,\varphi)=\mu(f;\vartheta)(1+\Delta\mathcal{A}_s(f;\{f_k,\alpha_k\}))\exp[i\Delta\phi_s(f;\{f_k,\varphi_k\})]

To sample waveform uncertainty, include all necessary "dA" and "dphi" parameters in the prior.

.. code-block:: python

   __init__(duration=None,sampling_frequency=None,start_time=0,
                                                  frequency_domain_source_model=None,
                                                  time_domain_source_model=None,parameters=None,
                                                  parameter_conversion=None,waveform_arguments=None,
                                                  frequency_nodes=None,correct_amplitude=False,
                                                  correct_phase=False,geometrized=True)

New Parameters:
---------------
frequency_nodes: numpy.ndarray, optional, (None)
   array of frequency nodes to be used in generating the dA and dphi splines
correct_amplitude: bool, optional, (False)
   toggle for the amplitude correction
correct_phase: bool, optional, (False)
   toggle for the phase correction
geometrized: bool, optional, (True)
   whether or not the frequency nodes are in dimensionless frequency or not
