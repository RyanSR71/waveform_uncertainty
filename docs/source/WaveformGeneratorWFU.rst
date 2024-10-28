WaveformUncertainty.WaveformGeneratorWFU
========================================

.. code-block:: python

   class WaveformUncertainty.WaveformGeneratorWFU(duration=None,sampling_frequency=None,start_time=0,
                                                  frequency_domain_source_model=None,
                                                  time_domain_source_model=None,parameters=None,
                                                  parameter_conversion=None,waveform_arguments=None,
                                                  frequency_nodes=None,indexes=None)

Bases: ``object``

Modified WaveformGenerator object from bilby.gw to include waveform uncertainty corrections in the strain calculation.

.. math::

   \mu_\mathcal{C}(f;\theta,\alpha,\varphi)=\mu(f;\theta)(1+\Delta\mathcal{A}_s(f;\{f_k,\alpha_k\}))\exp[i\Delta\phi_s(f;\{f_k,\varphi_k\})]

To sample waveform uncertainty, include all necessary "dA" and "dphi" parameters in the prior.

.. note::

   See bilby's documentation for `bilby.gw.WaveformGenerator <https://lscsoft.docs.ligo.org/bilby/api/bilby.gw.waveform_generator.WaveformGenerator.html#bilby.gw.waveform_generator.WaveformGenerator>`

.. code-block:: python

   __init__(duration=None,sampling_frequency=None,start_time=0,frequency_domain_source_model=None,
                                            time_domain_source_model=None,parameters=None,
                                            parameter_conversion=None,waveform_arguments=None,
                                            frequency_nodes=None,indexes=None)

New Parameters:
---------------
frequency_nodes: numpy.ndarray, optional, (None)
   array of frequency nodes to be used in generating the dA and dphi splines
indexes: numpy.ndarray, optional, (None)
   list of dA and/or dphi parameters present in the prior; i.e. if you have dphi_1-dphi_5, indexes=[1,2,3,4,5]
