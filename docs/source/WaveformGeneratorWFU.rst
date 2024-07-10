WaveformUncertainty.WaveformGeneratorWFU
========================================

.. code-block:: python

   WaveformUncertainty.WaveformGeneratorWFU(duration=None,sampling_frequency=None,start_time=0,
                                            frequency_domain_source_model=None,
                                            time_domain_source_model=None,parameters=None,
                                            parameter_conversion=None,waveform_arguments=None,
                                            waveform_uncertainty_nodes=None,dA_sampling=False,
                                            dphi_sampling=False)

Bases: ``object``

Modified WaveformGenerator object from bilby.gw to include waveform uncertainty corrections in the strain calculation;
To sample waveform uncertainty, include all relevant "alpha" and "beta" parameters in the prior.
See bilby's documentation for 'bilby.gw.WaveformGenerator'_.

.. _bilby.gw.WaveformGenerator: https://lscsoft.docs.ligo.org/bilby/api/bilby.gw.waveform_generator.WaveformGenerator.html#bilby.gw.waveform_generator.WaveformGenerator

.. note::

  Make sure the number of alphas, betas, and waveform_uncertainty_nodes are the same!

.. code-block:: python

   __init__(duration=None,sampling_frequency=None,start_time=0,frequency_domain_source_model=None,
                                            time_domain_source_model=None,parameters=None,
                                            parameter_conversion=None,waveform_arguments=None,
                                            waveform_uncertainty_nodes=None,dA_sampling=False,
                                            dphi_sampling=False)

New Parameters:
---------------
waveform_uncertainty_nodes: numpy.ndarray, optional, (None)
    array of frequency nodes to be used in generating the dA and dphi splines
dA_sampling: bool, optional, (None)
    if True, the waveform generator will attempt to pull alpha parameters from the parameter dictionary (either an injection or the prior)
dphi_sampling: bool, optional, (None)
    if True, the waveform generator will attempt to pull beta parameters from the parameter dictionary (either an injection or the prior)
