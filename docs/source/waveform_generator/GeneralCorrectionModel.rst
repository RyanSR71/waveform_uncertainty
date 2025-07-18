GeneralCorrectionModel
======================

.. code-block:: python

   class WaveformUncertainty.waveform_generator.GeneralCorrectionModel(duration=None,sampling_frequency=None,start_time=0,
                                                                       frequency_domain_source_model=None,
                                                                       time_domain_source_model=None,
                                                                       parameters=None,parameter_conversion=None,
                                                                       correction_arguments=None,waveform_arguments=None)

Bases: ``object``

Modified WaveformGenerator object from `bilby.gw.WaveformGenerator <https://lscsoft.docs.ligo.org/bilby/api/bilby.gw.waveform_generator.WaveformGenerator.html#bilby.gw.waveform_generator.WaveformGenerator>`_ to include waveform uncertainty corrections in the strain calculation.

.. math::

   \mu_\mathrm{GC}(f;\vartheta,\mathrm{A},\Phi)=\mu(f;\vartheta)(1+\Delta\mathcal{A}_\mathrm{SI}(f;\{\mathrm{f}_k,\tilde\alpha_k\sigma_{\alpha,k}\}))\exp(i\Delta\phi_\mathrm{SI}(f;\{\mathrm{f}_k,\tilde\varphi_k\sigma_{\varphi,k}\}))

Make sure to include "xi_0", "delta_xi_tilde", and all necessary "dA" and "dphi" parameters in the prior. To change smoothing parameter, gamma, include it in the prior.

.. code-block:: python

   __init__(duration=None,sampling_frequency=None,start_time=0,
                                                  frequency_domain_source_model=None,
                                                  time_domain_source_model=None,
                                                  parameters=None,parameter_conversion=None,
                                                  correction_arguments=None,waveform_arguments=None)

New Parameters:
---------------
correction_arguments: dict, optional
  dictionary containing arguments for the waveform correction
  default: None
  
  contents:
      correct_amplitude: bool
          whether or not to attempt an amplitude correction
          include all dA parameters in the prior if True
      sigma_dA: numpy.ndarry
          if correct_amplitude is True, this is the array of standard deviations for the dA priors
      correct_phase: bool
          whether or not to attempt a phase correction
          include all dphi parameters in the prior if True
      sigma_dphi: numpy.ndarry
          if correct_phase is True, this is the array of standard deviations for the dphi priors
      nodes: int
          number of frequency nodes desired
      xi_high: float, optional
          absolute lower bound on dimensionless frequency, xi
          default: 1/pi
