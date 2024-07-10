WaveformUncertainty.WFU_prior
=============================

.. code-block:: python

   WaveformUncertainty.WFU_prior(mean_amplitude_difference,amplitude_uncertainty,
                                 mean_phase_difference,phase_uncertainty,frequency_grid,
                                 nnodes,prior=None,spacing='linear')

Automatically generates a bilby prior object containing Gaussian waveform uncertainty parameter priors (alphas and betas). If given a pre-existing prior object, the waveform uncertainty parameters will be added to it

Parameters:
-----------
mean_amplitude_difference: numpy.ndarray
    array of mean amplitude difference values
amplitude_uncertainty: numpy.ndarray
    array of the amplitude uncertainty; defined as the standard deviation of amplitude differences
mean_phase_difference: numpy.ndarray
    array of mean phase difference values
phase_uncertainty: numpy.ndarray
    array of the phase uncertainty; defined as the standard deviation of phase differences
frequency_grid: numpy.ndarray
    frequency grid corresponding to the waveform uncertainties
nnodes: int
    number of frequency nodes desired
prior: bilby.core.prior.PriorDict, optional, (None)
    if given, the output prior will simply be added to this
spacing: string, optional, ('linear')
    dictates the type of progression the frequency nodes will have; {'linear', 'geometric'}
      
Returns:
--------
prior: bilby.core.prior.PriorDict
    prior containing the waveform uncertainty parameters (alphas and betas)
frequency_nodes: numpy.ndarray
    frequency nodes used by __WaveformGeneratorWFU() to generate waveform difference splines
