WaveformUncertainty.injection
=============================

.. code-block:: python

   WaveformUncertainty.injection(data,index=randint(),precession=False,tides=True)

Pulls a sample out of a parameter dictionary

Parameters:
-----------
data: dictionary
    Dictionary of parameter samples (from paramater_dict_from_prior())
index: int, optional, (randint())
    position within the dict to pull the sample
precession: bool, optional, (False)
    only True if waveform approximant supports precessing spins;
    if False, precessing spin parameters will be removed/replaced with non-precessing parameters
tides: bool, optional, (True)
    only True if waveform approximant supports tidal deformabilities;
    if False, tidal parameters will be removed
      
Returns:
--------
injection: dictionary
    dictionary of injection parameters
