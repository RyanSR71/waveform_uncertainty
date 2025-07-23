binary_black_hole_correction_conversion
=======================================

.. code-block:: python

   def GWCorrect.waveform_generator.binary_black_hole_correction_conversion(parameters)

Converts set of parameters into the parameters needed by the source model. 

Parameters:
---------------
parameters: dict
    set of binary black hole parameters and correction parameters (xi_0, delta_xi_tilde, dA_0, dphi_0, etc.)

Returns:
--------
converted_parameters: dict
    input parameters, but with added parameters needed by the source model, such as 'dAs' and 'dphis'
added_keys: list
    list of parameter names that were added to the output parameter dictionary
