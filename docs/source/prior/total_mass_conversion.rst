total_mass_conversion
=====================

.. code-block:: python

  GWCorrect.prior.total_mass_conversion(parameters)

Conversion function between any mass parameters and the total system mass. Used with the bilby.core.prior.PriorDict() function

Parameters:
-----------
parameters: dictionary
  dictionary with binary black hole parameters

Returns:
--------
parameters: dictionary
  input dictionary, but with the total system mass added
