ppECorrectionModel
==================

.. code-block:: python

   GWCorrect.ppE.ppECorrectionModel(parameters,**kwargs)

Modified lal binary black hole frequency domain source model to include ppE correction in the strain calculation.

.. math::

   \mu_\mathrm{ppE}(f;\vartheta,\Phi)=\mu(f;\vartheta)\exp(i\Delta\phi_\mathrm{ppE}(f;\tilde\beta,\delta\tilde\epsilon,b)

Make sure to include "beta_tilde", "delta_epsilon_tilde", and "b" in the prior.

New Parameters:
---------------
beta_tilde: float
   rescaled phase correction amplitude parameter
delta_epsilon_tilde: float
   rescaled ringdown frequency correction parameter
b: float
   PN order parameter

Returns:
--------
model_strain: dict
   plus and cross polarizations of the frequency domain strain
