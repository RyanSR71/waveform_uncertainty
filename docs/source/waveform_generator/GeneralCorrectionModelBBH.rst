GeneralCorrectionModelBBH
=========================

.. code-block:: python

   def GWCorrect.waveform_generator.GeneralCorrectionModelBBH(parameters,**kwargs)

Modified lal binary black hole frequency domain source model to include waveform error corrections in the strain calculation.

.. math::

   \mu_\mathrm{GC}(f;\vartheta,\mathrm{A},\Phi)=\mu(f;\vartheta)(1+\Delta\mathcal{A}_\mathrm{SI}(f;\{\mathrm{f}_k,\tilde\alpha_k\sigma_{\alpha,k}\}))\exp(i\Delta\phi_\mathrm{SI}(f;\{\mathrm{f}_k,\tilde\varphi_k\sigma_{\varphi,k}\}))

Make sure to include "xi_0", "delta_xi_tilde", and all necessary "dA" and "dphi" parameters in the prior.

New Parameters:
---------------
xi_0: float
   position of the zeroth frequency node in dimensionless frequency
delta_xi_tilde: float
   parameter defining the position of the last frequency node in dimensionless frequency
dAs: list
   list of dA parameters
dphis: list
   list of dphi parameters
sigma_dA_spline: scipy.interpolate.CubicSpline, optional, (None)
   spline object from which to calculate the sigma values for the amplitude error correction
sigma_dphi_spline: scipy.interpolate.CubicSpline, optional, (None)
   spline object from which to calculate the sigma values for the phase error correction
xi_high: float, optional, (1/pi)
   absolute upper bound on dimensionless frequency

Returns:
--------
model_strain: dict
   plus and cross polarizations of the frequency domain strain
