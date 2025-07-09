TotalMassConstraint
===================

.. code-block:: python

  WaveformUncertainty.prior.TotalMassConstraint(*,name,f_low,f_high,latex_label=r'$M$',boundary=None,
                                                unit=r'$\mathrm{M}_\odot$',xi_low=0.018,xi_high=0.318)

Generates a bilby prior that constrains the total mass to ensure that the waveform correction is always within the frequency band.

.. math::

  \frac{c^3\xi_\mathrm{high}}{Gf_\mathrm{high}}<M<\frac{c^3\xi_\mathrm{low}}{Gf_\mathrm{low}}

Parameters:
-----------
name: string
  name of prior
f_low: float
  lower bound on frequency band in Hz
f_high: float
  upper bound on the frequency band in Hz
latex_label: string, optional, (r'$M$')
  label for the parameter in LaTeX
boundary: string, optional, (None)
  boundary condition for the prior
unit: string, optional, (r'$\mathrm{M}_\odot$')
  label for the unit of the parameter; default is solar mass
xi_low: float, optional, (0.018)
  lower bound on the dimensionless frequency band
xi_high: float, optional, (1/pi, 0.318...)
  upper bound on the dimensionless frequency band

Returns:
--------
total_mass_prior: bilby.core.prior.base.Constraint
  bilby constraint prior object for the total mass
