Likelihood Equation
===================
Here we show how the likelihood equation for sampling over waveform uncertainty was derived.

Following `(Payne et al Phys. Rev. D 102, 122004 (2020)) <https://arxiv.org/abs/2009.10193>`_ and `bilby documentation <https://lscsoft.docs.ligo.org/bilby/likelihood.html#the-simplest-likelihood>`_, we start from a standard likelihood equation that assumes no uncertainties:

.. math::

  \begin{equation}
      \mathcal{L}_{\varnothing}(h|\theta)=\prod_{j}\frac{1}{2\pi{P(f_{j})}}\mathrm{exp}\left(-2\Delta{f}\frac{|h(f_{j})-\mu(f_{j};\theta)|^{2}}{P(f_{j})}\right)
  \end{equation}
