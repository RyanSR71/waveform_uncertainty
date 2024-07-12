Likelihood Equation
===================
Here we show how the likelihood equation for sampling over waveform uncertainty was derived.

Following `(Payne et al Phys. Rev. D 102, 122004 (2020)) <https://arxiv.org/abs/2009.10193>`_ and `bilby documentation <https://lscsoft.docs.ligo.org/bilby/likelihood.html#the-simplest-likelihood>`_, we start from a standard likelihood equation that assumes no uncertainties:

.. math::

  \begin{equation}
      \mathcal{L}_{\varnothing}(h|\theta)=\prod_{j}\frac{1}{2\pi{P(f_{j})}}\mathrm{exp}\left(-2\Delta{f}\frac{|h(f_{j})-\mu(f_{j};\theta)|^{2}}{P(f_{j})}\right)
  \end{equation}

where h is our frequency domain stain and :math:`\mu(f;\theta)` is our waveform model. To sample over waveform uncertainty, we will need to add waveform differences to the model. To start, let's rewrite the waveform model in terms of its amplitude, A, and phase, :math:`\phi`:

.. math::

  \begin{equation}
      \mu=Ae^{i\phi}
  \end{equation}

Now we can modify the model's amplitude and phase. We will do this by adding small deviations to both:

.. math::

  \begin{equation}
      \mu^{\prime}=(A+dA)e^{i(\phi+d\phi)}
  \end{equation}

What are dA and :math:`d\phi` in terms of our waveform differences, :math:`\Delta{A}` and :math:`\Delta\Phi`? Let us look at the definitions of amplitude and phase difference:

.. note::

  We always use residual phase difference, so it will not be explicitly stated. If we use raw phase uncertainty, it will be make clear.

.. math::
    
    \begin{equation}
        \Delta{A}=\frac{|h_{1}|-|h_{2}|}{|h_{1}|}
    \end{equation} 

.. math::

    \begin{equation}
        \Delta\Phi=\mathrm{tan}^{-1}\left(\frac{\mathrm{Im}\left[h_{1}\right]}{\mathrm{Re}\left[h_{1}\right]}\right)-\mathrm{tan}^{-1}\left(\frac{\mathrm{Im}\left[h_{2}\right]}{\mathrm{Re}\left[h_{2}\right]}\right)-(2\pi{t_{c}}f+\phi_{c})
    \end{equation}

We see that our phase difference is itself a phase, and can therefore be directly added to our model's phase. In other words: :math:`d\phi=\Delta\Phi`:

.. math::

  \begin{equation}
      \mu^{\prime}=(A+dA)e^{i(\phi+\Delta\Phi)}
  \end{equation}








