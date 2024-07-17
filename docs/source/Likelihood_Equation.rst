Likelihood Derivation
=====================
Here we show how the likelihood equation for sampling over waveform uncertainty was derived.

Following `Payne et al Phys. Rev. D 102, 122004 (2020) <https://arxiv.org/abs/2009.10193>`_ and `Bilby's Documentation <https://lscsoft.docs.ligo.org/bilby/likelihood.html#the-simplest-likelihood>`_, we start from a standard likelihood equation that assumes no uncertainties:

.. math::

  \begin{equation}
      \mathcal{L}_{\varnothing}(h|\theta)=\prod_{j}\frac{1}{2\pi{P(f_{j})}}\mathrm{exp}\left(-2\Delta{f}\frac{|h(f_{j})-\mu(f_{j};\theta)|^{2}}{P(f_{j})}\right),
  \end{equation}

where :math:`h` is frequency domain gravitational wave strain, :math:`\theta` is a set of source parameters for the waveform approximants, :math:`j` is an index corresponding to frequency bins, :math:`\Delta{f}` is the distance between frequency bins, :math:`P` is power spectral density data, and :math:`\mu` is a frequency domain waveform model. To sample over waveform uncertainty, we will need to add waveform differences to the model. To start, we rewrite the waveform model in terms of its amplitude, :math:`A`, and phase, :math:`\phi`:

.. math::

  \begin{equation}
      \mu=Ae^{i\phi}.
  \end{equation}

Now we can modify the model's amplitude and phase. We will do this by adding small deviations to both:

.. math::

  \begin{equation}
      \mu^{\prime}=(A+dA)e^{i(\phi+d\phi)}.
  \end{equation}

Now we consider what :math:`dA` and :math:`d\phi` are in terms of our waveform differences, :math:`\Delta{A}` and :math:`\Delta\Phi`. We look at the definitions of amplitude and phase difference (see `Equations and Notation <https://waveformuncertainty.readthedocs.io/en/latest/WFU_Equations.html>`_ for details):

.. note::

  We always use residual phase difference, so it will not be explicitly stated. If we use raw phase difference, it will be make clear.

.. math::
    
    \begin{equation}
        \Delta{A}=\frac{|h_{1}|-|h_{2}|}{|h_{1}|},
    \end{equation} 

.. math::

    \begin{equation}
        \Delta\Phi=\mathrm{tan}^{-1}\left(\frac{\mathrm{Im}\left[h_{1}\right]}{\mathrm{Re}\left[h_{1}\right]}\right)-\mathrm{tan}^{-1}\left(\frac{\mathrm{Im}\left[h_{2}\right]}{\mathrm{Re}\left[h_{2}\right]}\right)-(2\pi{t_{c}}f+\phi_{c}).
    \end{equation}

We see that our phase difference, :math:`\Delta\Phi`, is itself a phase, and can therefore be directly added to our model's phase. In other words: :math:`d\phi=\Delta\Phi`:

.. math::

  \begin{equation}
      \mu^{\prime}=(A+dA)e^{i(\phi+\Delta\Phi)}.
  \end{equation}

Amplitude difference, :math:`\Delta{A}`, however, cannot be directly added to our model's amplitude. Amplitude difference is a relative error, so to get it in a form that we can use, we need to multiply it by the model's amplitude, or :math:`dA=A\Delta{A}`:

.. math::

  \begin{equation}
      \mu^{\prime}=(A+A\Delta{A})e^{i(\phi+\Delta\Phi)}.
  \end{equation}

Now we can rearrange this equation using basic algebra; we will pull out the common factor, A, and split the exponential term into two exponentials:

.. math::

  \begin{equation}
      \mu^{\prime}=A(1+\Delta{A})e^{i\phi}e^{i\Delta\Phi}.
  \end{equation}

Using the definition of our model, :math:`\mu=Ae^{i\phi}`, we can rewrite our new model in terms of :math:`\mu`:

.. math::

  \begin{equation}
      \mu^{\prime}=\mu(1+\Delta{A})e^{i\Delta\Phi}=\mu(1+\Delta{A})\mathrm{exp}[i\Delta\Phi].
  \end{equation}

.. note:: 

  This final expression we get matches Equation 6 in `Jocelyn Read 2023 Class. Quantum Grav. 40 135002 <https://arxiv.org/abs/2301.06630v2>`_.

Using our new waveform model, :math:`\mu^{\prime}`, in the likelihood function, we arrive at the following equation:

.. math::

  \begin{equation}
      \mathcal{L}(h|\theta)=\prod_{j}\frac{1}{2\pi{P(f_{j})}}\mathrm{exp}\left(-2\Delta{f}\frac{|h(f_{j})-\mu(f_{j};\theta)(1+\Delta{A})\mathrm{exp}\left[i\Delta\Phi\right]|^{2}}{P(f_{j})}\right).
  \end{equation}

To sample these waveform differences, we need to express :math:`\Delta{A}` and :math:`\Delta\Phi` in terms of a small number of parameters. The simplest way to do this is to use cubic splines. Cubic splines take points, or nodes, and fill the space between them with cubic functions. We redefine our waveform differences to be spline functions:

.. math:: 

  \begin{equation}
      \Delta{A}\rightarrow\Delta{A}(f;\{f_{n},\alpha_{n}\}),
  \end{equation}

.. math:: 

  \begin{equation}
      \Delta\Phi\rightarrow\Delta\Phi(f;\{f_{n},\beta_{n}\}),
  \end{equation}

where the :math:`\alpha` and :math:`\beta` parameters can be varied by bilby's samplers and :math:`f_{n}` are fixed frequency nodes.

.. note::

  Cubic functions are chosen to ensure that the whole function is continuous and that the function's first derivative is also continuous, which makes the function smooth. It takes four degrees of freedom to connect two points this way, which is the number of degrees of freedom that cubic functions have.

We now need to know what prior distribution we are going to use to draw the :math:`\alpha` and :math:`\beta` parameters from. For this, we choose a Gaussian, or normal, distribution around zero with our waveform uncertainty model as the standard deviation:

.. math::

    \begin{equation}
        \alpha_{n}\sim\mathcal{N}(0,\delta{A}_{\mu}(f_{n})),
    \end{equation}

.. math::

    \begin{equation}
        \beta_{n}\sim\mathcal{N}(0,\delta\Phi_{\mu}(f_{n})).
    \end{equation}

Plugging these spline functions into the likelihood function gives the final form of the likelihood equation:

.. math::

    \small \begin{equation}
        \mathcal{L}(h|\theta,\alpha,\beta)=\prod_{j}\frac{1}{2\pi{P(f_{j})}}\mathrm{exp}\left(-2\Delta{f}\frac{|h(f_{j})-\mu(f_{j};\theta)\left(1+\Delta{A}_{\delta}(f_{j};\{f_{n},\alpha_{n}\})\right)\mathrm{exp}\left[i\Delta\Phi_{\delta}(f_{j};\{f_{n},\beta_{n}\})\right]|^{2}}{P(f_{j})}\right).
    \end{equation}

.. note::

  We give :math:`\Delta{A}` and :math:`\Delta\Phi` the subscript, :math:`\delta`, to denote that these waveform difference models are drawn from waveform uncertainties :math:`\delta{A}` and :math:`\delta\Phi`.









