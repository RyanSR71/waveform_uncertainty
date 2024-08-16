Equations and Notation
======================
This work is in collaboration with and derived from past work by Dr. Jocelyn Read. See `Jocelyn Read 2023 Class. Quantum Grav. 40 135002 <https://arxiv.org/abs/2301.06630v2>`_.

Introduction
------------
Compact Binary Inspiral Gravitional Waves, as the name suggests, are gravitational waves that result from the inspiral of two compact objects, such as black holes and neutron stars. As the two objects fall together, they orbit around each other faster and faster, sending out more intense gravitational waves. This reaches a peak amplitude at the moment the objects come together. We can detect these gravitational waves using large Michelson interferometers, such as LIGO Hanford and LIGO Livingston. The waveform we receive is in the form of gravitational wave strain, :math:`h`, which is complex and a function of time:

.. math::

    \begin{equation}
        h(t)=\mathcal{A}(t)\mathrm{exp}\left({i\psi(t)}\right),
    \end{equation}

where :math:`\mathcal{A}(t)` is the time domain amplitude and :math:`\psi(t)` is the time domain phase. 

Typically, we work with waveforms in the frequency domain, which is the Fourier transform of the time domain waveform:

.. math::

    \begin{equation}
        \tilde{h}(f)=\int_{-\infty}^{\infty}h(t)\mathrm{e}^{-i2\pi{ft}}dt.
    \end{equation}

In the frequency domain, we represent the waveform in terms of its amplitude, :math:`{A}`, and phase, :math:`\phi`:

.. math::

    \begin{equation}
        \tilde{h}(f)={A}(f)\mathrm{exp}\left({i\phi(f)}\right).
    \end{equation}

From the waveform itself, we can find the amplitude and phase with

.. math::

    \begin{equation}
        {A}(f)=\sqrt{\tilde{h}^{*}(f)\tilde{h}(f)}=|\tilde{h}(f)|,
    \end{equation}

and

.. math::

    \begin{equation}
        \phi(f)=\mathrm{tan}^{-1}\left(\frac{\mathrm{Im}[\tilde{h}(f)]}{\mathrm{Re}[\tilde{h}(f)]}\right).
    \end{equation}

Using these waveforms and our best understanding and models of general relativity, we can work backwards from the signal to infer the properties of the two objects that created the gravitational waves. We do this with parameter estimation, which uses Bayesian statistics alongside a waveform model, or waveform approximant, to show a probability distribution for the parameters in the system; such as source masses, spins, and tidal deformabilties; which is known as the posterior. The waveform approximants we use, such as ``IMRPhenomPv2_NRTidalv2`` and ``SEOBNRv4T_surrogate``, approximate the system well, but are not perfect. The error between the approximated waveform and the true waveform is waveform uncertainty, which is defined in the amplitude and phase of the waveform as :math:`\delta\mathcal{A}` and :math:`\delta\Phi` respectively. 

``WaveformUncertainty`` introduces a method of correcting for these uncertainties by adding waveform difference parameters to the waveform model, allowing the waveform differences to be sampled during a parameter estimation run. This process molds the approximated waveform to the true waveform, eliminating the error and returning a more accurate posterior. The following sections will introduce the necessary equations and notation behind the functions and usage of this package.

Frequency Domain Waveform Differences
-------------------------------------
Amplitude difference, :math:`\Delta\mathcal{A}`, is defined generally as

.. math::
    
    \begin{equation}
        \Delta{A}(f)=\frac{|\tilde{h}_{2}(f)|-|\tilde{h}_{1}(f)|}{|\tilde{h}_{1}(f)|},
    \end{equation} 

where :math:`h_{1}` and :math:`h_{2}` are sets of frequency domain gravitational wave strain, which are complex. Amplitude difference is a relative error between the two waveforms relative to :math:`h_{1}`.

Raw phase difference, :math:`\Delta\phi`, is then defined as

.. math::

    \begin{equation}
        \Delta\phi(f)=\mathrm{tan}^{-1}\left(\frac{\mathrm{Im}\left[\tilde{h}_{2}(f)\right]}{\mathrm{Re}\left[\tilde{h}_{2}(f)\right]}\right)-\mathrm{tan}^{-1}\left(\frac{\mathrm{Im}\left[\tilde{h}_{1}(f)\right]}{\mathrm{Re}\left[\tilde{h}_{1}(f)\right]}\right),
    \end{equation} 
    
Raw phase difference contains overall phase and time shifts which need to be removed for analysis. To do this, we find residual phase difference, :math:`\Delta\Phi`:

.. math::

    \begin{equation}
        \Delta\Phi(f)=\Delta\phi(f)-(2\pi{t_{c}}f+\phi_{c}),
    \end{equation}

where :math:`t_{c}` is coalescence time and :math:`\phi_{c}` is coalescence phase. We can find these values by fitting a line to :math:`\Delta\phi` weighted by :math:`\frac{S_{n}(f)}{{A}(f)^{2}}`, where :math:`S_{n}(f)` is power spectral density (PSD) data and :math:`{A}(f)` is the amplitude of the waveform. The PSD data tells us the variance of the signal at each frequency point. We then subtract this line from raw phase difference to get residual phase difference.

Waveform Model Differences
--------------------------
When finding the waveform differences between two waveform models, :math:`\mu`, we need to add some slight deviatations from the general forms listed above. If the amplitude of the first waveform, :math:`|h_{1}|`, went to zero faster than the second waveform, :math:`|h_{2}|`, then we would get a discontinuity. This is shown with the following limit:

.. math::

    \begin{equation}
        \lim_{|\tilde{h}_{1}(f)|\to{0}}\left(\frac{|\tilde{h}_{2}(f)|-|\tilde{h}_{1}(f)|}{|\tilde{h}_{1}(f)|}\right)=\infty.
    \end{equation}

This results in the :math:`\Delta{A}` values abruptly increasing to infinity. When using waveform approximants, especially those that handle tidal deformabilities, this discontinuity is extremely common. If :math:`|\tilde{h}_{2}(f)|` instead goes to zero faster than :math:`|\tilde{h}_{1}(f)|`, we observe the following:

.. math::

    \begin{equation}
        \lim_{|\tilde{h}_{2}(f)|\to{0}}\left(\frac{|\tilde{h}_{2}(f)|-|\tilde{h}_{1}(f)|}{|\tilde{h}_{1}(f)|}\right)=-1.
    \end{equation}

This results in the amplitude difference abruptly cutting to :math:`-1`. While this is not necessarily a discontinuity, it will be treated as such.

.. note::

    Both of these discontinuities are related to each other in that swapping :math:`h_{1}` and :math:`h_{2}` in the amplitude difference equation switches the nature of the discontinuity from one to the other.

To deal with these discontinuities, we cut off the waveform differences at the discontinuity and hold it constant afterwards. We then define the following waveform model differences:

.. math::

    \begin{equation}
        \Delta{A}_{\mu}(f;\theta)= \begin{cases} 
          \frac{|\mu_{2}(f;\theta)|-|\mu_{1}(f;\theta)|}{|\mu_{1}(f;\theta)|} & f \leq f_{COR} \\
          \Delta{A}_{\mu}(f_{COR};\theta) & f > f_{COR}
       \end{cases}\hspace{0.2cm},
    \end{equation}

.. math::

    \begin{equation}
        \Delta\phi_{\mu}(f;\theta)= \begin{cases} 
          \mathrm{tan}^{-1}\left(\frac{\mathrm{Im}[\mu_{2}(f;\theta)]}{\mathrm{Re}[\mu_{2}(f;\theta)]}\right)-\mathrm{tan}^{-1}\left(\frac{\mathrm{Im}[\mu_{1}(f;\theta)]}{\mathrm{Re}[\mu_{1}(f;\theta)]}\right) & f \leq f_{COR} \\
          \Delta\phi_{\mu}(f_{COR};\theta) & f > f_{COR}
       \end{cases}\hspace{0.2cm},
    \end{equation}

.. math::

    \begin{equation}
        \Delta\Phi_{\mu}(f;\theta)= \begin{cases} 
          \Delta\phi_{\mu}(f;\theta)-(2\pi{t}_{c}{f}+\phi_{c}) & f \leq f_{COR} \\
          \Delta\Phi_{\mu}(f_{COR};\theta) & f > f_{COR} 
       \end{cases}\hspace{0.2cm},
    \end{equation}

where :math:`\mu_{1}` and :math:`\mu_{2}` are waveform models and :math:`\theta` is a set of source parameters used by the models to generate the waveform. The discontinuity correction frequency, :math:`f_{COR}`, is the frequency value at which the discontinuity occurs. The model waveform differences are defined in such a way that adding them to :math:`\mu_{1}` will match it's amplitude and phase to :math:`\mu_{2}`. For that reason, we often call :math:`\mu_{2}` the reference waveform.

Waveform Uncertainty
--------------------
Waveform uncertainties are the variabilities of the waveform's amplitude and phase at a given frequency. We can find a model for waveform uncertainty by taking the standard deviation of many sample sets of waveform difference. We define our model amplitude uncertainty and phase uncertainty in this way:

.. math::

    \begin{equation}
        \delta{A}_{\mu}(f)=\sqrt{\frac{\sum_{i=1}^{N}\left(\Delta{A}_{\mu}(f;\theta_{i})-\overline{\Delta{A}_{\mu}}(f)\right)^{2}}{N}},
    \end{equation}

.. math::

    \begin{equation}
        \delta\Phi_{\mu}(f)=\sqrt{\frac{\sum_{i=1}^{N}\left(\Delta\Phi_{\mu}(f;\theta_{i})-\overline{\Delta\Phi_{\mu}}(f)\right)^{2}}{N}},
    \end{equation}

where :math:`\theta_{i}` is a set of source parameters, :math:`N` is the number of draws of waveform difference, and :math:`\overline{\Delta{A}_{\mu}}` and :math:`\overline{\Delta\Phi_{\mu}}` are the mean waveform differences.

.. note::

    We will be using residual phase uncertainty, :math:`\Delta\Phi`, as our phase uncertainty from now on.

The mean waveform difference are defined in amplitude and phase, :math:`\overline{\Delta{A}_{\mu}}` and :math:`\overline{\Delta\Phi_{\mu}}` respectively, and are defined as:

.. math::

    \begin{equation}
        \overline{\Delta{A}_{\mu}}(f)=\frac{\sum_{i=1}^{N}(\Delta{A}_{\mu}(f;\theta_{i}))}{N},
    \end{equation}

and

.. math::

    \begin{equation}
        \overline{\Delta\Phi_{\mu}}(f)=\frac{\sum_{i=1}^{N}(\Delta\Phi_{\mu}(f;\theta_{i}))}{N}.
    \end{equation} 

.. note::

    For both the waveform uncertainties (:math:`\delta{A}_{\mu}` and :math:`\delta\Phi_{\mu}`) and the mean waveform differences (:math:`\overline{\Delta{A}_{\mu}}` and :math:`\overline{\Delta\Phi_{\mu}}`), each draw has different source parameters, denoted by :math:`\theta_{i}`.

Likelihood and Sampling
-----------------------
Parameter estimation, in the context of gravitational waves, is a process that utilizes Bayes' Theorem and Bayesian statistics to infer the properties of the objects that created the gravitational waves. Given a waveform model and the gravitational wave data, a sampler, such as ``nestle`` or ``dynesty`` can choose random samples for each parameter in the system. This random draw is then put into the model, which is then compared to the gravitational wave data. This comparison is done using a likelihood function, which peaks when the model and the data match. Repeating this process many times maps out the likelihood for each parameter. 

The likelihood function we use to sample over waveform uncertainty is

.. math::

    \begin{equation}
        \mathcal{L}(h|\theta,\alpha,\varphi)=\prod_{j}\frac{1}{2\pi{S_{n}(f_{j})}}\mathrm{exp}\left(-2\Delta{f}\frac{|\tilde{h}(f_{j})-\mu(f_{j};\theta)\cdot\nu(f_{j};\alpha,\varphi)|^{2}}{S_{n}(f_{j})}\right),
    \end{equation}

where :math:`h` is frequency domain gravitational wave strain, :math:`\theta` is a set of source parameters for the waveform approximants, :math:`\alpha` and :math:`\varphi` parameters are spline parameters corresponding to frequency nodes :math:`f_{k}`, :math:`j` is an index corresponding to frequency bins, :math:`\Delta{f}` is the distance between frequency bins (step frequency), :math:`S_{n}` is power spectral density data, :math:`\mu` is a frequency domain waveform model, and :math:`\nu` is a function of waveform differences known as the model correction function. The model correction function serves to match the waveform model to the data by taking into account waveform uncertainty. It is defined as

.. math::

    \begin{equation}
        \nu(f;\alpha,\varphi)=(1+\Delta{A}_{\delta}(f;\{f_{k},\alpha_{k}\})\mathrm{exp}[i\Delta\Phi_{\delta}(f;\{f_{k},\varphi_{k}\})],
    \end{equation}

where :math:`\Delta{A}_{\delta}` is an amplitude difference function defined by waveform uncertainty, :math:`f_{k}` is a set of frequency nodes, :math:`\alpha` is a set of amplitude difference spline nodes, :math:`\Delta\Phi_{\delta}` is a phase difference function, and :math:`\varphi` is a set of phase difference spline nodes. Each :math:`\alpha` and :math:`\varphi` parameter is a draw from a Gaussian distribution. Their priors are defined as

.. math::

    \begin{equation}
        P(\alpha_{k})=\frac{(2\pi)^{-\frac{1}{2}}}{\delta{A}_{\mu}(f_{k})}\mathrm{exp}\left[-\frac{1}{2}\left(\frac{\Delta{A}-\overline{\Delta{A}_{\mu}}(f_{k})}{\delta{A}_{\mu}(f_{k})}\right)^{2}\right]
    \end{equation},

and

.. math::

    \begin{equation}
        P(\varphi_{k})=\frac{(2\pi)^{-\frac{1}{2}}}{\delta\Phi_{\mu}(f_{k})}\mathrm{exp}\left[-\frac{1}{2}\left(\frac{\Delta\Phi-\overline{\Delta\Phi_{\mu}}(f_{k})}{\delta\Phi_{\mu}(f_{k})}\right)^{2}\right]
    \end{equation},

where :math:`\delta{A}_{\mu}(f)` and :math:`\delta\Phi_{\mu}(f)` are amplitude and phase uncertainty respectively, :math:`\overline{\Delta{A}_{\mu}}(f)` and :math:`\overline{\Delta\Phi_{\mu}}(f)` are mean amplitude and phase differences respectively, and :math:`\Delta{A}_{\mu}(f;\theta)` and :math:`\Delta\Phi_{\mu}(f;\theta)` are model amplitude and phase difference respectively.

Parameterizing Waveform Differences
-----------------------------------
Computationally, generating individual waveform differences is a simple and quick task. However, to generate waveform uncertainty, we need many sets of waveform differences; at least 1000 for a decent model. Generating this number of waveform differences can take a lot of time and is generally tedious to do every time we want waveform uncertainty. To solve this issue, we can parameterize each waveform difference curve and save the parameters in a file. That way, we can generate all of our draws of waveform differences once and can simply load in the data in seconds next time we need them. This is achieved using Chebyshev polynomial series up to the discontinuity, as shown here:

.. math:: 

    \begin{equation}
        \Delta{A}_{\mu}(f;\theta)\approx\Delta{A}_{T}(f;a,f_{COR},\Delta{A}_{\mu}(f_{COR};\theta))= \begin{cases} 
          \sum_{i=0}^{N-1}a_{i}T_{i}(f) & f \leq f_{COR} \\
          \Delta{A}_{\mu}(f_{COR};\theta) & f > f_{COR}
       \end{cases}\hspace{0.2cm},
    \end{equation}

.. math::

    \begin{equation}
       \Delta\Phi_{\mu}(f;\theta)\approx\Delta\Phi_{T}(f;b,f_{COR},\Delta\Phi_{\mu}(f_{COR};\theta))= \begin{cases} 
          \sum_{i=0}^{N-1}b_{i}T_{i}(f) & f \leq f_{COR} \\
          \Delta\Phi_{\mu}(f_{COR};\theta) & f > f_{COR} 
       \end{cases}\hspace{0.2cm},
    \end{equation}

where :math:`T_{n}` are Chebyshev polynomials of the first kind. In a file, we store the Chebyshev coefficients, :math:`a` and :math:`b`; the discontinuity correction frequency, :math:`f_{COR}`; the values of the waveform differences at :math:`f_{COR}`, :math:`\Delta{A}_{\mu}(f_{COR};\theta)` and :math:`\Delta\Phi_{\mu}(f_{COR};\theta)`; and other parameters needed to store the data. With these parameters, we can reconstruct the original waveform differences within 2% in :math:`\Delta{A}` and :math:`2^{\circ}` in :math:`\Delta\Phi`. 

.. note::

    The error margins on :math:`\Delta{A}_{T}` and :math:`\Delta\Phi_{T}` can be adjusted in this package's functions. See ``max_ampltitude_error`` and ``max_phase_error`` in `WaveformUncertainty.parameterization <https://waveformuncertainty.readthedocs.io/en/latest/parameterization.html>`_.
