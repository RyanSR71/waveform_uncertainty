Equations and Notation
======================
This work is in collaboration with and derived from past work by Dr. Jocelyn Read. See `Jocelyn Read 2023 Class. Quantum Grav. 40 135002 <https://arxiv.org/abs/2301.06630v2>`_.

Frequency Domain Waveform Differences
-------------------------------------
Amplitude difference is defined generally as:

.. math::
    
    \begin{equation}
        \Delta{A}=\frac{|h_{1}|-|h_{2}|}{|h_{1}|}
    \end{equation} 

where :math:`h_{1}` and :math:`h_{2}` are sets of frequency domain gravitational wave strain, which are complex. Amplitude difference is a relative error between the two waveforms relative to :math:`h_{1}`.

Raw phase difference is then defined as:

.. math::

    \begin{equation}
        \Delta\phi=\mathrm{tan}^{-1}\left(\frac{\mathrm{Im}\left[h_{1}\right]}{\mathrm{Re}\left[h_{1}\right]}\right)-\mathrm{tan}^{-1}\left(\frac{\mathrm{Im}\left[h_{2}\right]}{\mathrm{Re}\left[h_{2}\right]}\right)
    \end{equation} 
    
Raw phase difference contains overall phase and time shifts which need to be removed for analysis. To do this, we find residual phase difference:

.. math::

    \begin{equation}
        \Delta\Phi=\Delta\phi-(2\pi{t_{c}}f+\phi_{c})
    \end{equation}

where :math:`t_{c}` is an overall shift in coalescence time and :math:`\phi_{c}` is an overall phase shift. We can find these values by fitting a line to :math:`\Delta\phi` weighted by the power spectral density (PSD) data. The PSD data tells us the variance of the signal at each frequency point. This gives us information on :math:`t_{c}`. The 'y-intercept' of this line is :math:`\phi_{c}`. We then subtract this line from raw phase difference to get residual phase difference.

Waveform Model Differences
--------------------------
When finding the waveform differences between two waveform models, :math:`\mu`, we need to add some slight devitations from the general forms listed above. If the amplitude of the first waveform, :math:`|h_{1}|`, went to zero faster than the second waveform, :math:`|h_{2}|`, then we would get a discontinuity. This is shown with the following limit:

.. math::

    \begin{equation}
        \lim_{|h_{1}|\to{0}}\left(\frac{|h_{1}|-|h_{2}|}{|h_{1}|}\right)=-\infty
    \end{equation}

This results in the :math:`\Delta{A}` curve abruptly going down to negative infinity. When dealing with waveform approximants, especially those that handle tidal defomabilities, this discontinuity is extremely common. To deal with these discontinuities, we simply cut off the curve at the discontinuity and hold it constant afterwards. We then define the following waveform model differences:

.. math::

    \begin{equation}
        \Delta{A}_{\mu}(f;\theta)= \begin{cases} 
          \frac{|\mu_{IMR}(f;\theta)|-|\mu_{EOB}(f;\theta)|}{|\mu_{IMR}(f;\theta)|} & f \leq f_{COR} \\
          \Delta{A}_{\mu}(f_{COR};\theta) & f > f_{COR} 
       \end{cases}
    \end{equation}

.. math::

    \begin{equation}
        \Delta\phi_{\mu}(f;\theta)= \begin{cases} 
          \mathrm{tan}^{-1}\left(\frac{\mathrm{Im}[\mu_{IMR}(f;\theta)]}{\mathrm{Re}[\mu_{IMR}(f;\theta)]}\right)-\mathrm{tan}^{-1}\left(\frac{\mathrm{Im}[\mu_{EOB}(f;\theta)]}{\mathrm{Re}[\mu_{EOB}(f;\theta)]}\right) & f \leq f_{COR} \\
          \Delta\phi_{\mu}(f_{COR};\theta) & f > f_{COR} 
       \end{cases}
    \end{equation}

.. math::

    \begin{equation}
        \Delta\Phi_{\mu}(f;\theta)= \begin{cases} 
          \Delta\phi_{\mu}(f;\theta)-(2\pi{t}_{0}{f}+\phi_{0}) & f \leq f_{COR} \\
          \Delta\Phi_{\mu}(f_{COR};\theta) & f > f_{COR} 
       \end{cases}
    \end{equation}

where :math:`\mu_{IMR}` is an ``IMRPhenomPv2_NRTidalv2`` waveform and :math:`\mu_{EOB}` is an ``SEOBNRv4T_surrogate`` waveform. :math:`f_{COR}`, the discontinuity correction frequency, is defined the following way:

.. math::

    \begin{equation}
        \frac{\partial^{2}}{\partial{f}^{2}}\left(\Delta{A}_{\mu}(f_{COR};\theta)\right)=-10^{-5}
    \end{equation}

This definition comes from utilizing the fact that the slopes in the :math:`\Delta{A}_{\mu}` curve are very small up until the discontinuity.

.. note::

    The right hand side of the discontinuity correction frequency equation can be adjusted in this package's functions.

Waveform Uncertainty
--------------------
Waveform uncertainties are the variabilities of the waveform's amplitude and phase at a given frequency; they are enveloped. We can find a model for waveform uncertainty by taking the standard deviation of many draws of waveform difference. We define our model amplitude uncertainty and phase uncertainty in this way:

.. math::

    \begin{equation}
        \delta{A}_{\mu}(f)=\sqrt{\frac{\sum_{i=1}^{N}\left(\Delta{A}_{\mu}(f;\theta_{i})-\overline{\Delta{A}_{\mu}}(f)\right)}{N}}
    \end{equation}

.. math::

    \begin{equation}
        \delta\Phi_{\mu}(f)=\sqrt{\frac{\sum_{i=1}^{N}\left(\Delta\Phi_{\mu}(f;\theta_{i})-\overline{\Delta\Phi_{\mu}}(f)\right)}{N}}
    \end{equation}

.. note::

    We will be using residual phase uncertainty, :math:`\Delta\Phi`, as our phase uncertainty from now on.

Parameterization
----------------
Computationally, generating individual waveform differences is a simple and quick task. However, to generate waveform uncertainty, we need many draws of waveform differences; at least 1000 for a decent model. Generating this number of waveform differences can take a lot of time and is generally tedious to do every time we want waveform uncertainty. To solve this issue, we can parameterize each waveform difference curve and save the parameters in a file. That way, we can generate all of our draws of waveform differences once and can simply load in the data in seconds next time we need them. This is achieved using Chebyshev polynomial series:

.. math:: 

    \begin{equation}
        \Delta{A}_{\mu}(f;\theta)\approx\Delta{A}_{T}(f;a,f_{COR},\Delta{A}_{\mu}(f_{COR};\theta))= \begin{cases} 
          \sum_{i=0}^{N-1}a_{i}T_{i}(f) & f \leq f_{COR} \\
          \Delta{A}_{\mu}(f_{COR};\theta) & f > f_{COR} 
       \end{cases}
    \end{equation}

.. math::

    \begin{equation}
       \Delta\Phi_{\mu}(f;\theta)\approx\Delta\Phi_{T}(f;b,f_{COR},\Delta\Phi_{\mu}(f_{COR};\theta))= \begin{cases} 
          \sum_{i=0}^{N-1}b_{i}T_{i}(f) & f \leq f_{COR} \\
          \Delta\Phi_{\mu}(f_{COR};\theta) & f > f_{COR} 
       \end{cases}
    \end{equation}

where :math:`T_{n}` are Chebyshev polynomials of the first kind. We see that instead of trying to carry around waveform models, which do not have simple functional forms, we can carry around a handful of coefficients, discontinuity correction frequencies, and the values the waveform differences level off at. With these parameters, we can reconstruct the original waveform differences within 2% in :math:`\Delta{A}` and :math:`2^{\circ}` in :math:`\Delta\Phi`. 

.. note::

    The error margins on :math:`\Delta{A}_{T}` and :math:`\Delta\Phi_{T}` can be adjusted in this package's functions.

Likelihood
----------
The likelihood function we use to sample over waveform uncertainty is

.. math::

    \small \begin{equation}
        \mathcal{L}(h|\theta,\alpha,\beta)=\prod_{j}\frac{1}{2\pi{P(f_{j})}}\mathrm{exp}\left(-2\Delta{f}\frac{|h(f_{j})-\mu(f_{j};\theta)\left(1+\Delta{A}_{\delta}(f_{j};\{f_{n},\alpha_{n}\})\right)\mathrm{exp}\left[i\Delta\Phi_{\delta}(f_{j};\{f_{n},\beta_{n}\})\right]|^{2}}{P(f_{j})}\right)
    \end{equation}

where the :math:`\alpha` and :math:`\beta` parameters are spline parameters corresponding to frequency nodes :math:`f_{n}`. These parameters are defined as being draws from a normal distribution around zero with their standard deviations being our waveform uncertainties:

.. math::

    \begin{equation}
        \alpha_{n}\sim\mathcal{N}(0,\delta{A}_{\mu}(f_{n}))\hspace{1cm}\beta_{n}\sim\mathcal{N}(0,\delta\Phi_{\mu}(f_{n}))
    \end{equation}




















