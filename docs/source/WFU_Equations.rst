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

Where :math:`h_{1}` and :math:`h_{2}` are sets of frequency domain gravitational wave strain, which are complex. Amplitude difference is a relative error between the two waveforms relative to :math:`h_{1}`.

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

Where :math:`t_{c}` is an overall shift in coalescence time and :math:`\phi_{c}` is an overall phase shift. We can find these values by fitting a line to :math:`\Delta\phi` weighted by the power spectral density (PSD) data. The PSD data tells us the variance of the signal at each frequency point. This gives us information on :math:`t_{c}`. The 'y-intercept' of this line is :math:`\phi_{c}`.

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

Where :math:`\mu_{IMR}` is an ``IMRPhenomPv2_NRTidalv2`` waveform and :math:`\mu_{EOB}` is an ``SEOBNRv4T_surrogate`` waveform. :math:`f_{COR}`, the discontinuity correction frequency, is defined the following way:

.. math::

    \begin{equation}
        \frac{\partial^{2}}{\partial{f}^{2}}\left(\Delta{A}_{\mu}(f_{COR};\theta)\right)=-10^{-5}
    \end{equation}

This definition comes from utilizing the fact that the slopes in the :math:`\Delta{A}_{\mu}` curve are very small up until the discontinuity.

Waveform Uncertainty
--------------------
Waveform uncertainties are the variabilities of the waveform's amplitude and phase at a given frequency; they are enveloped. We can find a model for waveform uncertainty by taking the standard deviation of many draws of waveform difference. We define amplitude uncertainty and phase uncertainty in this way:

.. math::

    \begin{equation}
        \delta{A}(f)=\sqrt{\frac{\sum_{i=1}^{N}\left(\Delta{A}_{i}(f)-\overline{\Delta{A}}(f)\right)}{N}}
    \end{equation}

.. math::

    \begin{equation}
        \delta\Phi(f)=\sqrt{\frac{\sum_{i=1}^{N}\left(\Delta\Phi_{i}(f)-\overline{\Delta\Phi}(f)\right)}{N}}
    \end{equation}

.. note::

    We will be using residual phase uncertainty, :math:`\Delta\Phi`, as our phase uncertainty from now on.




















