WFU Equations
=============

Frequency Domain Waveform Differences
--------------------
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
