=================================
EIS and Equivalent Circuit Models
=================================

.. toctree::
   :maxdepth: 2
   :hidden:

Abstract
========

The simulations tools that are used to monitor battery performances include varying degree of physical
complexity. On one side, we have P2D models based on the Doyle-Fuller-Newman approach but it can be very
computationally demanding and with many parameters difficult to parameterize. On the other side, we have the
Equivalent Circuit Models (ECM) which are used for real-time Battery Management Systems (BMS) due to
their low computational cost.

In this work, we aim at bridging the two models, that is developing a methodology which allows us to
transfer the parameters between the models. We present an automated calibration workflow in BattMo to extract
2-RC Thevenin parameters from synthetic (P2D) and experimental Electrochemical Impedance Spectroscopy (EIS) data.

Objective
=========

Effective everyday operation of a battery requires the Battery Management
System (BMS) to continuously monitor its real-time state (State of Charge,
State of Health). Given the practical challenges of measuring this state via
chemical methods, modeling offers a useful alternative. Physics-based models
offer very good approximations, but need thousands of parameters. By
identifying the very few parameters of the Equivalent Circuit Model (ECM) in
real-time, valuable information about the system can be extracted. Because each
set of parameters maps to a specific physical state, the first approach could
be referencing a lookup table to link a particular state to a set of
parameters. Another one, more accurate, could be to train a neural network to
develop a strong state-prediction model.

Theory
======

The electrical activity of a battery can be modelled by an Open Circuit Voltage
(OCV), which is a perfect source of current. This is the simplest way to
simulate it. One can realize that it is necessary to go a bit deeper in the
analysis, with the addition of a resistor, which symbolizes the heat
dissipation. This model remains very simple. It is possible to add the
diffusion phenomena with a series of :math:`n` parallel resistor-capacitance
pairs. This is the :math:`n`-Thevenin model.

Example of a 2-Thevenin Model (`Sun, comprehensive review of ECM <https://www.mdpi.com/2079-9292/15/9/1968>`_)

It is a good compromise to have 2 RC networks, because the model remains simple
but allows two different time-scaled phenomena. However, another flaw of this
model is the symmetry of the two pairs, which makes identifiability much
harder.

The impedance of this circuit is

.. math::

   Z = R_0 + \frac{R_1}{1 + j R_1 C_1 \omega}
   + \frac{R_2}{1 + j R_2 C_2 \omega}.

Another way to model diffusion is with a Warburg component, which has the
impedance :math:`Z_W = \frac{1}{Y_0 \sqrt{j\omega}}` and adds to the Nyquist
diagram a 45° slope for low frequencies. A Warburg element is the equivalent of
an infinity of capacitance/resistor pairs.

Warburg circuit model (`Analog Devices <https://www.analog.com/media/en/reference-design-documentation/reference-designs/cn0510.pdf>`_, 2021)

Some researchers, as Santoni in his `paper <https://www.sciencedirect.com/science/article/pii/S2352152X2303788X>`_,
have tried some other similar circuits with more accurate tools, but further
from basic electrical elements. He replaced capacitors by constant phase
elements (CPE), adding an :math:`\alpha` parameter to quantify non-ideality.
The impedance is :math:`Z = \frac{1}{Q\left(j\omega\right)^{\alpha}}`.
:math:`\alpha = 1` would make it an ideal capacitor, while
:math:`\alpha = 0.5` corresponds to a Warburg element.

A 2nd-Order model with CPE used by Santoni

Physical Meaning of ECM components (`Analog Devices <https://www.analog.com/media/en/reference-design-documentation/reference-designs/cn0510.pdf>`_, 2021)
==========================================================================================================================================================

ECM are built after Occam's razor principles, where each component has to be
linked to a physical phenomenon of the battery (`Fletcher <https://link.springer.com/article/10.1007/s10008-013-2328-4>`_, 2013):

*The ohmic resistance* (:math:`R_\Omega`) represents the battery's internal
resistance (of the electrolyte and of the electrode interface for instance) and
accounts for heat generation during operation. This resistance naturally
increases as the battery ages. Experimentally, it can be determined by
observing the battery's voltage drop response to a load. On a Nyquist diagram,
it corresponds to the real part of the impedance where the imaginary component
is zero.

*The first RC network* is associated with the charge-transfer process at the
electrode-electrolyte interface. The resistance :math:`R_{\textrm{CT}}`
reflects the kinetic barrier electrons must overcome to move from the solid
electrode into the liquid electrolyte. The capacitance
:math:`C_{\textrm{DL}}` represents the double-layer capacitance, formed by two
parallel layers of opposing charges. This RC pair governs the mid-frequency
behavior (typically between 1 Hz and 1 kHz) and manifests as the first
semi-circular arc on the Nyquist plot.

*The second RC network*, or Warburg element, is associated with a much slower
phenomenon characterized by a significantly larger magnitude: the mass
transport driven by lithium ion diffusion. This dictates the low-frequency
response, appearing as a characteristic 45° phase shift in the Nyquist diagram,
which reflects the fact that diffusion becomes the rate-determining step in
this range.

Electrochemical processes linked to EIS regions (`Analog Devices <https://www.analog.com/media/en/reference-design-documentation/reference-designs/cn0510.pdf>`_, 2021)

Impact on Impedance of physical parameters from PXD model
=========================================================

The PXD BattMo model was developed using parameters derived from the study by
Chen et al. (2020). These parameters were subsequently modified to better
understand their influence on the impedance diagram.

Experimental impedance data were acquired using Electrochemical Impedance
Spectroscopy (EIS). EIS is a highly sensitive, non-destructive analytical
technique utilized to investigate the internal dynamics of electrochemical
systems, such as lithium-ion batteries. The method involves applying a
small-amplitude alternating current (AC) excitation across a broad frequency
spectrum, ranging from thousands of hertz down to fractions of a hertz, and
recording the subsequent response and phase shift.

Because the applied variations remain small, EIS assumes the system remains in
a pseudo-steady state, thereby enabling the use of Fourier transform analysis.

By analyzing the battery's response across distinct frequency bands, its
internal mechanisms can be diagnosed. Specifically, high frequencies reveal the
pure ohmic resistance of the components, mid-frequencies capture the
charge-transfer resistance, and low frequencies reflect mass transport and
diffusion processes, providing insight into the migration rate of lithium ions
through the solid electrode particles.

.. grid:: 2

   .. grid-item-card::
      :padding: 2

      :ref:`Battery Impedance Computation <runImpedanceScript>`
           
   .. grid-item-card::
      :padding: 2

      :ref:`ECM calibration <runCalibration>`
