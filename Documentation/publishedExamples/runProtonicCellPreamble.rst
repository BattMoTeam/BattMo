Model coupling
==============

We couple the :ref:`gas supply layer <runGasSupply>` with the anode of the :ref:`proton membrane
<runProtonicMembrane>`. The coupling occurs at the interface of the two regions. Along the anode, we assume a constant
voltage, that is, we neglect the potential loss in this part.

The chemical reaction at the interface is 

.. math::

   \ce{1/2 H2O <-> H^+ + e^- + 1/4 O2}

It relates the mass fluxes in the gas supply layer and inside the anode in the following way

.. math::

   \begin{array}{rcl}
    j_\ce{H2O}\cdot n &=& 2j_{\ce{H^+}}\cdot n\\
    j_{\ce{O2}}\cdot n &=& -4j_{\ce{H^+}}\cdot n
   \end{array}

Here, :math:`n` denotes the normal at the interface (same orientation in the whole expression), so that, for example,
:math:`j_\ce{H2O}\cdot n` denotes the mass of :math:`\ce{H2O}` leaving the gas layer in :math:`kg/m^2/s` (for a 3D
model). For the gas layer, the fluxes are given by the sum of the convection and diffusion fluxes, so that, for
:math:`\ce{H2O}`, we have

.. math::

    j_\ce{H2O}\cdot n = -D \nabla(\rho x_{\ce{H2O}})\cdot n - \frac{K}{\mu}\rho x_i \nabla p\cdot n

and a similar expression for the flux of :math:`\ce{O2}` at the interface.



