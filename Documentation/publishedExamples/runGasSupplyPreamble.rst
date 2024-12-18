Model overview
==============

We present here a simulation model for a single gas layer. The complete electroyser has two such layers, see
illustration below.

.. figure:: /img/gas-layer-illustration.png  
   :target: ../_images/gas-layer-illustration.png
   :align: center
           
   Illustration of the electrolyser with the two gas layer domains

We consider the gas layer alone, without coupling with the electrodes. We will choose here the gas layer used at the
anode. We inject a mixture of two gases (H2O and O2). We model the gas supply layer as a porous media. The flow is then driven by a pressure gradient and the effect of diffusion.
   
Governing equations
-------------------

The governing equations are given by mass conservation. We have

.. math::
   \frac{\partial\rho x_i}{\partial t} + \nabla\cdot\left(j_{\text{conv},i} + j_{\text{diff}, i}\right) = 0.

where :math:`\rho` denotes the density, :math:`x_i` the mass fraction of the component :math:`i`, for
:math:`i=\{1,2\}`. The convection and diffusion fluxes are denoted with :math:`j_{\text{conv},i}` and
:math:`j_{\text{diff}, i}`, respectively.

We use a simple equation of state for the gases, as we assume that they behave as ideal gas.

For the momentum equation, we use the Darcy approximation which covers the case of a flow in a porous media. The Darcy velocity is given by

.. math::
   u = -\frac{K}{\mu}\nabla p,

where :math:`K` denotes the permeability and :math:`mu` the viscosity. We are working on a model which uses the Stokes equations, which can be used for a free flow.

The convection term for each species is given by

.. math::
   u_{\text{conv}, i} =  \rho x_i u.

For the diffusion term, we use a Fick diffusion model and the diffusion flux for each component is given by

.. math::
   \hat j_{\text{diff}, i} = -D_i \nabla c_i

where :math:`c_i` is the concentration of the specie :math:`i` and :math:`\hat j_{\text{diff}, i}` is the flux in
:math:`mol/m^2/s`. We convert this flux in a mass flux, using the assumption that the gas are ideal gas, to obtain

.. math::
   j_{\text{diff}, i} = -D_i \nabla (\rho x_i)

Our governing equations can now be rewritten using the explicit expression for the fluxes,

.. math::
   \frac{\partial\rho x_i}{\partial t} + \nabla\cdot\left(\rho x_i u - D_i\nabla (\rho x_i)\right) = 0

This system of equation is discretized using finite volume and solved in time using backward Euler, for stability. The
unknowns are the pressure and one of the mass fractions.
   
Boundary conditions
-------------------

We have implemented boundary equations of three types, which can be used for both inlet and outlet:

* total mass rate and composition (mass fractions)
* pressure and composition (mass fractions)
* pressure and zero mass fraction derivative condition

The boundary conditions are applied on a subset of external faces in the grid (see drawing below). The two first types
of boundary conditions are standard.

When we impose the rate and composition, the given rate corresponds to the sum of all the input mass rates of every face
which belongs to the boundary condition (It means that each boundary face individually can see different mass
rates). This condition is use for the inlet in the example below, see json input :battmofile:`here <ProtonicMembrane/jsonfiles/gas-supply.json#16>`.

The pressure and composition boundary condition corresponds to a standard Dirichlet condition. We would typically use
such condition at the outlet. Even for a convection dominated flow, the composition at the boundary is given. This could
be seen as an artifact because, in an experimental setup, we do not control the composition directly at the outlet, and
we do not want to model and simulate the flow outside the gas layer region. Therefore, we introduced the third boundary
condition.

When we use a *zero mass fraction derivative* boundary condition, we impose that the *derivative* of the mass fraction
in the normal direction at the boundary is zero, :math:`\frac{\partial x_i}{\partial n} = 0`. In this way, we model a
situation where the outlet is connected to a long pipe which is itself connected to the ambient room, where we have a
given composition. The pipe is long enough to contain the diffusion effects which


