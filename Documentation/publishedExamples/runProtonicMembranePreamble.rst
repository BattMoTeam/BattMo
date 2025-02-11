Model overview
==============

We consider the model of a mixed proton and electron conducting membrane, as described in :cite:`V_llestad_2019` from 2019.
        
Governing equations
-------------------

In the membrane, we have three components or *species* given by the proton (:math:`H^+`), and the :math:`p` and :math:`n`
type charge carriers. We need expressions for the fluxes for each of those. The governing equations will then be given
by charge and mass conservation equations.

We denote by :math:`\phi` the electrostatic potential. For each of the components :math:`\alpha=\{H^+, p, n\}`, we
introduce the electrochemical potential denoted :math:`\bar\mu_\alpha` and the chemical potential denoted
:math:`\mu_\alpha`. The potentials are related with each other through the relation

.. math::

   \bar\mu_\alpha = \mu_\alpha + z_\alpha F \phi.

The fluxes are governed by the gradient of the electrochemical potential. We have

.. math::
   
   j_{\alpha} = -k_\alpha\nabla\bar\mu_\alpha,

for some coefficient :math:`k_\alpha` which is not necessarily a constant.

We assume that the chemical potential of :math:`H^+`, that is :math:`\mu_{H^+}`, is constant in the electrolyte and that
:math:`k_{H^+}` is constant. Hence, the current density of :math:`H^+` is given

.. math::          

    i_{H^+} = -\sigma_{H^+} \nabla \phi,

for a constant conductivity :math:`\sigma_{H^+}`. We denote by :math:`E = -\frac{\mu_n}{F}` the electronic chemical potential and :math:`\pi = E + \phi` the electromotive potential. Since we have equilibrium for the n-p reaction, we have the identity

.. math::

   d\mu_p + d\mu_n = 0.

For the :math:`p` and :math:`n` type conductivities,  we use the empirical relations

.. math::

   \sigma_p(E) = \sigma_p^0\exp\left(\frac{F(E - E_{\text{ref},p})}{RT}\right)\quad\text{ and }\quad\sigma_n(E) = \sigma_n^0\exp\left(\frac{-F(E - E_{\text{ref},n})}{RT}\right).

Here :math:`E_{\text{ref},p}` and :math:`E_{\text{ref},n}` are two reference potentials. 

We consider the steady state. We could introduce later charge and mass capacitors. The unknowns are the functions :math:`\phi(x)` and :math:`E(x)` in the electrolyte. The governing equations are given by the mass conservation for the proton and the charge conservation.

The **mass conservation equation** for :math:`H^+` is given by

.. math::

   \nabla\cdot j_{H^+} = 0.

The total current density is given by :math:`i = F (j_{H^+} + j_p - j_n)` and the **charge conservation equation** is

.. math::

   \nabla\cdot i = 0.
   
We introduce the electronic current density :math:`i_{\text{el}} = F(j_p - j_n)` and we get

.. math::
   
   \nabla\cdot i_{\text{el}} = 0.

We can rewrite :math:`i_{\text{el}}` as

.. math::
   
   i_{\text{el}} = - (\sigma_p(E) + \sigma_n(E))\nabla ( E + \phi ).

We finally obtain the governing equations for :math:`\phi(x)` and :math:`\pi(x)` as the following system of differential
equations

.. math::

   \begin{align}
   -\nabla\cdot(\sigma_{H^+}\nabla\phi) &= 0,\\
   -\nabla\cdot((\sigma_p(E) + \sigma_n(E))\nabla \pi) &= 0.
   \end{align}   

Boundary conditions
-------------------

We define the over-potential :math:`\eta` as

.. math::
    
   \eta_\text{elde} = \pi_\text{elde} - \phi_\text{elde} - \text{OCP}_\text{elde}

where :math:`\text{OCP}_\text{elde}` is the open-circuit potential for the given electrode. The value of the
OCP at each electrode depend on the composition at the electrode (see :cite:`V_llestad_2019` for the expressions).

At the anode, we imposte that the proton current  is given through the following Buttler-Volmer type expression

.. math::
   i_{H^+, \text{an}} = -i_0\frac{e^{-\beta\frac{ z F \eta_{\text{an}}}{RT}} - e^{( 1- \beta)\frac{ z F \eta_{\text{an}}}{RT}}}{ 1+ \frac{i_0}{i_{l,c}}e^{-\beta\frac{ z F \eta_{\text{an}}}{RT}} -  \frac{i_0}{i_{l,a}}e^{( 1- \beta)\frac{ z F \eta_{\text{an}}}{RT}}}


Here, :math:`i_{l,c}` and :math:`i_{l,a}` are given constants. The value of the reference current density :math:`i_0` is
also constant.

The total current is given by :math:`i_{\text{an}} = I` for some constant current :math:`I`.

At the cathode, we impose that the electrostatic potential is equal to zero and 
a relation between the :math:`H^+` current and the over-potential that takes a linear form,

.. math::
   
   R_\text{CT} \ i_{H^+, \text{ct}} = \eta_\text{ct},

for a given charge transfer constant :math:`R_\text{CT}`.
   
