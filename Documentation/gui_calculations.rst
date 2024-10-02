================
Calculations
================
Multiple quick calculations are done within the app to calculate certain parameters.
In this section you can find the formulas that belong to these calculations. 

**Effective density**
   The effective density of an electrode coating:

   .. math::

      \rho = \sum_i \omega_{i} \rho_i

   Here is :math:`\omega_{i}` the mass fraction of an electrode coating material and :math:`\rho_{i}` the density of the same coating material.

**Mass loading**
   The mass loading of an electrode:

   .. math::

      \rho A = d \rho(1-\varepsilon)

   Here is :math:`d` the thickness of the electrode coating, :math:`\rho` is the effective density, and :math:`\varepsilon` is the porosity.

**Specific capacity**
   The specific capacity of an active material:

   .. math::

      Q_s = c_{max} |x_{max} - x_{min}| \frac{nF}{\rho}

   Here is :math:`c_{max}` the maximum concentration, :math:`x_{max}` and :math:`x_{min}` are the maximum and minimum stoichiometry,
   :math:`F` is the Faraday constant, :math:`\rho` the density of the active material, and :math:`n` is the number of electrons transfered.


**Capacity**
   The capacity of an electrode:

   .. math::

      Q =  \omega Q_s \rho V (1-\varepsilon)

   Here is :math:`\omega` the mass fraction, :math:`Q_s` is the specific capacity, :math:`\rho` the effective density, :math:`V` is the volume of the electrode, and :math:`\varepsilon` is the porosity.  

   The capacity of the cell is the limiting capacity out of the thwo electrode capacities.

**Cell mass**
   The cell mass:

   .. math::

      m = \sum_i m_{i}

   :math:`m_i` represents here the mass of the several battery components: negative and positive electrode, separator, electrolyte, current collector, and the packing mass.
   
   