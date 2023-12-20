==================
Battery Geometries
==================

Description of the models that are used to setup the geometry. Here, we include a short list of the geometrical
parameters. We refer to the source code for more details. We plan to include more geometries.


:battmo:`BatteryGeneratorP2D`
-----------------------------

.. image:: img/1dmodel.png
           :width: 80%
                   

.. list-table:: Parameters
   :header-rows: 1

   * - Parameter
     - Name
     - Default value
   * - length of negative current collector
     - x(1)
     - 25 micro meter
   * - length of negative active material
     - x(2)
     - 64 micro meter
   * - length of separator
     - x(3)
     - 15 micro meter
   * - length of positive active material
     - x(4)
     - 57 micro meter
   * - length of positive current collector
     - x(5)
     - 15 micro meter
   * - discretization number for negative current collector
     - ccnenx   
     - 10
   * - discretization number for negative active material
     - nenx    
     - 10
   * - discretization number for separator
     - sepnx    
     - 10
   * - discretization number for positive current collector
     - penx  
     - 10
   * - discretization number for positive active material
     - ccpenx
     - 10
   * - Face area
     - faceArea
     - 2*1.6387e-04;

.. _2dgeometry:
         
:battmo:`BatteryGeneratorP3D`
-----------------------------

.. image:: img/runbattery2d.png
           :width: 80%
                   

.. list-table:: Parameters
   :header-rows: 1

   * - Parameter
     - Name
     - Default value
   * - length of negative current collector
     - xlength(1)
     - 10 micro meter
   * - length of negative active material
     - xlength(2)
     - 100 micro meter
   * - length of separator
     - xlength(3)
     - 50 micro meter
   * - length of positive active material
     - xlength(4)
     - 80 micro meter
   * - length of positive current collector
     - xlength(5)
     - 10 micro meter
   * - length in y direction
     - ylength
     - 1 centi meter
   * - discretization number for negative current collector
     - ccnenx   
     - 10
   * - discretization number for negative active material
     - nenx    
     - 10
   * - discretization number for separator
     - sepnx    
     - 10
   * - discretization number for positive current collector
     - penx  
     - 10
   * - discretization number for positive active material
     - ccpenx
     - 10
   * - discretization number in y direction
     - ny
     - 10   

                   
.. _3dgeometry:
      
:battmo:`BatteryGeneratorP4D`
-----------------------------

.. image:: img/runbattery3d.png
           :width: 80%
                   
.. list-table:: Parameters
   :header-rows: 1

   * - Parameter
     - Name
     - Default value
   * - x-length of first tab
     - x(1)
     - 4 cm
   * - x-length between the tabs
     - x(2)
     - 2 cm
   * - x-length of last tab
     - x(3)
     - 4cm
   * - y-length of first tab
     - y(1)
     - 1mn
   * - y-length between the tabs
     - y(2)
     - 2 cm
   * - y-length of last tab
     - y(3)
     - 1 mm        
   * - length of negative current collector
     - z(1)
     - 10 μm
   * - length of negative active material
     - z(2)
     - 100 μm
   * - length of separator
     - z(3)
     - 50 μm
   * - length of positive active material
     - z(4)
     - 80 μm
   * - length of positive current collector
     - z(5)
     - 10 μm
   * - discretization number in z-direction for separator
     - sep_nz
     - 3 
   * - discretization number in z-direction for positive active material
     - ne_am_nz
     - 3 
   * - discretization number in z-direction for negative active material
     - pe_am_nz
     - 3 
   * - discretization number in z-direction for negative current collector
     - ne_cc_nz
     - 2 
   * - discretization number in z-direction for positive current collector
     - pe_cc_nz
     - 2 
   * - discretization number in x-direction interior region
     - int_elyte_nx
     - 3 
   * - discretization number in x-direction negative tab region
     - ne_cc_nx
     - 3 
   * - discretization number in x-direction positive tab region
     - pe_cc_nx
     - 3 
   * - discretization number in y-direction interior region
     - elyte_ny
     - 4 
   * - discretization number in y-direction negative tab region
     - ne_cc_ny
     - 2 
   * - discretization number in y-direction positive tab region
     - pe_cc_ny
     - 2 
                   
.. _jellyroll:
      
:battmo:`SpiralBatteryGenerator`
--------------------------------

.. image:: img/jellyrollmodel.png
           :width: 80%
                   

.. _coincell:
      
:battmo:`CoinCellBatteryGenerator`
----------------------------------

.. image:: img/coincell.png
           :width: 80%

         
:battmo:`Base class<BatteryGenerator>`
--------------------------------------

This is the base class that gather the methods to setup the different grid. This class will be usefull if you want to
setup your own tailored grid.

         
