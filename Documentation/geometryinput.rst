==================
Battery Geometries
==================

We support multidimensional geometries (from 1D to 3D). The geometry is part of the input data and must be setup up
before a simulation. The base class :battmo:`BatteryGenerator` provides a template to construct the battery geometry,
which includes a mesh for all the components (Coatings, sepators, current collectors...) and the coupling between those
as the main output. Standard geometries can often be parameterized, meaning that a small set parameters is enough to
construct the whole geometrical model. We have implemented some standard geometries and provide here examples with an
illustration and a list of the parameters used to produce this illustration. The parameters can be passed through json
interface and are described in some more detail in the :battmofile:`geometry json
schema<Utilities/JsonSchemas/Geometry.schema.json>`.


:battmo:`BatteryGeneratorP2D`
-----------------------------

.. image:: img/1dmodel.png
           :width: 80%
                   

.. list-table:: Parameters with values used in plot above
   :header-rows: 1

   * - Parameter
     - Name
     - Value
   * - length of negative current collector
     - :code:`x(1)`
     - 25 micro meter
   * - length of negative active material
     - :code:`x(2)`
     - 64 micro meter
   * - length of separator
     - :code:`x(3)`
     - 15 micro meter
   * - length of positive active material
     - :code:`x(4)`
     - 57 micro meter
   * - length of positive current collector
     - :code:`x(5)`
     - 15 micro meter
   * - discretization number for negative current collector
     - :code:`ccnen`
     - 10
   * - discretization number for negative active material
     - :code:`nenx`
     - 10
   * - discretization number for separator
     - :code:`sepnx`
     - 10
   * - discretization number for positive current collector
     - :code:`penx`
     - 10
   * - discretization number for positive active material
     - :code:`ccpenx`
     - 10
   * - Face area
     - :code:`faceArea`
     - 2*1.6387e-04;

.. _2dgeometry:
         
:battmo:`BatteryGeneratorP3D`
-----------------------------

.. image:: img/runbattery2d.png
           :width: 80%
                   

.. list-table:: Parameters with values used in plot above
   :header-rows: 1

   * - Parameter
     - Name
     - Default value
   * - length of negative current collector
     - :code:`xlength(1)`
     - 10 micro meter
   * - length of negative active material
     - :code:`xlength(2)`
     - 100 micro meter
   * - length of separator
     - :code:`xlength(3)`
     - 50 micro meter
   * - length of positive active material
     - :code:`xlength(4)`
     - 80 micro meter
   * - length of positive current collector
     - :code:`xlength(5)`
     - 10 micro meter
   * - length in y direction
     - :code:`ylength`
     - 1 centi meter
   * - discretization number for negative current collector
     - :code:`ccnenx`
     - 10
   * - discretization number for negative active material
     - :code:`nenx`
     - 10
   * - discretization number for separator
     - :code:`sepnx`
     - 10
   * - discretization number for positive current collector
     - :code:`penx`
     - 10
   * - discretization number for positive active material
     - :code:`ccpenx`
     - 10
   * - discretization number in y direction
     - :code:`ny`
     - 10   

                   
.. _3dgeometry:
      
:battmo:`BatteryGeneratorP4D`
-----------------------------

.. image:: img/runbattery3d.png
           :width: 80%
                   
.. list-table:: Parameters with values used in plot above
   :header-rows: 1

   * - Parameter
     - Name
     - Default value
   * - x-length of first tab
     - :code:`x(1)`
     - 4 cm
   * - x-length between the tabs
     - :code:`x(2)`
     - 2 cm
   * - x-length of last tab
     - :code:`x(3)`
     - 4cm
   * - y-length of first tab
     - :code:`y(1)`
     - 1mn
   * - y-length between the tabs
     - :code:`y(2)`
     - 2 cm
   * - y-length of last tab
     - :code:`y(3)`
     - 1 mm        
   * - length of negative current collector
     - :code:`z(1)`
     - 10 μm
   * - length of negative active material
     - :code:`z(2)`
     - 100 μm
   * - length of separator
     - :code:`z(3)`
     - 50 μm
   * - length of positive active material
     - :code:`z(4)`
     - 80 μm
   * - length of positive current collector
     - :code:`z(5)`
     - 10 μm
   * - discretization number in z-direction for separator
     - :code:`sep_nz`
     - 3 
   * - discretization number in z-direction for positive active material
     - :code:`ne_am_nz`
     - 3 
   * - discretization number in z-direction for negative active material
     - :code:`pe_am_nz`
     - 3 
   * - discretization number in z-direction for negative current collector
     - :code:`ne_cc_nz`
     - 2 
   * - discretization number in z-direction for positive current collector
     - :code:`pe_cc_nz`
     - 2 
   * - discretization number in x-direction interior region
     - :code:`int_elyte_nx`
     - 3 
   * - discretization number in x-direction negative tab region
     - :code:`ne_cc_nx`
     - 3 
   * - discretization number in x-direction positive tab region
     - :code:`pe_cc_nx`
     - 3 
   * - discretization number in y-direction interior region
     - :code:`elyte_ny`
     - 4 
   * - discretization number in y-direction negative tab region
     - :code:`ne_cc_ny`
     - 2 
   * - discretization number in y-direction positive tab region
     - :code:`pe_cc_ny`
     - 2 
                   
.. _jellyroll:
      
:battmo:`SpiralBatteryGenerator`
--------------------------------

.. image:: img/jellyrollmodel.png
           :width: 80%

.. list-table:: Thickness and discretization number (Nr) are passed through the dictionaries :code:`widthDict` and :code:`nrDict`, with values used in plot above.
   :header-rows: 1

   * - Component
     - Key name
     - Length 
     - Nr
   * - Separator
     - :code:`Separator`
     - 1 mm
     - 1 mm
   * - Negative Electrode Coating
     - :code:`NegativeCoating`
     - 1 mm
     - 1 mm
   * - Negative Electrode Current Collector
     - :code:`NegativeCurrentCollector`
     - 1 mm
     - 1 mm
   * - Positive Electrode Coating
     - :code:`PositiveCoating`
     - 1 mm
     - 1 mm
   * - Positive Electrode Current Collector
     - :code:`PositiveCurrentCollector`
     - 1mm
     - 1mm
          
.. list-table:: Other parameters, with values used in plot above.
   :header-rows: 1
                 
   * - Parameter
     - Name
     - Value
   * - number of windings in the spiral
     - :code:`nwindings`
     - 1
   * - Inner Radius correspoding to the empty space in the middle
     - :code:`rInner`
     - 1
   * - length of the battery
     - :code:`L`
     - 1
   * - number of cells in the angular direction
     - :code:`nas`
     - 1
   * - number of discretization cells in the longitudonal
     - :code:`nL`
     - 1

There are parameters for the tabs that we do not detail here, see :battmofile:`json schema<Utilities/JsonSchemas/Geometry.schema.json#242>`

.. _coincell:
      
:battmo:`CoinCellBatteryGenerator`
----------------------------------

.. image:: img/coincell.png
           :width: 80%


.. list-table:: Parameters for each component : thicknes, diameter, number of cell layers (Nl), with the values used in the plot above
   :header-rows: 1
                 
   * - Component
     - Name
     - Thickness
     - Diameter
     - Nl
   * - Negative electrode current collector
     - :code:`NegativeCurrentCollector`
     - 1
     - 1
     - 1
   * - Negative electrode coating
     - :code:`NegativeCoating`
     - 1
     - 1
     - 1
   * - Separator
     - :code:`Separator`
     - 1
     - 1
     - 1
   * - Positive electrode coating
     - :code:`PositiveCoating`
     - 1
     - 1
     - 1
   * - Positive electrode current collector
     - :code:`PositiveCurrentCollector`
     - 1
     - 1
     - 1
   
.. list-table:: Other parameters, with values used in plot above
   :header-rows: 1
                 
   * - Parameter
     - Name
     - Value
   * - Discretization number in radial direction
     - :code:`numRadial`
     - 1
   * - Discretization number in angular direction
     - :code:`numAngular`
     - 1       


         
