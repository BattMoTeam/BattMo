==================
Battery Geometries
==================

We support multidimensional geometries (from 1D to 3D). The geometry is part of the input data and must be setup up
before a simulation. The base class :battmo:`BatteryGenerator` provides a template to construct the battery geometry,
which includes a mesh for all the components (Coatings, sepators, current collectors...) and the coupling between those
as the main output.

To create our grids, we rely essentially on the grid generation functions provided by `MRST
<https://www.sintef.no/Projectweb/MRST/>`_. We use also the visualization tools available there, see :ref:`here <visualization>`.

Standard geometries can often be parameterized, meaning that a small set parameters is enough to construct the whole
geometrical model. We have implemented some standard geometries and provide here examples with an illustration and a
list of the parameters used to produce this illustration. The parameters can be passed through json interface and aref
described in some more detail in the :battmofile:`geometry json schema<Utilities/JsonSchemas/Geometry.schema.json>`.


:battmo:`BatteryGeneratorP2D`
-----------------------------

The geometrical model is 1D. Here, for the sake of the illustration, we plot a corresponding 3D model.

.. image:: img/1dmodel.png
           :width: 80%
                   

.. list-table:: Parameters with values used in plot above
   :header-rows: 1

   * - Parameter
     - Name
     - Value
   * - length of negative current collector
     - :code:`x(1)`
     - 25 μm
   * - length of negative active material
     - :code:`x(2)`
     - 64 μm
   * - length of separator
     - :code:`x(3)`
     - 15 μm
   * - length of positive active material
     - :code:`x(4)`
     - 57 μm
   * - length of positive current collector
     - :code:`x(5)`
     - 15 μm
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
     - 3 cm^2

.. _2dgeometry:
         
:battmo:`BatteryGeneratorP3D`
-----------------------------

The geometrical model is 2D.

.. image:: img/runbattery2d.png
           :width: 80%
                   

.. list-table:: Parameters with values used in plot above
   :header-rows: 1

   * - Parameter
     - Name
     - Default value
   * - length of negative current collector
     - :code:`xlength(1)`
     - 10 μm
   * - length of negative active material
     - :code:`xlength(2)`
     - 100 μm
   * - length of separator
     - :code:`xlength(3)`
     - 50 μm
   * - length of positive active material
     - :code:`xlength(4)`
     - 80 μm
   * - length of positive current collector
     - :code:`xlength(5)`
     - 10 μm
   * - length in y direction
     - :code:`ylength`
     - 1 cm
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

A 3D geometrical model consising only of one layer with one tab on each current collector.

.. image:: img/runbattery3d.png
           :width: 80%

An illustration with different scalings for each axes which shows the different component:

.. image:: img/3dmodel.png
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

A geometry model for jelly rolls. Here, we have used parameters corresponding to th 4680 model.

.. image:: img/jellyrollmodel.png
           :width: 80%

.. list-table:: Thickness and discretization number (**N**) are passed through the dictionaries :code:`widthDict` and :code:`nrDict`, with values used in plot above.
   :header-rows: 1

   * - Component
     - Key name
     - Length 
     - N
   * - Separator
     - :code:`Separator`
     - 50 μm
     - 3
   * - Negative Electrode Coating
     - :code:`NegativeCoating`
     - 94 μm
     - 3
   * - Negative Electrode Current Collector
     - :code:`NegativeCurrentCollector`
     - 25 μm
     - 3
   * - Positive Electrode Coating
     - :code:`PositiveCoating`
     - 84 μm
     - 3
   * - Positive Electrode Current Collector
     - :code:`PositiveCurrentCollector`
     - 10 μm
     - 3
          
.. list-table:: Other parameters, with values used in plot above.
   :header-rows: 1
                 
   * - Parameter
     - Name
     - Value
   * - number of windings in the spiral
     - :code:`nwindings`
     - 52
   * - Inner Radius correspoding to the empty space in the middle
     - :code:`rInner`
     - 2 mm
   * - Height of the battery
     - :code:`L`
     - 70 mm
   * - number of cells in the angular direction
     - :code:`nas`
     - 20
   * - number of discretization cells in the longitudonal
     - :code:`nL`
     - 10

There are parameters for the tabs that we do not detail here, see :battmofile:`json
schema<Utilities/JsonSchemas/Geometry.schema.json#242>`. The json source is available
:battmofile:`here<Documentation/scripts/jsonfiles/4680-geometry.json>`.

.. _coincell:
      
:battmo:`CoinCellBatteryGenerator`
----------------------------------

A geometrical model for a coin cell.

.. image:: img/coincell.png
           :width: 80%


.. list-table:: Parameters for each component : thickness, diameter, number of cell layers (Nl), with the values used in the plot above (a CR 2016 coin cell)
   :header-rows: 1
                 
   * - Component
     - Name
     - Thickness
     - Diameter
     - Nl
   * - Negative electrode current collector
     - :code:`NegativeCurrentCollector`
     - 0.73 mm
     - 20 mm
     - 9
   * - Negative electrode coating
     - :code:`NegativeCoating`
     - 50 μm
     - 16 mm
     - 3
   * - Separator
     - :code:`Separator`
     - 20 μm
     - 18 mm
     - 2
   * - Positive electrode coating
     - :code:`PositiveCoating`
     - 67 μm
     - 16 mm
     - 3
   * - Positive electrode current collector
     - :code:`PositiveCurrentCollector`
     - 0.73 mm
     - 20 mm
     - 9
   
.. list-table:: Other parameters, with values used in plot above
   :header-rows: 1
                 
   * - Parameter
     - Name
     - Value
   * - Discretization number in radial direction
     - :code:`numRadial`
     - 7
   * - Discretization number in angular direction
     - :code:`numAngular`
     - 30
