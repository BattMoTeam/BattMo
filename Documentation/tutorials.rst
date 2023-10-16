===============
Getting Started
===============
           
Here, we present some guidelines and tutorials to get started with and give a short description of the example scripts
that are available in |battmo|

* |battmo| supports **json inputs**. The json input is specified in schemas, which also can be used to validate your input
  data. :todo:`Add page for schema and link it here`. Some standard data sets are included :todo:`add page on that and link it here`
  and add much more

  First, we load and parse the json input file using the commands

  .. code:: matlab

     filename = fullfile('ParameterData', 'BatteryCellParameters', 'LithiumIonBatteryCell', 'lithium_ion_battery_nmc_graphite.json')
     jsonstruct = parseBattmoJson(filename)

  Then, we run it as follows
  
  .. code:: matlab

     runBatteryJson(filename)

  To visualize the results, we refer to the :todo:`add link`
     
* In this :ref:`json example <runJsonScript>`, we show how to combine different json files
  
* The :ref:`battMoTutorial <battMoTutorial>` describes the different step of a simulation when it is set manually beyond the options that
  are given by the json input format. They consist of

  * Setup model parameters (physical parameters and modeling choice)
  * Setup the geometry (BattMo supports a variety of geometry)
  * Initialise the model
  * Setup a simulation schedule (time stepping and control)
  * Setup the initial state
  * Run the simulation
  * Visualize the results
    

  .. image:: img/battMoTutorial.png
             :width: 80%

* In this :ref:`example <runBattery1D>`, more options are used such as some settings of the non-linear solver.

* Composite material example

  :todo:`add here or maybe in gallery`

* More Examples

  In the :code:`Examples` directory examples, the following examples can be found. They should run out of the box (write
  the name of script after :ref:`installing BattMo<installation>`). If not, tell us!
  
  - :code:`runBattery2D` : 2D example using :ref:`2D model geometry<2dgeometry>`
  - :code:`runBattery3D` : 3D example using :ref:`3D model geometry<3dgeometry>`
  - :code:`runChen2020` : Example using Chen data including a comparison with `pybamm <https://www.pybamm.org/>`_
  - :code:`runCR` : Example running for a coin cell
  - :code:`runGittTest` : Example running a Gitt test
  - :code:`runJellyRoll` : Example running a :ref:`jelly roll geometry <jellyroll>`


.. toctree::
   :hidden:

   publishedExamples/runJsonScript   
   publishedExamples/battMoTutorial
   publishedExamples/runBattery1D
   publishedExamples/runSEIActiveMaterial
   
