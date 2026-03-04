============
BattMo Input
============

.. toctree::
   :maxdepth: 2
   :hidden:

   json
   mergejsonstruct
   functioninterface
   controlinput
   parsets

BattMo uses a simple dictionary structure as input for a simulation. In matlab, it means that we use a simple `struct`.

The `json format <https://en.wikipedia.org/wiki/JSON>`_ provides a basic but yet extremly flexible way to provide a
dictionary like data

As seen in :ref:`Your First Example<You First Example>`, you can simply provide the path of your input data using the
command :code:`parseBattmoJson`

.. code:: matlab

   jsonstruct = parseBattmoJson('Examples/JsonDataFiles/sample_input.json')

Then, you can modify this structure as you want within Matlab before sending it to the simulator

.. code:: matlab

   output = runBatteryJson(jsonstruct)

We use json schema to describe the expected keys or field names in the input, see :ref:`JSON input
specification<json:JSON input specification>`. All the schemas are collected under the directory
:battmofile:`JsonSchemas<Utilities/JsonSchemas>`. However, the easiest way to discover the available input parameters is
probably to look at examples and we invite you to look at the list below. Most of the parameters in the list are
components of a complete inputs, which can be aggregated using the function :battmo:`mergeJsonStructs` to setup a full
input, see :ref:`Merging Parameters<mergejsonstruct:Merging parameters>`.

.. list-table::
   :header-rows: 1
   :width: 90%
   :align:  center
            
   * - Full Cell Parameters
   * - :battmofile:`Graphite-NMC cell <ParameterData/BatteryCellParameters/LithiumIonBatteryCell/lithium_ion_battery_nmc_graphite.json>`   
   * - :battmofile:`Xu 2015 cell <ParameterData/ParameterSets/Xu2015/lfp.json>`
   * - :battmofile:`Chen 2020 cell <ParameterData/ParameterSets/Chen2020/chen2020_lithium_ion_battery.json>`
   * - :battmofile:`Graphite-Silicon cell <ParameterData/BatteryCellParameters/LithiumIonBatteryCell/lithium_ion_battery_nmc_silicon.json>`
     
.. list-table::
   :header-rows: 1
   :width: 90%
   :align:  center
            
   * - Geometry Parameters
   * - :battmofile:`1D model <Examples/JsonDataFiles/geometry1d.json>`
   * - :battmofile:`3D model <Examples/JsonDataFiles/geometry3d.json>`
   * - :battmofile:`Multi-layer pouch cell <Examples/JsonDataFiles/geometryMultiLayerPouch.json>`

.. list-table::
   :header-rows: 1
   :width: 90%
   :align:  center
            
   * - Material Parameters
   * - :battmofile:`Graphite <ParameterData/MaterialProperties/Graphite/graphite.json>`
   * - :battmofile:`LFP (Xu 2015) <ParameterData/MaterialProperties/LFP/LFP_Xu2015.json>`
   * - :battmofile:`Celgard2500 Electrolyte <ParameterData/BatteryComponentParameters/celgard2500.json>`
     
.. grid:: 2

   .. grid-item-card::
      :padding: 2
      
      :ref:`JSON input specification<json:JSON input specification>`
          
   .. grid-item-card::
      :padding: 2
      
      :ref:`Merging Parameters<mergejsonstruct:Merging parameters>`

   .. grid-item-card::
      :padding: 2
      
      :ref:`Functional Parameters<functioninterface:Functional parameters>`

   .. grid-item-card::
      :padding: 2
      
      :ref:`Control models<controlinput:Control models>`

   .. grid-item-card::
      :padding: 2

      :ref:`Parameter Sets <parsets:Parameter sets>`
           
           
