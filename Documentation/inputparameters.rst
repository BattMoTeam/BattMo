============
BattMo Input
============

.. toctree::
   :maxdepth: 2
   :hidden:

   Introduction <self>
   jsonexample
   json
   mergejsonstruct
   jsonfiles
   functioninterface
   controlinput
   parsets


BattMo uses a simple object dictionary structure as unique input for a whole simulation. The structure is easily
editable programmatically in matlab. As a file, it is written using the `json format
<https://en.wikipedia.org/wiki/JSON>`_ provides a basic but yet extremly flexible way to provide a dictionary-like data.

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
probably to look at examples

In :ref:`My first Json input file<jsonexample:My first Json Input>`, we guide you through some selected examples. In
:ref:`List of Json File Examples<jsonfiles:List of Json File Examples>`, we provide links to json datasets available in
your BattMo installationa and corresponding to different components in the simulation (material properties, geometry,
controls ...).

Combining different inputs is a standard way to create you own input, for the model you are targetting, and is is easily
achieved using the function `mergeJsonStructs`, see :ref:`Merging Parameters<mergejsonstruct:Merging parameters>`.

.. grid:: 2

   .. grid-item-card::
      :padding: 2
      
      :ref:`My first Json input file<jsonexample:My first Json Input>`

   .. grid-item-card::
      :padding: 2
      
      :ref:`JSON input specification<json:JSON input specification>`
          
   .. grid-item-card::
      :padding: 2
      
      :ref:`Merging Parameters<mergejsonstruct:Merging parameters>`

   .. grid-item-card::
      :padding: 2
      
      :ref:`List of Json File Examples<jsonfiles:List of Json File Examples>`
           
   .. grid-item-card::
      :padding: 2
      
      :ref:`Functional Parameters<functioninterface:Functional parameters>`

   .. grid-item-card::
      :padding: 2
      
      :ref:`Control models<controlinput:Control models>`

   .. grid-item-card::
      :padding: 2

      :ref:`Parameter Sets <parsets:Parameter sets>`
           
           


