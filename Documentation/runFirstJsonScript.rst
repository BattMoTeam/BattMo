
.. _runFirstJsonScript:
   
=============
First Example
=============

:todo:`Add page for schema and link it here`. Some standard data sets are included :todo:`add page on that and link it here`
  
First, we load and parse the json input file using the commands

.. code:: matlab

   jsonstruct = parseBattmoJson('input.json')

Then, we run it as follows
  
.. code:: matlab

   output = runBatteryJson(filename)

The json file :code:`input.json` in this case describes the whole simulation. To get an idea of how such file looks like, get a look :ref:`here<jsoninputfile>`.

:todo:`Describe the output structure (what it contains)`
      
To visualize the results, we refer to the :todo:`add link`
     
