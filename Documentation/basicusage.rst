===========
Basic Usage
===========

In this section, we describe how to run |battmo| using json inputs. The specification of json inputs is described in :ref:`json:Json input specification`

First Run
=========

.. container:: toggle

    .. container:: header

       **Show/Hide Code**

    .. code-block:: xml
       :linenos:

       from plone import api

First, we load and parse the json input file using the commands

.. code:: matlab

   jsonstruct = parseBattmoJson('input.json')

Then, we run it as follows
  
.. code:: matlab

   output = runBatteryJson(filename)

.. collapse:: The json file :code:`input.json` provides all the data need to run the whole simulation.

   .. literalinclude:: inputfile.json
      :language: json
                 
:todo:`Describe the output structure (what it contains)`
      
To visualize the results, we refer to the :todo:`add link`
     

