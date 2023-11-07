===========
Basic Usage
===========

In this section, we describe how to run |battmo| using json inputs. The specification of json inputs is described in :ref:`json:Json input specification`

First Run
=========

First, we load and parse the json input file using the commands

.. code:: matlab

   jsonstruct = parseBattmoJson('input.json')

Then, we run it as follows
  
.. code:: matlab

   output = runBatteryJson(jsonstruct)
                      
The json file :code:`input.json` given :ref:`here <jsoninputfile>` provides the needed data for the simulation. The
section :ref:`json:Json input specification` gives a complete description of each of the field. We have used long and
explicit names for a good readability.

The :code:`output` structure returns among other thing the model and the states. 

.. code:: matlab

   model  : [1x1 Battery]
   states : [1x1 struct]
          
We can plot the results using :battmo:`plotDashboard`

.. code:: matlab

   plotDashboard(model, states)

This is a snapshot of the result


   





                 
:todo:`Describe the output structure (what it contains)`
      
To visualize the results, we refer to the :todo:`add link`
     

