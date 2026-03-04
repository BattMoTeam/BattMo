==============
Control models
==============

.. toctree::
   :maxdepth: 2
   :hidden:
      
   runControlExamples.nblink
   runTimeControlExample.nblink
   runGenericStepControlSimple.nblink
   runGenericStepControlCycle.nblink
      
The battery :battmo:`model <Battery#19>` contains a *Control* sub-model, see the overall model structure
:ref:`description<architecture:BattMo Model Architecture>`. The Control model determine the control *values* and *type*.

At a **given time**, there are in general two control types:

* Total current
* Voltage  

Given a control type, the control *value* can change in time. For the most standard control types, the json interface
can be used. We plan to include there more control models. Below, we give some short explanations on how a control model
can be implemented.

The most standard controls can be called from BattMo using the json interface:

* Constant Current Discharge (CCDischarge)
* Constant Current Charge (CCCharge)
* Constant Current Constant Voltage (CCCV)

We also have a generic control setup, see example below.

The parameters for each of the model are described in the json schema :battmofile:`ControlModel.schema.json
<Utilities/JsonSchemas/ControlModel.schema.json>`. 
   
.. grid:: 2

   .. grid-item-card::
      :padding: 2

      :ref:`Control Examples <runControlExamples>`
           
   .. grid-item-card::
      :padding: 2

      :ref:`Time Functional Control Examples <runTimeControlExample>`

   .. grid-item-card::
      :padding: 2

      :ref:`Generic Control Simple Example <runGenericStepControlSimple>`
              
   .. grid-item-card::
      :padding: 2

      :ref:`Generic Control Cycle Example <runGenericStepControlCycle>`
              
