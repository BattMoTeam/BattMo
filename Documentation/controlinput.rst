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
      
The battery model contains a *Control* sub-model, see the overall model structure :ref:`description<architecture:BattMo
Model Architecture>`. The Control model determine the control *values* and *type*.

At a **given time**, there are in general two control types current or voltage. We have also include a power control

Given a control type, the control *value* can change in time. For the most standard control types, the json interface
can be used. 

The parameters for each of the model are described in the schema, see :ref:`Control Parameters <json:Control Parameters>`.

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
              
