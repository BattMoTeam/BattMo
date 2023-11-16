:orphan:
   
===============
Getting Started
===============
           
Here, we present some guidelines and tutorials to get started with and give a short description of the example scripts
that are available in |battmo|

* |battmo| is straightforward to run from a single json, as described in this :ref:`Example<runFirstJsonScript>`. 

* In this :ref:`json example <runJsonScript>`, we show how to combine different json files from different sources
  
* The :ref:`battMoTutorial <battMoTutorial>` describes the different step of a simulation when it is set manually beyond the options that
  are given by the json input format. They consist of

* In this :ref:`example <runBatteryP2D>`, more options are used such as some settings of the non-linear solver.

* More Examples

  :todo:`fix this section`
        
  In the :code:`Examples` directory examples, the following examples can be found. They should run out of the box (write
  the name of script after :ref:`installing BattMo<installation>`). If not, tell us!
  
  - :code:`runBatteryP3D` : 2D example using :ref:`2D model geometry<2dgeometry>`
  - :code:`runBatteryP4D` : 3D example using :ref:`3D model geometry<3dgeometry>`
  - :code:`runChen2020` : Example using Chen data including a comparison with `pybamm <https://www.pybamm.org/>`_
  - :code:`runCR` : Example running for a coin cell
  - :code:`runGittTest` : Example running a Gitt test
  - :code:`runJellyRoll` : Example running a :ref:`jelly roll geometry <jellyroll>`


.. toctree::
   :hidden:

   runFirstJsonScript
   publishedExamples/runJsonScript   
   publishedExamples/battMoTutorial
   publishedExamples/runBatteryP2D
   publishedExamples/runSEIActiveMaterial
   publishedExamples/runSiliconGraphiteBattery
      
