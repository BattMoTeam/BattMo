==================
Battery Geometries
==================

.. automodule:: Battery.BatteryGeometry

Description of the models that are used to setup the geometry


BatteryGenerator1D
------------------

.. image:: img/1dmodel.png
           :width: 80%
                   
.. autoclass:: BatteryGenerator1D
   :members:
      
.. _2dgeometry:
         
BatteryGenerator2D
------------------

.. image:: img/runbattery2d.png
           :width: 80%
                   
.. autoclass:: BatteryGenerator2D
   :members:
         
.. _3dgeometry:
      
BatteryGenerator3D
------------------

.. image:: img/runbattery3d.png
           :width: 80%
                   
.. autoclass:: BatteryGenerator3D
   :members:

.. _jellyroll:
      
SpiralBatteryGenerator
----------------------

.. image:: img/jellyrollmodel.png
           :width: 80%
                   
.. autoclass:: SpiralBatteryGenerator
   :members:

.. _coincell:
      
CoinCellBatteryGenerator
------------------------

.. image:: img/coincell.png
           :width: 80%
                   
.. autoclass:: CoinCellBatteryGenerator
   :members:
      
         
BlockBatteryGenerator
---------------------

not yet documented
         
CoinCellSectorBatteryGenerator
------------------------------

not yet documented

         
FlatBatteryGenerator
--------------------

not yet documented

         
SectorBatteryGenerator
----------------------

not yet documented

         
Base class
----------

This is the base class that gather the methods to setup the different grid. This class will be usefull if you want to
setup your own tailored grid.

.. autoclass:: Battery.BatteryGeometry.BatteryGenerator
   :members:
         
