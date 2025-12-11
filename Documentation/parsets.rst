==============
Parameter sets
==============

The following parameter sets are available directly through the installation. They are convenient for testing. They are
listed in the directory :battmofile:`ParameterSets<ParameterData/ParameterSets>`.

.. list-table:: Parameter sets
   :header-rows: 1

   * - Authors
     - Year
     - Reference
     - Directory
   * - Chen et al
     - 2020
     - :cite:`Chen2020DevelopmentModels`
     - :battmofile:`Chen2020<ParameterData/ParameterSets/Chen2020>`
   * - Safari et al
     - 2009
     - :cite:`Safari_2009`
     - :battmofile:`Safari2009<ParameterData/ParameterSets/Safari2009>`
   * - Xu et al
     - 2015
     - :cite:`XU2015303`
     - :battmofile:`Xu2015<ParameterData/ParameterSets/Xu2015>`


Otherwise, separate input json files corresponding to different part of a battery model are included in the directory
:battmofile:`ParameterData`. They are organised in sub directories, whose name describe the type of data they
correspond to

  * :battmofile:`BatteryCellParameters<ParameterData/BatteryCellParameters>`
  * :battmofile:`BatteryComponentParameters<ParameterData/BatteryComponentParameters>`
  * :battmofile:`MaterialProperties<ParameterData/MaterialProperties>`
