=====
Units
=====

BattMo is designed to use SI units in all its calculations. However, to make things easier for the user we provide utilities for converting various common units to SI. The full list of common units can be found here: `Supported Units <https://github.com/SINTEF-AppliedCompSci/MRST/tree/main/core/utils/units>`_ .

The list also includes common unit prefixes such as milli, kilo etc.

Converting units in BattMo
--------------------------

Some examples of using unit conversions in BattMo are shown below.

.. code:: matlab

  >> cElectrolyte   = 5e-1*mol/litre

  cElectrolyte =

    500.0000

  >> I = 0.62*ampere/(1*centi*meter)^3

  I =

    6.2000e+05


Units and JSON input
--------------------

When a value is given in a json input file without a corresponding unit string, or a value is given directly in BattMo, the value is assumed to be in SI units. 

Values can be given in the json input as a value and unit pair as shown below. In this case the value will be converted to SI before use.

.. code:: JSON

        "density": {
          "value": 1.1,
          "unit": "gram/((centi*meter)^3)"
        }

Unit conversion is demonstrated below. 

.. note:: This conversion happens automatically inside :battmo:`BatteryInputParams`, the code below is just for demonstration purposes.

.. code:: matlab

  >> json = '{"density": {"value": 1.1,"unit": "gram/((centi*meter)^3)"}}';

  >> jsonstruct = jsondecode(json);

  >> jsonstruct.density

  ans = 

    struct with fields:

      value: 1.1000
      unit: 'gram/((centi*meter)^3)'

  >> [val,isconverted] = convertUnitBattMo(jsonstruct.density)


  val =

    1.1000e+03


  isconverted =

    logical

    1

.. note:: If a string is given which is not supported, BattMo will return an error. However, care must be taken as a string which is supported but used in error will not throw an error.



