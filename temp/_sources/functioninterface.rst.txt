=====================
Functional parameters
=====================

Battery models contain functional parameters. The most important one is probably the OCV function, but we have also the
ionic conductivity and diffusion in the electrolyte. The function definitions can be entered in the input json
file using one of the following format

1. Function name
2. Function formula
3. Tabulated function

The json schema that describes the interface is available :battmofile:`here <Utilities/JsonSchemas/Function.schema.json>`.

The keyword to choose the format is given by :code:`functionFormat`.

Let us go through each different format separatly. We start with en example for each

Function name
=============

.. literalinclude:: exampleCodeSnippets/function_interface_name.json
   :language: json

For a formula, the function format is :code:`"named function"`. With the key :code:`"functionName"`, we give the name of
the function. In matlab, this function should then be found in the function path (see `matlab doc <https://se.mathworks.com/help/matlab/ref/addpath.html>`__). The key 

              
Function formula
================

.. literalinclude:: exampleCodeSnippets/function_interface_string.json
   :language: json

For a formula, the function format is :code:`"string expression"`. With the key :code:`"function"`, we give the
expression which, once evaluated, will give the value of the function. With the key :code:`"argumentList"`, we declare
the arguments of the function.





   


            
