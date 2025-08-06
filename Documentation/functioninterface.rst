=====================
Functional parameters
=====================

Battery models contain functional parameters. The most important one is probably the OCV function, but we have also the
ionic conductivity and diffusion in the electrolyte.

The function definitions are given in the input json file, using one of the following formats:

1. Function name
2. Function formula
3. Tabulated function

The json schema that describes the interface is available :battmofile:`here <Utilities/JsonSchemas/Function.schema.json>`.

The keyword to choose the format is given by :code:`functionFormat`.

The keyword :code:`argumentList` is common to all the formats and is used to describe the list of arguments that will be
sent to the function. The argument names given there are independent to the implementation in itself.

Let us go through each different format separatly. For each format, we start by giving an example.

Function name
=============

.. literalinclude:: exampleCodeSnippets/function_interface_name.json
   :language: json

For a formula, the function format is set to :code:`"named function"`.

With the key :code:`"functionName"`, we give the name of the function. In matlab, this function should then be found in
the function path (more on the in  `matlab documentation <https://se.mathworks.com/help/matlab/ref/addpath.html>`__).

The key :code:`argumentList` is used to describe the list of arguments that will be sent to the function. Note that the
names do not necessarily to the one used in the implementation itself, as we can check in the implementation of
:battmo:`computeOCP_Graphite_Torchio` where the argument there is set to the variable name :code:`theta`.

              
Function formula
================

.. literalinclude:: exampleCodeSnippets/function_interface_string.json
   :language: json

For a formula, the function format is set to :code:`"string expression"`.

The key :code:`argumentList` is used to describe the list of arguments that will be sent to the functions. Again, it is not related to the actual implementation described under the key :code:`expression`.

With the key :code:`"expression"`, we give the expression which, once evaluated, will give the value of the
function.

The key :code:`language` is used to provide the programming language that the expression should be evaluated with. The
default value is :code:`matlab` but we aim at a language agnostic interface and we open for the possibility to include
an expression using an other language, even if it is not supported in the matlab code.

The key :code:`formula` is used to give the expression that we will be run and return the value. The variable names
should be given using the key :code:`variableNames`.

You can find an other example 



   


            
