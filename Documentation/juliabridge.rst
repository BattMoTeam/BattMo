============
Julia bridge
============

:todo:`add explaination on installation, examples and some time comparison`

.. note::

   There is a problem in the newest versions of Matlab we have reported and which prevents the use of
   :code:`JuliaBridge`. For the moment, use Matlab versions prior to 2022a.

|battmo| has a julia implementation for the standard PXD model, called **BattMo.jl**, see  `here
<https://github.com/BattMoTeam/BattMo.jl>`_ .

For small problems, the simulator implemented in Julia is *significantly* faster than the matlab solver (several order
of magnitude of speed-up). The reason is that the assembly of the residual equations and their derivatives can be
optimized in completly different ways in a compiled code such as Julia.

The *Julia Bridge* is a set of functionalities which make possible to run the julia solver and retrieve the results
directly from matlab.

The BattMo.jl package is registered in the General Julia registry and the Julia Bridge uses julia package management
system to install it without further user intervention.

.. note::

   Julia uses *just-in-time* compilation. It means that when we run a function for the first time in a julia session,
   the code is going to be compiled. The compilation may take time (more than 10s). When using Julia Bridge, you will
   notice it at the start. But, by starting a Julia server in the background where the compilation is done only once, Julia
   Bridge makes it possible to use the |battmo| julia solver at expected speed.




