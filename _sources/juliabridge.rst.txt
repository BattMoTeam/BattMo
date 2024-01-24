===================
BattMo Julia bridge
===================

Introduction
============


|battmo| has a julia implementation for the standard PXD model: **BattMo.jl**. This package is registered in the General
Julia registry and available from `github <https://github.com/BattMoTeam/BattMo.jl>`_ . For small problems, the
simulator implemented in Julia is *significantly* faster than the matlab solver (several order of magnitude of
speed-up). The reason is that the assembly of the residual equations and their derivatives can be optimized in completly
different ways in a compiled code such as Julia. Our julia implementation relies on the package `Jutul
<https://github.com/sintefmath/Jutul.jl>`_.

We have implemented in matlab a *Julia Bridge* which provides a set of functionalities to run the julia solver and
retrieve the results directly from within matlab. The Julia Bridge uses julia package management system to install
BattMo.jl without further user intervention.




Start Server
============

To run simulation with Julia, we use a server. To communicate with the server we create a server manager from Matlab

.. code:: matlab

   man = ServerManager();          

This command starts the Julia server and return a manager :code:`man` which we use to pass data and run simulations.

.. important::

   Julia uses *just-in-time* compilation. It means that when we run a function for the first time in a julia session,
   the code is going to be compiled. The compilation may take time (more than 10s). When using Julia Bridge, you will
   notice it at the start.

   But, **by starting a Julia server in the background** where the compilation is done only once, Julia
   Bridge makes it possible to use the |battmo| julia solver at expected speed.

   It is important to start the server only once and reuse it for further simulation, as described below

.. attention::

   We could not start a persistent process from Matlab in Windows. The step must then be done manually, as follows

   * Launch a command prompt window in Windows (NOT the powershell)
   * Run the following command:

     .. code:: matlab

        julia --startup-file=no --project=\path\to\RunMatlab\directory -e "using Revise, DaemonMode; serve(3000, true, call_stack=true, async=true)"

     Replace :code:`/path/to/RunMatlab/directory` with the path to the
     :battmofile:`RunFromMatlab<Utilities/JuliaBridge/JuliaInterface/RunFromMatlab/>` directory, that is
     :code:`Utilities\\JuliaBridge\\JuliaInterface\\RunFromMatlab` if your current directory is BattMo root installation
     directory.
     
   * Running this command will block the command prompt. The server will remain active until the window is closed or it
     is deactivated in any other way. Calls to the server can now be made using the :code:`ServerManager` class.

Send simulation parameters
==========================
   
We pass data to the server using the :code:`load` method.

.. code:: matlab

   inputFileName = fullfile(battmoDir(), ...
                            'Examples' , ...
                            'julia'    , ...
                            'jsonfiles', ...
                            'p2d_40_jl_ud.json')

   man.load('inputType'    , 'JSON', ...
            'inputFileName', inputFileName);
          
The :code:`inputType` can be either :code:`JSON` as in example above or :code:`Matlab`. For JSON input type, the input
is fully processed in julia but we only support for the moment 1D geometry in the Julia code. With the Matlab input
type, we can use any of the :ref:`geometries<geometryinput:Battery Geometries>` available in |battmo| and pass it as
an input to the julia solver.

The json input :code:`p2d_40_jl_ud.json` contains functional parameters. For example, we have

.. code:: json

   { "Electrolyte": {
        "ionicConductivity": {
           "type": "function",
           "function": "1e-4*c*((-10.5 + 0.668e-3*c + 0.494e-6*c^2) + (0.074 - 1.78e-5*c - 8.86e-10*c^2)*T + (-6.96e-5 + 2.80e-8*c)*T^2)^2",
           "argumentlist": [
             "c",
             "T"]}} 

The full listing is available :battmofile:`here<Examples/julia/jsonfiles/p2d_40_jl_ud.json>`. For the :code:`function`
property, the string that is given to compute the corresponding value (ionic conductivity in the electrolyte in the
snippet above) should be written with a Julia syntax, as it is passed directly to the julia solver. This should not be a
big issue for Matlab users because the Julia syntax is very close to Matlab for such arithmetic expressions. We plan to
implement json support for tabulated data. Tabulated data give more flexibility and could be used both in Matlab and
Julia.

Run the simulation
==================

Finally, we run the simulation

.. code:: matlab

   result = man.run();

Post process the output
=======================
   
The matlab structure contains the simulation output, which can be processed in Matlab. For example,

.. code:: matlab

   voltage = cellfun(@(x) x.Phi, {result.states.BPP});
   time    = cumsum(result.extra.timesteps);
   plot(time/hour, voltage, "BattMo Julia", LineWidth = 2)
   legend
   grid on
   xlabel('Time / h')
   ylabel('Voltage / h')

.. figure:: img/juliarun.png
   :target: _images/juliarun.png
   :width: 70%
   :align: center









