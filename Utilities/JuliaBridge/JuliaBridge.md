# JuliaBridge

The purpose of JuliaBridge is to act as an interface between BattMo Matlab and BattMo Julia. 

## Table of contents

1. [How to use JuliaBridge](#how-to-use-juliabridge)
    1. [Windows users](#additional-steps-for-windows-users)
2. [How JuliaBridge works](#how-juliabridge-works)
    1. [Step-by-step](#matlab-↔-julia-communication)
    2. [ServerManager](#servermanager)
    3. [DaemonMode](#deamonmode)
    4. [RunFromMatlab](#runfrommatlab)
3. [Application: parameter sweeps](#parameter-sweeps)
4. [Improving JuliaBridge](#possible-improvements)

## How to use JuliaBridge

In order for JuliaBridge to work properly, the following must be done:
<ul>
    <li> The RunFromMatlab project must be instantiated in Julia</li>
    <li> The GenerateModel and JuliaInterface folders must be on the Matlab path </li>
</ul>

The first point will be automatically done by the <code>ServerManager</code> constructor, while the second point is also
automatically taken care of by the setup script <code>statupBattMo</code>

For specific use cases please the files provided in the Examples folder. Julia uses just-in-time compilation. It means
that the first time you run a simulation, the code is going to be compiled and it will take relatively a lot of
time. However, we use a persistent server model so that this is done only once on the server. The subsequent simulations
will not suffer from this compilation slow-down.


### Additional steps for Windows users:

JuliaBridge uses <code>DeamonMode.jl</code> to launch a persistent Julia server running in the background using the
<code>system</code> command in Matlab. Launching a background process in this manner does not work on Windows. The user
will therefore have to launch the server themselves. This can be done in the following way:

1. Launch a command prompt window in Windows (NOT the powershell)
1.5. (Recommended) Use the <code>cd</code> command to move to current directory to the place where RunFromMatlab is stored
2. Run the following command (example):
<code> julia --startup-file=no --project=. -e "using Revise, DaemonMode; serve(3000, true, call_stack=true, async=true)" </code>
3. Running this command will block the command prompt. The server will remain active until the window is closed or it is deactivated in any other way. Calls to the server can now be made using the <code>ServerManager</code> class.

## Example

```matlab
   
    man = ServerManager();
   
    inputFileName = fullfile(battmoDir()    , ...
                             'Examples'     , ...
                             'JsonDataFiles', ...
                             'p2d_40_jl.json')

    man.load('inputType'    , 'JSON', ...
             'inputFileName', inputFileName);
            
    result = man.run();
```

see more in the [online documentation](https://battmoteam.github.io/BattMo/juliabridge.html ) 


## How JuliaBridge works

### Matlab &harr; Julia communication
 
In concrete terms JuliaBridge communicates between Matlab and Julia in the following way:

1. Generate a model, initial state and parameters (and optionally a reference solution) in Matlab. This data is written in a file
2. Launch a DaemonMode server
3. Call Julia scripts contained in RunFromMatlab/api/DaemonHandler.jl on the DaemonMode server
4. The Julia code reads the file containing data from Matlab and performs a simulation. The output of this simulation is saved to a file.
5. The contents of the file created by Julia is read in Matlab for post processing.

### ServerManager

The <code>ServerManager</code> class in Matlab has two main roles:
<ul>
    <li> Launch and shut down the DeamonMode server as necessary</li>
    <li> Act as an interface for running Julia scripts</li>
</ul>

#### Notable class methods:
<ul>
    <li> <b>Constructor</b>: Instantiates the RunFromMatlab project and launches a server</li>
    <li> <b>Load</b>: Loads data to the server. If <code>shared=true</code> variables can be accessed across calls to the server. It is recommended to use the function to minimize number of times data is transferred by reading and writing files</li>
    <li><b>Run</b>: Use the <code>run</code> function in Matlab with either preloaded data or data passed as an argument</li>
</ul>

**Please Note**: The DeamonMode server exists as an entity separate from any <code>ServerManager</code> object. Deleting the object will not shut down the server. However, the <code>shutdown</code> method provides a way to close the server. 

### DeamonMode

Source: https://github.com/dmolina/DaemonMode.jl

Below we summarize the two main features of DaemonMode used by JuliaBridge-

#### Launching the server

A DeamonMode server can be launched from the shell using the following command:
<code>julia <options\> "using DeamonMode; serve(port, shared, call_stack, async)"</code>

With arguments:
<ul>
    <li> <b>port</b>: The id of the server, default 3000</li>
    <li> <b>shared</b>: If true, variables will be shared across calls to the server</li>
    <li> <b> call_stack </b>: Displays entire call stack if the Julia code fails to execute</li>
    <li> <b> async </b>: Allows running clients in parallel </li>
</ul>

#### Running a script

We can run a script on the DeamonMode server using the following command:
<code>julia <options\> "using DeamonMode; runargs(port)" program.jl <arguments\></code>

If several servers are active, the port id should be specified. Default is 3000.

### RunFromMatlab

#### DeamonHandler

The DeamonHandler.jl file contains a variety of scripts which can be called by <code>ServerManager</code>. The idea behind this file is the DeamonHandler.jl should be called with a flag that indicates which process to run and then a range of inputs such as the name of the file from which should be loaded or written.

#### src

RunFromMatlab is meant to act as a wrapper for running functions from BattMo. Functionality can easily be extended by adding functions as necessary.

## Parameter sweeps

An important motivation for an interface between BattMo Julia and BattMo Matlab is to be able to quickly perform parameter of battery configurations. The <code>ParameterSweepControl.jl</code> script provides a template for how parameter sweeps can be run using the <code>Distributed</code> package in Julia. 

### API

#### Running a sweep

Parameter sweeps can be launched using the <code>sweep</code> method in <code>ServerManager</code>. This functions takes the following inputs:
<ul>
    <li> <b>Experiment</b>: The experiment tag defines which parameters are to be iterated over. This tag will be used by <code>ParameterSweepControl.jl</code> set the correct values. This tag must correspond to one of the options in <code>src/parameter_sweep_utils.jl</code></li>
    <li><b>Values</b>: Values of the parameters which are to be iterated over. The parameter sweep will take all combinations of the input vectors</li>
    <li><b>Name</b>: Name tag for the experiment being run. Meant to be used as an identifier. A folder will be created with this name and this is were the outputs from the simulations will be stored</li>
</ul>

The parameter input options will be saved as a <code>.mat</code> in the <code>name</code> folder inside the server folder. This file is then read by Julia which launches the specified amount of processors (see the <code>procs</code> option in the <code>ServerManager</code> constructor). The <code>ParameterSweepControl.jl</code> is launched from Matlab as a <b>detached</b> programme. For each parameter the output will be saved in a <code>.json</code> file with a randomly generated tag inside the <code>name</code> folder. Once all runs are complete the <code>/name/\<experiment name\>_output.json</code> file will be created. This file contains the file location of the experiment corresponding to each parameter.

#### Collecting results

Because the paramter sweeps run in the background it can be hard to determine when the simulation is complete. The <code>sweep</code> method in <code>ServerManager</code> returns a [Matlab Futures](https://se.mathworks.com/help/matlab/ref/parallel.future.html) which runs in the background, continously checking if <code>/name/\<experiment name\>_output.json</code> exists. When this file exists the Futures will read it. The <code>collect_results</code> method is a wrapper that simplifies collecting results based on this object.

## Making improvements

JuliaBridge is NOT a generic interface (see rather something like [MATDaemon.jl](https://github.com/jondeuce/MATDaemon.jl)). However, it is easily extendable. Adding functionalities can be done in the following steps:
<ul>
    <li>Add or adapt a method in <code>ServerManager</code></li>
    <li>Add a new script option in <code>DeamonHandler.jl</code></li>
    <li>Build an appropritate call to the Daemon. Please take a look at the <code>DeamonCall</code> method. Set <code>debug</code> to <code>true</code> to get some examples.</li>
</ul>

The idea is the JuliaBridge should evolve as the demands on BattMo evolve.

### Parameter sweeps

What follows is a (non-exhaustive) list of things that could be done to improve the parameter sweep functionality in JuliaBridge.
<ul>
    <li>Create unique name tags for experiment runs. Use functions like <code>md5</code> to create names</li>
    <li>Add objects/checks to ensure that experiments and parameter values are consistent. (Example: Matlab class or other procedural definition of a "valid" experiment)</li>
    <li>Reuse simulation objects across parameters. Currently the simulation object is recreated on each run which is inefficient.</li>
    <li>Add option to provide file to specify parameter tests + specify filenames of each parameter run</li>
    
</ul>
