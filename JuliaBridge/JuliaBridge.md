# JuliaBridge

The purpose of JuliaBridge is to act as an interface between BattMo Matlab and BattMo Julia. 

## How to use JuliaBridge

**Please Note**: JuliaBridge does currently **not work with Matlab 2023a** due to an issue with the <code>system</code> command in Matlab. (We recommend Matlab 2022a)

In order for JuliaBridge to work properly, the following must be done:
<ul>
    <li> The RunFromMatlab project must be instantiated in Julia</li>
    <li> The GenerateModel and JuliaInterface folders must be on the Matlab path </li>
</ul>

The first point will be automatically done by the <code>ServerManager</code> constructor, while the second point can be completed by running <code>startupJuliaBridge.m</code>.

For specific use cases please the files provided in the Examples folder. 


### Additional steps for Windows users:

JuliaBridge uses <code>DeamonMode.jl</code> to launch a persistent Julia server running in the background using the <code>system</code> command in Matlab. Launching a background process in this manner does not work on Windows. The user will therefore have to launch the server themselves. This can be done in the following way:

1. Launch a command prompt window in Windows (NOT the powershell)
1.5. (Recommended) Use the <code>cd</code> command to the place where RunFromMatlab is stored
2. Run the following command (example):
<code> julia --startup-file=no --project=. -e "using Revise, DaemonMode; serve(3000, true, call_stack=true, async=true)" </code>
3. Running this command will block the command prompt. The server will remain active until the window is closed or it is deactivated in any other way. Calls to the server can now be made using the <code>ServerManager</code> class.

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
    <li><b>Run_battery</b>: Runs the <code>run_battery</code> function in Matlab using either preloaded data or data passed as an argument</li>
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
