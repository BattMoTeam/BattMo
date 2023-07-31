# JuliaBridge

The purpose of JuliaBridge is to act as a way to run BattMo.jl from matlab, a sort of bridge between Matlab and Julia.

To run BattMo.jl from matlab we will use the package [MATDaemon.jl](https://github.com/jondeuce/MATDaemon.jl). This package works by creating a persistent Julia server using Daemon.jl and then sending information between matlab and julia through reading and writing .mat files. I have found that writing into .mat files in Matlab and then reading them in  Julia generally goes smoothly. The other way does however tend to be more tricky. MATDaemon uses the MAT.jl package which requires the input to be a dict when writing into .mat files. This can be quite problematic when using custom types and ordered objects. I have therefore made some small tweaks to this package where I instead use .json files when sending information from Julia to Matlab. This implementation will therefore use [Erasdna/MATDaemon](https://github.com/Erasdna/MATDaemon.jl).

JuliaBridge also only works using the [matlab-bridge](https://github.com/BattMoTeam/BattMo.jl/tree/matlab-bridge) branch of BattMo.jl. The important difference here from the master branch is that the mrst_utils.jl files has been modified to allow inputs from either .json or .mat files or directly from Matlab.

## How to use JuliaBridge

To run JuliaBridge we need:

<ul>
    <li> BattMo must be on the Matlab path</li>
    <li> The file <code>jlcall.m</code> needs to be on the Matlab path </li>
    <li> The <code>matlab-bridge</code> branch of BattMo.jl needs to be installed </li>
</ul>

Examples can be found in the Examples folder

## The <code> jlcall </code> function

A complete guide to the <code> jlcall </code> function can also be found in the [MATDaemon readme](https://github.com/jondeuce/MATDaemon.jl/blob/master/README.md). I will briefly summarize some of the most important parts of the present implementation:

### Activating an environment

We can use customized code by adding a project using the following option: <br>
<code>
    jlcall('project', 'path/to/MyProject')
</code>

In practice think of this as activating the environment <code>MyProject</code>. We can now use custom modules defined in this project. The Julia_utils/RunFromMatlab repo provides an example with a wrapper method for running <code>run_battery</code>.

### Setup file

We can add a setup file using the following line: <br>
<code>
    jlcall('setup', 'path/to/setup.jl')
</code>
<br>
The setup file is especially useful when managing the communication between matlab and Julia. In particular the setup file inside Julia_utils provides specific methods for serializing custom types into json objects. In the event of a serialization issue another such function should be added in this file.

### Import modules

We can import modules using the following line: <br>
<code>
    jlcall('modules', {'module1', 'module2'})
</code>
<br>
We can think of this as being equivalent to writing <code> import module1 module2 </code> in Julia. Note that these modules can also come from an activated project.

### Restarting the server

We can decided whether the server should remain active or be restarted each time the code is run in matlab using the following line: <br>

<code>
    jlcall('restart', true/false)
</code>

<br>

If restart is set to <code>true</code> the Julia code will not be recompiled. 

### Using the correct MATDaemon



### Running code from julia

Having imported a module we can now call functions from the module by: <br>

<code>
    jlcall('module.MyFunction',{args,kwargs})
</code>

<br>

Note that <code>args</code> should be a cell, while <code>kwargs</code> should be a dict. 

## Possible improvements, ideas

<ul>
    <li> Figure out how you can use Revise on the persistent Julia server to make it easier to modify Julia code</li>
    <li> Enable a non-blocking version of jlcall. Would make it possible to build the matlab model while Julia is compiling </li>
    <li> Somehow be able to open Julia REPL of the active environment</li>
</ul>