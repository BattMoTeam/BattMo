# JuliaBridge

The purpose of JuliaBridge is to act as a way to run BattMo.jl from BattMo, a sort of bridge between Matlab and Julia.

To run BattMo.jl from matlab we will use the package [MATDaemon.jl](https://github.com/jondeuce/MATDaemon.jl). This package works by creating a persistent Julia server using Daemon.jl and then sending information between matlab and julia through reading and writing .mat files. I have found that writing into .mat files in Matlab and then reading them in  Julia generally goes smoothly. The other way does however tend to be more tricky. MATDaemon uses the MAT.jl package which requires the input to be a dict when writing into .mat files. This can be quite problematic when using custom types and ordered objects. I have therefore made some small tweaks to this package where I instead use .json files when sending information from Julia to Matlab. This implementation will therefore use [Erasdna/MATDaemon](https://github.com/Erasdna/MATDaemon.jl).

JuliaBridge also only works using the [matlab-bridge](https://github.com/BattMoTeam/BattMo.jl/tree/matlab-bridge) branch of BattMo.jl. The important difference here from the master branch is that the mrst_utils.jl files has been modified to allow inputs from either .json or .mat files or directly from Matlab.

## The <code> jlcall </code> function

Define environment

Setup file

Modules to import

Restarting


