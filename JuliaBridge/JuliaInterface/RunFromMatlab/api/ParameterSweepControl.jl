@everywhere using Distributed, RunFromMatlab, Random,BattMo
@everywhere import ParallelDataTransfer
import MAT, JSON 

#Debug!
# dat = MAT.matread("/home/solheim/Documents/SINTEF/BattMo/JuliaBridge/JuliaInterface/RunFromMatlab/example_load.mat")
# inputType = dat["inputType"]
# juliafy_kwargs(xs::Dict{String, Any}) = Pair{Symbol, Any}[Symbol(k) => v for (k, v) in xs]

# kwargs    = juliafy_kwargs(dat["kwargs"])

# inputFileName = dat["inputFileName"]

# if inputType == "Matlab"
#     inputobj = BattMo.MatlabFile(inputFileName, dat["data"], use_state_ref=dat["use_state_ref"])
# elseif inputType == "JSON"
#     inputobj = BattMo.JSONFile("/home/solheim/Documents/SINTEF/BattMo/JuliaBridge/Examples/jsonfiles/p2d_40_jl.json")
# else
#     println("Invalid input type. Input data could not be read correctly")
# end

sweep_options = MAT.matread(ARGS[1])
ParallelDataTransfer.sendto(workers(), init=inputobj, experiment=sweep_options["experiment"], kwargs=kwargs, folder=sweep_options["save_folder"])


@everywhere function run_parameter(parameters)
    RunFromMatlab.setParameters!(init, parameters, experiment) #NB: Init is modified, but values are replaced in the same fields each time!
    output = RunFromMatlab.run_battery_from_matlab(init; info_level=-1, kwargs ...)
    save_location = folder * '/'*randstring(12)*".json";
    RunFromMatlab.save_output(output, save_location)
    return Dict("parameters" => parameters, "states_location" => save_location)
end

vals=sweep_options["values"]
experiment_results = pmap((params) -> run_parameter(params), Iterators.product(vals...))
exp_size = size(experiment_results)
save_output(Dict("Experiment" => sweep_options["experiment"], "Results" => reshape(experiment_results, (exp_size[2]*exp_size[4],))), sweep_options["save_folder"]*"/"*sweep_options["experiment"]*"_output.json")