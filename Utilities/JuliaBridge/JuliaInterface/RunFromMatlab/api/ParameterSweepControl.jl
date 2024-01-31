using Distributed
@everywhere using RunFromMatlab, Random, BattMo
@everywhere import ParallelDataTransfer
import MAT, JSON

sweep_options = MAT.matread(ARGS[1])
ParallelDataTransfer.sendto(workers(),
                            init=inputobj,
                            experiment=sweep_options["experiment"],
                            kwargs=kwargs,
                            folder=sweep_options["save_folder"])

@everywhere function run_parameter(parameters)

    # NB: Init is modified, but values are replaced in the same fields each time!
    RunFromMatlab.setParameters!(init, parameters, experiment)

    output = RunFromMatlab.run_from_matlab(init; info_level=-1, kwargs ...)

    save_location = folder * '/' *randstring(12) * ".json";
    RunFromMatlab.save_output(output, save_location)

    return Dict("parameters" => parameters, "states_location" => save_location)

end

# Values to sweep over
vals = sweep_options["values"]

# Perform the sweep using pmap
experiment_results = pmap((params) -> run_parameter(params), Iterators.product(vals...))
exp_size           = size(experiment_results)

save_output(Dict("Experiment" => sweep_options["experiment"],
                 "Results" => reshape(experiment_results, (exp_size[2]*exp_size[4], ))),
            sweep_options["save_folder"] * "/" * sweep_options["experiment"] * "_output.json")
