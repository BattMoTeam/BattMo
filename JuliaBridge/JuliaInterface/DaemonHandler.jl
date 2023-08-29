import MAT, JSON, RunFromMatlab

include((@__DIR__)*"/JuliaBridge/JuliaInterface/JSON_lower_overloads.jl")

juliafy_kwargs(xs::Dict{String, Any}) = Pair{Symbol, Any}[Symbol(k) => v for (k, v) in xs]

script = ARGS[1]
CALLABLE = ["-load","-load_options","-run_battery"];

if script=="-load"
    load_file = ARGS[2]
    dat=MAT.matread(load_file)
    data = dat["data"]
    kwargs = juliafy_kwargs(dat["kwargs"])
    if opts["debug"]
        rm(load_file)
    end
elseif script=="-load_options"
    load_file=ARGS[2]
    dat=MAT.matread(load_file)
    opts=dat["op"]
elseif script ==  "-run_battery"
    save_file = ARGS[2]
    output=RunFromMatlab.run_battery_from_matlab(data, load_file, true; kwargs ... )
    stringdata = JSON.json(output)
    # write the file with the stringdata variable information
    open(save_file, "w") do f
        write(f, stringdata)
    end
end