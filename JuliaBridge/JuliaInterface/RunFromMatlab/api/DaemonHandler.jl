import MAT, JSON 
using RunFromMatlab

juliafy_kwargs(xs::Dict{String, Any}) = Pair{Symbol, Any}[Symbol(k) => v for (k, v) in xs]

script = ARGS[1]
CALLABLE = ["-load","-load_options","-run_battery"];

if !in(script,CALLABLE)
    println("Error: ", script, " is not a valid call option")
    println("Currently available call options are: \n",CALLABLE)
else
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
        #Convert argument to Float
        use_state_ref = Bool(parse(Float64,ARGS[3]))
        println(use_state_ref)
        output=RunFromMatlab.run_battery_from_matlab(data, load_file, use_state_ref ; kwargs ... )
        stringdata = JSON.json(output)
        # write the file with the stringdata variable information
        open(save_file, "w") do f
            write(f, stringdata)
        end
    end
end