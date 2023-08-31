import MAT, JSON 
using RunFromMatlab

juliafy_kwargs(xs::Dict{String, Any}) = Pair{Symbol, Any}[Symbol(k) => v for (k, v) in xs]

script = ARGS[1]
CALLABLE = ["-load", "-load_options", "-run_battery"];

if !in(script,CALLABLE)
    
    println("Error: ", script, " is not a valid call option")
    println("Currently available call options are: \n",CALLABLE)
    
else
    
    if script == "-load"
        
        load_file = ARGS[2]
        dat       = MAT.matread(load_file)
        inputType = dat["inputType"]
        kwargs    = juliafy_kwargs(dat["kwargs"])

        if inputType == "Matlab"
            data          = dat["data"]
            use_state_ref = dat["use_state_ref"]
        end
        
        inputFileName = dat["inputFileName"]
        
        if opts["gc"]
            rm(load_file)
        end
        
    elseif script == "-load_options"
        
        load_file = ARGS[2]
        dat       = MAT.matread(load_file)
        opts      = dat["op"]
        
    elseif script == "-run_battery"

        if inputType == "Matlab"
            output = RunFromMatlab.run_battery_from_matlab(inputFileName, data, use_state_ref = use_state_ref; kwargs... )
        elseif inputType == "JSON"
            output = RunFromMatlab.run_battery_from_matlab(inputFileName; kwargs... )
        else
            error("inputType not recognized")
        end
        stringdata = JSON.json(output)
        
        # write the file with the stringdata variable information
        outputFileName = ARGS[2]
        open(outputFileName, "w") do f
            write(f, stringdata)
        end
        
    end
end
