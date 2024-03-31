import MAT, JSON, BattMo 
using RunFromMatlab

juliafy_kwargs(xs::Dict{String, Any}) = Pair{Symbol, Any}[Symbol(k) => v for (k, v) in xs]

script = ARGS[1]
CALLABLE = ["-load", "-load_options", "-run", "-matlab-sweep"]
# 
if !in(script, CALLABLE)
    
    println("Error: ", script, " is not a valid call option")
    println("Currently available call options are: \n", CALLABLE)
    
else
    
    if script == "-load"

        load_file = ARGS[2]
        dat       = MAT.matread(load_file)

        # Expected fields in dat
        # - "inputType"     : String, either "Matlab" or "Json"
        # - "inputFileName" : String, filename
        # - "kwargs"        : Struct, which is passed as key-valued argument command run_from_matlab
        #
        # Those are directly copied in the current workspace on the server and the object
        # inputobj is constructed, ready to be run by the command flag -run (see below)
        #
        
        inputType = dat["inputType"]
        
        if !isempty(dat["kwargs"])
            kwargs = juliafy_kwargs(dat["kwargs"])
        else
            kwargs = ()
        end

        inputFileName = dat["inputFileName"]

        if inputType == "Matlab"
            inputobj = BattMo.MatlabFile(inputFileName, dat["data"], use_state_ref=dat["use_state_ref"])
        elseif inputType == "JSON"
            inputobj = BattMo.JSONFile(inputFileName)
        else
            println("Invalid input type. Input data could not be read correctly")
        end
        
        if opts["gc"]
            rm(load_file)
        end
        
    elseif script == "-load_options"

        load_file = ARGS[2]
        dat       = MAT.matread(load_file)

        # The field "op" in dat is loaded in current workspace on the server, with the name opts
        
        opts = dat["op"]
        
    elseif script == "-run"

        output = RunFromMatlab.run_from_matlab(inputobj; kwargs...)
        stringdata = JSON.json(output)
        
        # write the file with the string data variable information
        outputFileName = ARGS[2]
        save_output(output, outputFileName);
        
    end

end


