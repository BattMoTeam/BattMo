
function convert_state_Matlab(states::AbstractVector, model, timesteps)
    ret = Array{Dict{String,Any}}(undef, length(states), 1)
    time = cumsum(timesteps)
    for (i,s) in enumerate(states)
        ret[i] = convert_state_Matlab(s, model)
        ret[i]["time"] = time[i];
    end
    return ret
end

# Modified of the setup_state function in BattMo
function convert_state_Matlab(julia_state, model)

    jsonNames = Dict(
        :NeCc => "NegativeElectrode",
        :NeAm => "NegativeElectrode",
        :PeAm => "PositiveElectrode",        
        :PeCc => "PositiveElectrode",
    )

    function convert_current_collector!(matlab_state, name::Symbol)
        """ initialize values for the current collector"""
        
        if haskey(model.models, name)
            use_cc = true
        else
            use_cc = false
        end
        
        if use_cc
            init = Dict("CurrentCollector" => Dict("phi" => julia_state[name][:Phi]))
            matlab_state[jsonNames[name]]= init
        end
        
    end

    function convert_active_material!(matlab_state, name::Symbol)
        """ initialize values for the active material"""

        jsonName = jsonNames[name]
        
        ccnames = Dict(
            :NeAm => :NeCc,
            :PeAm => :PeCc
        )

        if haskey(model.models, ccnames[name])
            use_cc = true
        else
            use_cc = false
        end

        # initialise NAM

        sys = model[name].system

        init = Dict()
        init2 = Dict("c" => collect(Iterators.flatten(julia_state[name][:Cp])), "cSurface" => julia_state[name][:Cs])
        init["ActiveMaterial"] = Dict("phi" => julia_state[name][:Phi], "SolidDiffusion" => init2)

        matlab_state[jsonName] = init
        
    end

    function convert_electrolyte!(matlab_state)

        init = Dict("phi" => julia_state[:Elyte][:Phi], "c" => julia_state[:Elyte][:C])
        
        matlab_state["Electrolyte"] = init

    end

    function convert_bpp!(matlab_state)
        init = Dict("E" => julia_state[:Control][:Phi], "I" => julia_state[:Control][:Current])
    
        matlab_state["Control"] = init
        
    end

    matlab_state= Dict()
    convert_current_collector!(matlab_state, :NeCc)
    convert_active_material!(matlab_state, :NeAm)
    convert_electrolyte!(matlab_state)
    convert_active_material!(matlab_state, :PeAm)
    convert_current_collector!(matlab_state, :PeCc)
    convert_bpp!(matlab_state)

    return matlab_state
    
end
