

function convert_state_Matlab(states::AbstractVector,model)
    ret = Array{Dict{String,Any}}(undef,length(states),1)
    for (i,s) in enumerate(states)
        ret[i]=convert_state_Matlab(s,model)
    end
    return ret
end

function convert_state_Matlab(julia_state,model)

    jsonNames = Dict(
        :CC  => "NegativeElectrode",
        :NAM => "NegativeElectrode",
        :PAM => "PositiveElectrode",        
        :PP  => "PositiveElectrode",
    )

    function convert_current_collector!(matlab_state, name::Symbol)
        """ initialize values for the current collector"""
        
        if haskey(model.models, name)
            use_cc = true
        else
            use_cc = false
        end
        
        if use_cc
            init =Dict("CurrentCollector" => Dict("phi" => julia_state[name][:Phi]))
            matlab_state[jsonNames[name]]= init
        end
        
    end

    function convert_active_material!(matlab_state, name::Symbol)
        """ initialize values for the active material"""

        jsonName = jsonNames[name]
        
        ccnames = Dict(
            :NAM => :CC,
            :PAM => :PP,
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
        init["ActiveMaterial"]=Dict("phi" => julia_state[name][:Phi], "SolidDiffusion" => init2)


        matlab_state[jsonName] = init
        
    end

    function convert_electrolyte!(matlab_state)

        init = Dict("phi" => julia_state[:ELYTE][:Phi], "c" => julia_state[:ELYTE][:C])
        
        matlab_state["Electrolyte"] = init

    end

    function convert_bpp!(matlab_state)
        init = Dict("E" => julia_state[:BPP][:Phi], "I" => julia_state[:BPP][:Current])
    
        matlab_state["Control"] = init
        
    end

    matlab_state= Dict()
    convert_current_collector!(matlab_state, :CC)
    convert_active_material!(matlab_state, :NAM)
    convert_electrolyte!(matlab_state)
    convert_active_material!(matlab_state, :PAM)
    convert_current_collector!(matlab_state, :PP)
    convert_bpp!(matlab_state)

    return matlab_state
end