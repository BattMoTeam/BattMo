import BattMo

function run_from_matlab(init::BattMo.BattMoFormattedInput;
                         kwargs...)

    output = BattMo.run_battery(init; kwargs...)

    # create output
    ret = Dict("states"             => output[:states],
               "cellSpecifications" => output[:cellSpecifications],
               "reports"            => output[:reports],
               "inputparams"        => output[:inputparams],
               "matlab states"      => convert_state_Matlab(output[:states],
                                                            output[:extra][:model],
                                                            output[:extra][:timesteps]))

    return ret

end
