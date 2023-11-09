import BattMo

function run_from_matlab(init::BattMo.MatlabFile; 
                                 kwargs...)
    """ 
        Summary: Wrapper method for running run_battery when launched from matlab. 
        
        Inputs:
            -data: Dict containing fields imported from matlab. 
            -inputFileName: Name of matlab file that generated data
            -use_state_ref: If true a previous simulation run in matlab is used as a starting point
        Output:
            -ret: Dict with fields:
                -states: States of the simulation
                -reports: Reports on each timesteps
                -extra: Miscellaneous other information on the simulation
                -data: Input received by Julia
    """

    # Create input
    # If use_state_ref is true the simulation will try to replicate the reference solution generated in matlab
    # by using the same timesteps, input current and cutoff voltage.
    states, reports, extra, exported = BattMo.run_battery(init; kwargs...);
    # create output
    ret = Dict("states" => states, "reports" => reports, "extra" => extra, "exported" => exported, "matlab states" => convert_state_Matlab(states,extra[:model],extra[:timesteps]))
    return ret
end

function run_from_matlab(init::BattMo.JSONFile; kwargs...)
    """ Json version"""
    states, reports, extra, exported = BattMo.run_battery(init; kwargs...);
    # create output
    ret = Dict("states" => states, "reports" => reports, "extra" => extra, "exported" => exported, "matlab states" => convert_state_Matlab(states,extra[:model],extra[:timesteps]))
    return ret
    
end

