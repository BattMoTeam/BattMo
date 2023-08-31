import BattMo

function run_battery_from_matlab(inputFileName::String; kwargs...)
    init = BattMo.JSONFile(inputFileName)
    states, reports, extra, exported = BattMo.run_battery(init; kwarg...);
    # create output
    ret = Dict("states" => states, "reports" => reports, "extra" => extra, "exported" => exported)
    return ret
end

function run_battery_from_matlab(inputFileName::String,
                                 data::Dict{String,Any}; 
                                 use_state_ref::Bool = false, 
                                 kwarg...)
    """ 
        Summary: Wrapper method for running run_battery when launched from matlab. 
        
        Inputs:
            -data: Dict containing fields imported from matlab. 
            -inputFileName: Name of matlab file that generated data
            -use_state_ref: If true run_battery will use a previous simulation run in matlab as a starting point
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
    init = BattMo.MatlabFile(inputFileName, data, use_state_ref = use_state_ref)
    states, reports, extra, exported = BattMo.run_battery(init; kwarg...);
    # create output
    ret = Dict("states" => states, "reports" => reports, "extra" => extra, "exported" => exported)
    return ret
end
