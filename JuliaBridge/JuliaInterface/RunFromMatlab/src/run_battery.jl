import BattMo

function run_battery_from_matlab(exported::Dict{String,Any}, 
                                 source_file::String, 
                                 use_reference::Bool; 
                                 kwarg...)
    """ 
        Summary: Wrapper method for running run_battery when launched from matlab. 
        
        Inputs:
            -exported: Dict containing fields imported from matlab. 
            -source_file: Name of matlab file that generated exported
            -use_reference: If true run_battery will use a previous simulation run in matlab as a starting point
        Output:
            -ret: Dict with fields:
                -states: States of the simulation
                -reports: Reports on each timesteps
                -extra: Miscellaneous other information on the simulation
                -exported: Input received by Julia
    """

    # Create input
    # If state_ref is true the simulation will try to replicate the reference solution generated in matlab
    # by using the same timesteps, input current and cutoff voltage.
    init = BattMo.MatlabFile(source_file, exported, state_ref=use_reference)
    states, reports, extra, exported = BattMo.run_battery(init; kwarg...);
    # create output
    ret = Dict("states" => states, "reports" => reports, "extra" => extra, "exported" => exported)
    return ret
end
