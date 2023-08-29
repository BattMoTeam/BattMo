# Dynamically include setup code, import modules, and finally evaluate the function expression passed by the user.
# User settings are loaded from the input file `MATDaemon.JL_OPTIONS` located in the jlcall.m workspace folder.
# The workspace folder is passed using the environment variable `MATDAEMON_WORKSPACE`.

# MATDaemon must be available
import MATDaemon

let
    # Load jlcall.m input parser results from workspace
    local workspace = ENV["MATDAEMON_WORKSPACE"]
    local opts = MATDaemon.load_options(workspace)

    # Initialize user project environment etc.
    MATDaemon.init_environment(opts)

    # Print environment for debugging
    if opts.debug
        println("* Environment for evaluating Julia expression:")
        println("*   Working dir: $(pwd())")
        println("*   Module: $(@__MODULE__)")
        println("*   Load path: $(LOAD_PATH)")
        println("*   Active project: $(Base.active_project())", "\n")
    end

    # Include setup code
    if !isempty(opts.setup)
        include(opts.setup)
    end

    # Load modules from strings; will fail if not installed
    for mod in opts.modules
        @eval import $(Meta.parse(mod))
    end

    # Parse and evaluate `f` from string
    local f_expr = :(let; $(Meta.parse(opts.f)); end)

    # If not a function call, return a thunk
    if opts.nofun
        f_expr = :(let; $(f_expr); (args...; kwargs...) -> nothing; end)
    end

    if opts.debug
        println("* Generated Julia function expression: ")
        println(string(MATDaemon.MacroTools.prettify(f_expr)), "\n")
    end

    local f = @eval $(f_expr)

    # Call `f`, loading MATLAB input arguments from `opts.infile`
    # and saving Julia outputs to `opts.outfile`
    local output = MATDaemon.jlcall(f, opts)
end
