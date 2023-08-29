function varargout = jlcall(varargin)
% JLCALL Call Julia from MATLAB using a Julia daemon launched by <a href="matlab: web('https://github.com/dmolina/DaemonMode.jl')">DaemonMode.jl</a>.
% Find the full documentation at the <a href="matlab: web('https://github.com/jondeuce/MATDaemon.jl')">MATDaemon.jl GitHub repository</a>.
% 
% ## Quickstart
% 
% Use the MATLAB function JLCALL to call Julia from MATLAB:
% 
%   >> JLCALL('sort', {rand(2,5)}, struct('dims', int64(2)))
%   
%   ans =
%   
%       0.1270    0.2785    0.6324    0.8147    0.9575
%       0.0975    0.5469    0.9058    0.9134    0.9649
% 
% The positional arguments passed to JLCALL are:
% 1. The Julia function to call, given as a MATLAB char array. This can be any Julia expression which evaluates to a function.
%    For example, 'a=2; b=3; x -> a*x+b'. For convenience, the empty string '' is interpreted as '(args...; kwargs...) -> nothing', returning nothing for any inputs.
%    NOTE: expressions are wrapped in a let block and evaluated in the global scope
% 2. Positional arguments, given as a MATLAB cell array. For example, args = {arg1, arg2, ...}
% 3. Keyword arguments, given as a MATLAB struct. For example, kwargs = struct('key1', value1, 'key2', value2, ...)
% 
% The first time JLCALL is invoked:
% 1. MATDaemon.jl will be installed into a local Julia project, if one does not already exist. By default, a folder .jlcall is created in the same folder as JLCALL
% 2. A Julia server will be started in the background using DaemonMode.jl
% 
% All subsequent calls to Julia are run on the Julia server.
% The server will be automatically killed when MATLAB exits.
% 
% ### Restarting the Julia server
% 
% In the event that the Julia server reaches an undesired state, the server can be restarted by passing the 'restart' flag with value true:
% 
%   >> JLCALL('', 'restart', true) % restarts the Julia server and returns nothing
% 
% Similarly, one can shutdown the Julia server without restarting it:
% 
%   >> JLCALL('', 'shutdown', true) % shuts down the Julia server and returns nothing
% 
% ### Setting up the Julia environment
% 
% Before calling Julia functions, it may be necessary or convenient to first set up the Julia environment.
% For example, one may wish to activate a local project environment, run setup scripts, import modules for later use,
% or set the number of threads for running multithreaded code.
% 
% This setup can be conveniently executed at the start of your MATLAB script with a single call to JLCALL as follows:
% 
%   >> JLCALL('', ...
%      'project', '/path/to/MyProject', ... % activate a local Julia Project
%      'setup', '/path/to/setup.jl', ... % run a setup script to load some custom Julia code
%      'modules', {'MyProject', 'LinearAlgebra', 'Statistics'}, ... % load a custom module and some modules from Base Julia
%      'threads', 'auto', ... % use the default number of Julia threads
%      'restart', true ... % start a fresh Julia server environment
%      )
% 
% See the corresponding sections below for more details about these flags.
% 
% ### Julia multithreading
% 
% The number of threads used by the Julia server can be set using the 'threads' flag:
% 
%   >> JLCALL('() -> Threads.nthreads()', 'threads', 8, 'restart', true)
%   
%   ans =
%   
%     int64
%   
%      8
% 
% The default value for 'threads' is 'auto', deferring to Julia to choose the number of threads.
% 
% NOTE: Julia cannot change the number of threads at runtime.
% In order for the 'threads' flag to take effect, the server must be restarted.
% 
% ### Loading modules
% 
% Julia modules can be loaded and used:
% 
%   >> JLCALL('LinearAlgebra.norm', {[3.0; 4.0]}, 'modules', {'LinearAlgebra'})
%   
%   ans =
%   
%        5
% 
% NOTE: modules are loaded using import, not using. Module symbols must therefore be fully qualified, e.g. LinearAlgebra.norm in the above example as opposed to norm.
% 
% ### Persistent environments
% 
% By default, previously loaded Julia code is available on subsequent calls to JLCALL.
% For example, following the above call to LinearAlgebra.norm, the LinearAlgebra.det function can be called without loading LinearAlgebra again:
% 
%   >> JLCALL('LinearAlgebra.det', {[1.0 2.0; 3.0 4.0]})
%   
%   ans =
%   
%       -2
% 
% ### Unique environments
% 
% Set the 'shared' flag to false in order to evaluate each Julia call in a separate namespace on the Julia server:
% 
%   % Restart the server, setting 'shared' to false
%   >> JLCALL('LinearAlgebra.norm', {[3.0; 4.0]}, 'modules', {'LinearAlgebra'}, 'restart', true, 'shared', false)
%   
%   ans =
%   
%        5
%   
%   % This call now errors, despite the above command loading the LinearAlgebra module, as LinearAlgebra.norm is evaluated in a new namespace
%   >> JLCALL('LinearAlgebra.norm', {[3.0; 4.0]}, 'shared', false)
%   ERROR: LoadError: UndefVarError: LinearAlgebra not defined
%   Stacktrace:
%    ...
% 
% ### Unique Julia instances
% 
% Instead of running Julia code on a persistent Julia server, unique Julia instances can be launched for each call to JLCALL by passing the 'server' flag with value false.
% 
% NOTE: this may cause significant overhead when repeatedly calling JLCALL due to Julia package precompilation and loading:
% 
%   >> tic; JLCALL('x -> sum(abs2, x)', {1:5}, 'server', false); toc
%   Elapsed time is 4.181178 seconds. % call unique Julia instance
%   
%   >> tic; JLCALL('x -> sum(abs2, x)', {1:5}, 'restart', true); toc
%   Elapsed time is 5.046929 seconds. % re-initialize Julia server
%   
%   >> tic; JLCALL('x -> sum(abs2, x)', {1:5}); toc
%   Elapsed time is 0.267088 seconds. % call server; significantly faster
% 
% ### Loading code from a local project
% 
% Code from a local Julia project can be loaded and called:
% 
%   >> JLCALL('MyProject.my_function', args, kwargs, ...
%       'project', '/path/to/MyProject', ...
%       'modules', {'MyProject'})
% 
% NOTE: the string passed via the 'project' flag is simply forwarded to Pkg.activate; it is the user's responsibility to ensure that the project's dependencies have been installed.
% 
% ### Loading setup code
% 
% Julia functions may require or return types which cannot be directly passed from or loaded into MATLAB.
% For example, suppose one would like to query Base.VERSION.
% Naively calling JLCALL('() -> Base.VERSION') would fail, as typeof(Base.VERSION) is not a String but a VersionNumber.
% 
% One possible remedy is to define a wrapper function in a Julia script:
% 
%   # setup.jl
%   julia_version() = string(Base.VERSION)
% 
% Then, use the 'setup' flag to pass the above script to JLCALL:
% 
%   >> JLCALL('julia_version', 'setup', '/path/to/setup.jl')
%   
%   ans =
%   
%       '1.6.1'
% 
% In this case, JLCALL('() -> string(Base.VERSION)') would work just as well.
% In general, however, interfacing with complex Julia libraries using MATLAB types may be nontrivial, and the 'setup' flag allows for the execution of arbitrary setup code.
% 
% NOTE: the setup script is loaded into the global scope using include; when using persistent environments, symbols defined in the setup script will be available on subsequent calls to JLCALL.
% 
% ### Handling Julia outputs
% 
% Output(s) from Julia are returned using the MATLAB cell array varargout, MATLAB's variable-length list of output arguments.
% A helper function MATDaemon.matlabify is used to convert Julia values into MATLAB-compatible values.
% Specifically, the following rules are used to populate varargout with the Julia output y:
% 
% 1. If y::Nothing, then varargout = {} and no outputs are returned to MATLAB
% 2. If y::Tuple, then length(y) outputs are returned, with varargout{i} given by matlabify(y[i])
% 3. Otherwise, one output is returned with varargout{1} given by matlabify(y)
% 
% The following matlabify methods are defined by default:
% 
%   matlabify(x) = x # default fallback
%   matlabify(::Union{Nothing, Missing}) = zeros(0,0) # equivalent to MATLAB's []
%   matlabify(x::Symbol) = string(x)
%   matlabify(xs::Tuple) = Any[matlabify(x) for x in xs] # matlabify values
%   matlabify(xs::Union{<:AbstractDict, <:NamedTuple, <:Base.Iterators.Pairs}) = Dict{String, Any}(string(k) => matlabify(v) for (k, v) in pairs(xs)) # convert keys to strings and matlabify values
% 
% NOTE: MATLAB cell and struct types correspond to Array{Any} and Dict{String, Any} in Julia.
% 
% Conversion via matlabify can easily be extended to additional types.
% Returning to the example from the above section, we can define a matlabify method for Base.VersionNumber:
% 
%   # setup.jl
%   MATDaemon.matlabify(v::Base.VersionNumber) = string(v)
% 
% Now, the return type will be automatically converted:
% 
%   >> JLCALL('() -> Base.VERSION', 'setup', '/path/to/setup.jl')
%   
%   ans =
%   
%       '1.6.1'
% 
% ### Performance
% 
% MATLAB inputs and Julia ouputs are passed back and forth between MATLAB and the DaemonMode.jl server by writing to temporary .mat files.
% The location of these files can be configured with the 'infile' and 'outfile' flags, respectively.
% Pointing these files to a ram-backed file system is recommended when possible (for example, the /tmp folder on Linux is usually ram-backed), as read/write speed will likely improve.
% This is now the default; 'infile' and 'outfile' are created via the MATLAB tempname function (thanks to @mauro3 for this tip).
% 
% Nevertheless, this naturally leads to some overhead when calling Julia, particularly when the MATLAB inputs and/or Julia outputs have large memory footprints.
% It is therefore not recommended to use JLCALL in performance critical loops.
% 
% ## MATLAB and Julia version compatibility
% 
% This package has been tested on a variety of MATLAB versions.
% However, for some versions of Julia and MATLAB, supported versions of external libraries may clash.
% For example, running JLCALL using Julia v1.6.1 and MATLAB R2015b gives the following error:
% 
%   >> JLCALL
%   
%   ERROR: Unable to load dependent library ~/.local/julia-1.6.1/bin/../lib/julia/libjulia-internal.so.1
%   
%   Message: /usr/local/MATLAB/R2015b/sys/os/glnxa64/libstdc++.so.6: version GLIBCXX_3.4.20' not found (required by ~/.local/julia-1.6.1/bin/../lib/julia/libjulia-internal.so.1)
% 
% This error results due to a clash of supported libstdc++ versions, and does not occur when using e.g. Julia v1.5.4 with MATLAB R2015b, or Julia v1.6.1 with MATLAB R2020b.
% 
% If you encounter this issue, see the <a href="matlab: web('https://github.com/JuliaLang/julia/blob/master/doc/build/build.md#required-build-tools-and-external-libraries')">Julia</a> and <a href="matlab: web('https://www.mathworks.com/support/requirements/supported-compilers.html')">MATLAB</a> documentation for information on mutually supported external libraries.
% 
% ## About this package
% 
% This repository contains utilities for parsing and running Julia code, passing MATLAB arguments to Julia, and retrieving Julia outputs from MATLAB.
% 
% The workhorse behind MATDaemon.jl and JLCALL is <a href="matlab: web('https://github.com/dmolina/DaemonMode.jl')">DaemonMode.jl</a> which is used to start a persistent Julia server in the background.

    % Parse inputs
    [f_args, opts] = parse_inputs(varargin{:});

    % Initialize workspace for communicating between MATLAB and Julia
    init_workspace(opts);

    % Optionally start persistent Julia server
    if opts.server
        manage_server(opts);
        if opts.shutdown
            return
        end
    end

    % Call Julia
    varargout = call_julia(f_args, opts);

end

function [f_args, opts] = parse_inputs(varargin)

    p = inputParser;

    addOptional(p, 'f', '(args...; kwargs...) -> nothing', @ischar);
    addOptional(p, 'args', {}, @iscell);
    addOptional(p, 'kwargs', struct, @isstruct);
    addParameter(p, 'infile', [tempname, '.mat'], @ischar);
    addParameter(p, 'outfile', [tempname, '.json'], @ischar);
    addParameter(p, 'runtime', try_find_julia_runtime, @ischar);
    addParameter(p, 'project', '', @ischar);
    addParameter(p, 'threads', 'auto', @validate_threads);
    addParameter(p, 'setup', '', @ischar);
    addParameter(p, 'nofun', false, @(x) validateattributes(x, {'logical'}, {'scalar'}));
    addParameter(p, 'modules', {}, @iscell);
    addParameter(p, 'cwd', pwd, @ischar);
    addParameter(p, 'workspace', relative_path('.jlcall'), @ischar);
    addParameter(p, 'server', true, @(x) validateattributes(x, {'logical'}, {'scalar'}));
    addParameter(p, 'port', 3000, @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'positive'}));
    addParameter(p, 'shared', true, @(x) validateattributes(x, {'logical'}, {'scalar'}));
    addParameter(p, 'restart', false, @(x) validateattributes(x, {'logical'}, {'scalar'}));
    addParameter(p, 'shutdown', false, @(x) validateattributes(x, {'logical'}, {'scalar'}));
    addParameter(p, 'gc', true, @(x) validateattributes(x, {'logical'}, {'scalar'}));
    addParameter(p, 'debug', false, @(x) validateattributes(x, {'logical'}, {'scalar'}));
    addParameter(p, 'source', 'MATDaemon', @ischar);

    parse(p, varargin{:});
    opts = p.Results;

    % Split Julia `f` args and kwargs out from options
    %   Note: must create `f_args` struct one field at a time as the `struct` constructor treats cell arrays specially;
    %   e.g. `s = struct('args', args)` would have `size(s) == size(args)` and `s(ii).args == args{ii}`, as opposed to
    %   the desired 1x1 struct `s` with `s.args == args`
    f_args = struct;
    f_args.args = opts.args;
    f_args.kwargs = opts.kwargs;
    opts = rmfield(opts, {'args', 'kwargs'});

    % For convenience, if `f` is empty replace it with a dummy function
    if isempty(opts.f)
        opts.f = '(args...; kwargs...) -> nothing';
    end

end

function init_workspace(opts)

    % Return if workspace is initialized
    if exist(opts.workspace, 'dir') && exist(fullfile(opts.workspace, 'Project.toml'), 'file')
        return
    end

    % Ignored outputs are needed to mute "folder exists" warning
    [~, ~] = mkdir(opts.workspace);

    % Install MATDaemon into workspace
    % install_script = build_julia_script(opts, 'Pkg', {
    %     'println("* Installing MATDaemon...\n")'
    %     sprintf('Pkg.add("MATDaemon"; io = %s)', jl_maybe_stdout(opts.debug))
    % });

    %TODO: Add option to provide MATDaemon path!!!!
    install_script = build_julia_script(opts, 'Pkg', {
        'println("* Installing MATDaemon...\n")'
        sprintf('Pkg.add(url=%s; io = %s)',opts.source, jl_maybe_stdout(opts.debug))
    });

    try_run(opts, install_script, 'client', 'Running `MATDaemon` install script');

end

function manage_server(opts)

    mlock % Prevent MATLAB from clearing persistent variables via e.g. `clear all`
    persistent cleanup_server % Julia server cleanup object

    if opts.restart || opts.shutdown
        cleanup_server = []; % triggers server cleanup, if server has been started
        if opts.shutdown
            return
        end
    end
    is_server_off = isempty(cleanup_server);

    if is_server_off
        % Initialize Julia server
        if opts.debug
            fprintf('* Starting Julia server\n\n');
        end

        % If shared is false, each Julia server call is executed in it's own Module to avoid namespace collisions, etc.
        start_script = build_julia_script(opts, 'MATDaemon', {
            sprintf('MATDaemon.start(%d; shared = %s, verbose = %s)', opts.port, jl_bool(opts.shared), jl_bool(opts.debug))
        });

        try_run(opts, start_script, 'server', 'Running `MATDaemon.start` script from Julia server');

        % Wait for server pong
        while ~ping_server(opts)
            pause(0.1);
        end

        % Kill server and collect garbage on MATLAB exit
        cleanup_server = onCleanup(@() kill_server(opts));
    end

end

function succ = ping_server(opts)

    try
        tcpclient('127.0.0.1', opts.port);
        succ = true;
    catch me
        if strcmpi(me.identifier, 'MATLAB:networklib:tcpclient:cannotCreateObject')
            succ = false;
        else
            rethrow(me)
        end
    end

end

function kill_server(opts)

    if opts.debug
        fprintf('* Killing Julia server\n\n');
    end

    kill_script = build_julia_script(opts, 'MATDaemon', {
        sprintf('MATDaemon.kill(%d; verbose = %s)', opts.port, jl_bool(opts.debug))
    });

    try_run(opts, kill_script, 'client', 'Sending kill script to Julia server');

    if opts.gc
        collect_garbage(opts);
    end

end

function output = call_julia(f_args, opts)

    % Save `f` arguments to `opts.infile`
    save(opts.infile, '-struct', 'f_args', '-v7.3');

    % Save input parser results to .mat file in workspace folder
    save(fullfile(opts.workspace, 'jlcall_opts.mat'), '-struct', 'opts', '-v7.3');

    % Script to run from Julia
    job_script = build_julia_script(opts, 'MATDaemon', {
        'include(MATDaemon.jlcall_script())'
    });

    if opts.server
        % Script to call the Julia server
        server_script = build_julia_script(opts, 'MATDaemon', {
            sprintf('MATDaemon.DaemonMode.runfile(raw"%s"; port = %d)', job_script, opts.port)
        });

        % Call out to Julia server
        try_run(opts, server_script, 'client', 'Sending `DaemonMode.runfile` script to Julia server');
    else
        % Call out to local Julia process
        try_run(opts, job_script, 'local', 'Calling `MATDaemon.jlcall` from local Julia process');
    end

    % Load outputs from disk
    if exist(opts.outfile, 'file')
        %output = load(opts.outfile);
        %output = output.output;
        fid = fopen(opts.outfile); 
        raw = fread(fid,inf); 
        str = char(raw'); 
        fclose(fid); 
        output = {jsondecode(str)};
    else
        % Throw error before garbage collecting below so that workspace folder can be inspected
        e.message = sprintf('Julia call failed to produce the expected output file:\n%s', opts.outfile);
        e.identifier = 'jlcall:fileNotFound';
        error(e)
    end

    % Collect temporary garbage
    if opts.gc
        collect_garbage(opts);
    end

end

function jl_script = build_julia_script(opts, pkgs, body)

    if nargin < 2; body = {}; end
    if nargin < 1; pkgs = {}; end

    if ischar(pkgs); pkgs = {pkgs}; end
    if ischar(body); body = {body}; end

    % Create temporary file for Julia script
    jl_script = [workspace_tempname(opts), '.jl'];
    fid = fopen(jl_script, 'w');
    cleanup_fid = onCleanup(@() fclose(fid));

    % MATDaemon workspace should always be on the Julia LOAD_PATH
    preamble = {
        'if !in(ENV["MATDAEMON_WORKSPACE"], LOAD_PATH)'
        '    pushfirst!(LOAD_PATH, ENV["MATDAEMON_WORKSPACE"])'
        'end'
    };

    % Build script
    for ii = 1:length(preamble)
        fprintf(fid, '%s\n', preamble{ii});
    end
    for ii = 1:length(pkgs)
        fprintf(fid, 'import %s\n', pkgs{ii});
    end
    for ii = 1:length(body)
        fprintf(fid, '%s\n', body{ii});
    end

end

function try_run(opts, script, mode, msg)

    if nargin < 4
        msg = 'Command';
    end

    % Set MATDaemon environment variables
    setenv('MATDAEMON_WORKSPACE', opts.workspace);

    % Set Julia binary path and flags
    flags = sprintf('--project=%s --threads=%s --startup-file=no', opts.workspace, jl_threads(opts.threads));
    switch mode
        case 'server'
            flags = [flags, ' --optimize=3'];
            detach = ' &';
        case 'client'
            flags = [flags, ' --optimize=0 --compile=min'];
            detach = '';
        case 'local'
            flags = [flags, ' --optimize=3'];
            detach = '';
        otherwise
            error('Unknown mode: ''%s''', mode)
    end

    % Build and run Julia command
    cmd = [opts.runtime, ' ', flags, ' ', script, detach];
    st = system(cmd);

    if opts.debug
        fprintf('* %s (status = %d):\n*   %s\n\n', msg, st, cmd);
    end

end

function runtime = try_find_julia_runtime()

    % Default value
    runtime = 'julia';

    try
        if isunix
            [st, res] = system('which julia');
        elseif ispc
            [st, res] = system('where julia');
        else
            return % default to 'julia'
        end
        if st == 0
            runtime = strtrim(res);
        end
    catch me
        % ignore error; default to 'julia'
    end

end

function tmp = workspace_tempname(opts)

    tmp_dir = fullfile(opts.workspace, 'tmp');
    [~, ~] = mkdir(tmp_dir); % ignore "folder exists" warning
    [dirname, filename] = fileparts(tempname(tmp_dir));

    persistent filecount
    if isempty(filecount)
        filecount = 0;
    else
        filecount = mod(filecount + 1, 10000);
    end

    tmp = fullfile(dirname, [pad_num(filecount), '_mat_', filename]);

end

function collect_garbage(opts)

    if exist(opts.infile, 'file'); delete(opts.infile); end
    if exist(opts.outfile, 'file'); delete(opts.outfile); end
    delete(fullfile(opts.workspace, 'tmp', '*'));
    delete(fullfile(opts.workspace, '*.mat'));

end

function path = relative_path(varargin)

    jlcall_dir = fileparts(mfilename('fullpath'));
    path = fullfile(jlcall_dir, varargin{:});

end

function str = jl_bool(bool)

    if bool
        str = 'true';
    else
        str = 'false';
    end

end

function str = jl_threads(threads)

    if strcmpi(threads, 'auto')
        str = 'auto';
    else
        str = num2str(threads);
    end

end

function validate_threads(x)

    if ischar(x)
        if ~strcmpi(x, 'auto')
            error('Received value: ''%s''. Expected either ''auto'' or a positive integer', x)
        end
    else
        validateattributes(x, {'numeric'}, {'scalar', 'integer', 'positive'})
    end

end

function str = jl_maybe_stdout(bool)

    if bool
        str = 'stdout';
    else
        str = 'devnull';
    end

end

function str = pad_num(s)

    if ~ischar(s)
        s = num2str(s);
    end

    str = [repmat('0', 1, max(4 - numel(s), 0)), s];

end
