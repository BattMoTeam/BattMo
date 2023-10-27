classdef ServerManager < handle

    properties

        options
        base_call
        use_folder
        call_history = {}
        file_version = '-v7'

    end

    properties (Hidden)

        cleanup

    end

    methods

        function manager    = ServerManager(varargin)
        % Parse inputs

            serverFolder = fileparts(mfilename('fullpath'));

            p = inputParser;
            addParameter(p, 'julia'        , try_find_julia_runtime, @ischar);
            addParameter(p, 'project'      , fullfile(serverFolder , 'RunFromMatlab') , @ischar);
            addParameter(p, 'script_source', fullfile(serverFolder , 'RunFromMatlab','api','DaemonHandler.jl'), @ischar);
            addParameter(p, 'startup_file' , 'no'                  , @ischar);
            addParameter(p, 'threads'      , 'auto'                , @validate_threads);
            addParameter(p, 'procs'        , 1                     , @(x) validateattributes(x, {'numeric'}, {'integer'}));
            addParameter(p, 'cwd'          , serverFolder          , @ischar);
            addParameter(p, 'port'         , 3000                  , @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'positive'}));
            addParameter(p, 'shared'       , true                  , @(x) validateattributes(x, {'logical'}, {'scalar'}));
            addParameter(p, 'print_stack'  , true                  , @(x) validateattributes(x, {'logical'}, {'scalar'}));
            addParameter(p, 'async'        , true                  , @(x) validateattributes(x, {'logical'}, {'scalar'}));
            addParameter(p, 'gc'           , true                  , @(x) validateattributes(x, {'logical'}, {'scalar'}));
            addParameter(p, 'debug'        , true                 , @(x) validateattributes(x, {'logical'}, {'scalar'}));

            parse(p, varargin{:});
            manager.options=p.Results;

            % We will save all files related to the execution of the server
            % in a folder on the form Server_id=<port>
            manager.use_folder = fullfile(manager.options.cwd, ['Server_id=', num2str(manager.options.port)]);
            [~, ~] = mkdir(manager.use_folder);

            %Save options in file to be read in Julia
            op = manager.options;
            opt_file = fullfile(manager.use_folder, 'options.mat');
            save(opt_file, 'op', manager.file_version);

            % Build DaemonMode call:

            % Manage options
            manager.base_call = [manager.options.julia, ' '                              , ...
                                 '--startup-file='    , manager.options.startup_file, ' ', ...
                                 '--project='         , manager.options.project     , ' ', ...
                                 '--procs='           , num2str(manager.options.procs), ' ', ...
                                 '--threads='         , manager.options.threads];

            %Ensure that project is instantiated
            instaniate_call = '"using Pkg; Pkg.instantiate();"';
            manager.DaemonCall(instaniate_call);

            % Start server
            manager.startup();

            % Read options in Julia
            option_call = ['"using Revise, DaemonMode; runargs(', ...
                           num2str(manager.options.port),')" ', manager.options.script_source, ...
                           ' -load_options ', opt_file];

            manager.DaemonCall(option_call);

            % Shut down server on cleanup
            manager.cleanup = onCleanup(@() delete(manager));

        end

        function load(manager, varargin)
        % Load data onto server (requires shared=true to make sense)
        % Add warning if shared is false

            opts = struct('data'         , []      , ...
                          'kwargs'       , []      , ...
                          'inputType'    , 'Matlab', ...
                          'use_state_ref', false   , ...
                          'inputFileName', []);
            opts = merge_options(opts, varargin{:});

            data          = opts.data;
            kwargs        = opts.kwargs;
            inputType     = opts.inputType;
            use_state_ref = opts.use_state_ref;

            loadingDataFilename = [tempname, '.mat'];

            switch opts.inputType
              case 'Matlab'
                inputFileName = loadingDataFilename;
              case 'JSON'
                inputFileName = opts.inputFileName;
              otherwise
                error('inputType not recognized');
            end

            save(loadingDataFilename, ...
                 'data'             , ...
                 'kwargs'           , ...
                 'inputType'        , ...
                 'use_state_ref'    , ...
                 'inputFileName'    , ...
                 manager.file_version);

            call_load = ['"using Revise, DaemonMode; runargs(', num2str(manager.options.port), ')" ', manager.options.script_source, ...
                         ' -load ', loadingDataFilename];
            if manager.options.debug
                fprintf("Loading data into Julia \n")
            end
            manager.DaemonCall(call_load);

        end

        function result = run_battery(manager)

            outputFileName = [tempname,'.json'];

            %Call DaemonMode.runargs
            call_battery = ['"using Revise, DaemonMode; runargs(',num2str(manager.options.port) ,')" ', manager.options.script_source, ...
                            ' -run_battery ', outputFileName];

            if manager.options.debug
                fprintf("Calling run battery \n")
            end

            st = manager.DaemonCall(call_battery);
            %Read only if system call completed succesfully
            if st
                % Read julia output
                fid = fopen(outputFileName);
                raw = fread(fid, inf);
                str = char(raw');
                fclose(fid);
                result = {jsondecode(str)};

                if manager.options.gc
                    delete(outputFileName);
                end
            else
                result=[];
            end

        end

        function startup(manager)

            if ~ping_server(manager.options) %If server is already live, do nothing
                                             % Create DaemonMode.serve call
                startup_call = ['"using Revise, DaemonMode; serve(', ...
                                num2str(manager.options.port), ', ', jl_bool(manager.options.shared), ...
                                ', print_stack=', jl_bool(manager.options.print_stack),')" &'];
                %, ...
                %               ', async=', jl_bool(manager.options.async),')" &'];

                if manager.options.debug
                    fprintf("Starting Julia server \n")
                end

                manager.DaemonCall(startup_call);

                % Check if server is active. Ensures that we do not make
                % calls to the server until it is ready
                while ~ping_server(manager.options)
                    pause(0.1);
                end
            end

        end

        function shutdown(manager)
        % Close server if active
            kill_call= ['"using Revise, DaemonMode; sendExitCode(', num2str(manager.options.port), ...
                        ');"'];
            cmd = [manager.base_call, ' -e ', kill_call];
            system(cmd)
            if manager.options.debug
                fprintf("Shutting down server \n");
            end

        end

        function restart(manager, varargin)
        % Close server and restart

            kill_call= ['"using Revise, DaemonMode; sendExitCode(', num2str(manager.options.port), ...
                        ');"'];
            cmd = [manager.base_call, ' -e ', kill_call];
            system(cmd)
            if manager.options.debug
                fprintf("Shutting down server \n");
            end
            manager.startup()

        end

        function success = DaemonCall(manager, call)
        %Call the DaemonMode server

            cmd = [manager.base_call, ' -e ', call];

            if manager.options.debug
                fprintf("Call to julia: %s \n", cmd);
            end

            try
                st=system(cmd);
                success=true;
            catch
                fprintf("System call failed: \n %s", cmd);
                success=false;
            end

            manager.call_history{end+1} = cmd;

        end

        function f = sweep(manager, experiment, values, name)
        % Perform a parameter sweep
            save_folder = fullfile(manager.use_folder, name);
            [~,~] = mkdir(save_folder);
            options_file = [save_folder, '/', name, '.mat'];
            save(options_file, 'save_folder', 'experiment', 'values', manager.file_version);
            sweep_source = fullfile(fileparts(mfilename('fullpath')), 'RunFromMatlab','api','ParameterSweepControl.jl');
            sweep_call = ['"using Revise, DaemonMode; runargs(', num2str(manager.options.port), ')" ', sweep_source, ' ', options_file, ' &'];
            manager.DaemonCall(sweep_call);

            if mrstPlatform('octave')
                f = checkSweep(experiment, save_folder);
            else
                f = parfeval(@checkSweep, 1, experiment, save_folder);
            end

        end

        function result = collect_results(manager, f, index)
        %Collect results from a Futures object at given indices

            while isprop(f, 'State') && f.State == "running"
                pause(0.1);
                if manager.options.debug
                    disp("Waiting for futures object to finish")
                end
            end
            if manager.options.debug
                disp("Loading data from object")
            end

            try
                result = f.OutputArguments{1}.Results(index);
            catch
                result = f.Results(index);
            end

            for i = 1:length(index)
                fid = fopen(result(i).states_location);
                raw = fread(fid, inf);
                str = char(raw');
                fclose(fid);
                states = jsondecode(str);
                result(i).states = states;
            end
        end
    end
end

function [result, count] = checkSweep(experiment, save_folder)
    output_file = fullfile(save_folder, [experiment,'_output.json']);
    while ~exist(output_file,'file')
        count = length(dir(fullfile(save_folder,'*.json')));
        pause(0.1);
    end

    %Read file
    fid = fopen(output_file);
    raw = fread(fid, inf);
    str = char(raw');
    fclose(fid);
    result = jsondecode(str);
end

% Determine correct julia source
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

function str = jl_bool(bool)

    if bool
        str = 'true';
    else
        str = 'false';
    end

end

function succ = ping_server(opts)

    try
        obj = tcpclient('127.0.0.1', opts.port); %#ok Octave requires return object
        succ = true;
    catch me
        if (mrstPlatform('octave') && strcmpi(me.message, 'tcpclient: error on connect : 111 - Connection refused')) || ...
                (mrstPlatform('matlab') && strcmpi(me.identifier, 'MATLAB:networklib:tcpclient:cannotCreateObject'))
            succ = false;
        else
            rethrow(me);
        end
    end

end
