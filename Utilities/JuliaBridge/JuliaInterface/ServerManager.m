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

        function manager = ServerManager(varargin)
        % Parse inputs

            serverFolder = fileparts(mfilename('fullpath'));

            p = inputParser;
            addParameter(p , 'julia'              , try_find_julia_runtime(), @ischar);
            addParameter(p , 'project'            , fullfile(serverFolder, 'RunFromMatlab')         , @ischar);
            addParameter(p , 'script_source'      , fullfile(serverFolder, 'RunFromMatlab','api','DaemonHandler.jl'), @ischar);
            addParameter(p , 'startup_file'       , 'no'        , @ischar);
            addParameter(p , 'threads'            , 'auto'      , @isnumeric);
            %addParameter(p, 'procs'              , 8           , @isnumeric);
            addParameter(p , 'cwd'                , serverFolder, @ischar);
            addParameter(p , 'port'               , 3000        , @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'positive'}));
            addParameter(p , 'shared'             , true        , @(x) validateattributes(x, {'logical'}, {'scalar'}));
            addParameter(p , 'verbose'            , true        , @(x) validateattributes(x, {'logical'}, {'scalar'}));
            addParameter(p , 'print_stack'        , true        , @(x) validateattributes(x, {'logical'}, {'scalar'}));
            addParameter(p , 'async'              , true        , @(x) validateattributes(x, {'logical'}, {'scalar'}));
            addParameter(p , 'gc'                 , true        , @(x) validateattributes(x, {'logical'}, {'scalar'}));
            addParameter(p , 'debug'              , false       , @(x) validateattributes(x, {'logical'}, {'scalar'}));
            addParameter(p , 'reset'              , false       , @(x) validateattributes(x, {'logical'}, {'scalar'}));
            addParameter(p , 'updateJuliaPackages', false       , @(x) validateattributes(x, {'logical'}, {'scalar'}));
            addParameter(p , 'extra_args'         , ''          , @ischar);

            parse(p, varargin{:});
            manager.options = p.Results;

            % We will save all files related to the execution of the server
            % in a folder on the form Server_id=<port>
            manager.use_folder = fullfile(manager.options.cwd, ['Server_id=', num2str(manager.options.port)]);

            if manager.options.reset
                manager.deleteTempFolder();
            end

            [st, result] = mkdir(manager.use_folder);
            assert(st == 1, 'Unable to use mkdir: %s', result);

            manager.options.extra_args = [' ', manager.options.extra_args, ' '];

            % Save options in file to be read in Julia
            op = manager.options;
            opt_file = fullfile(manager.use_folder, 'options.mat');
            save(opt_file, 'op', manager.file_version);

            % Build DaemonMode call:

            % Manage options
            manager.base_call = [manager.options.julia, ' '                              , ...
                                 '--startup-file='    , manager.options.startup_file, ' ', ...
                                 '--project='         , manager.options.project     , ' ', ...
                                 '--threads='         , num2str(manager.options.threads), ...
                                 manager.options.extra_args];

            if isunix
                
                manager.base_call = ['env -u LD_LIBRARY_PATH ', manager.base_call];
                
            end

            if ~isfile(fullfile(serverFolder, 'RunFromMatlab', 'Manifest.toml'))

                fprintf('Setting up Battmo Julia server (may take 1 minute) ... ');
                setup_battmojl_call = [''''                                     , ...
                                       'using Pkg;'                             , ...
                                       'Pkg.add(name = "BattMo", rev = "main");', ...
                                       'Pkg.update();'                          , ...
                                       ''''];
                manager.DaemonCall(setup_battmojl_call);
                fprintf('done\n');
                
            end
            
            if manager.options.updateJuliaPackages

                manager.updateJuliaPackages()
                
            end

            % Start server
            manager.startup();

            % Read options in Julia
            option_call = ['"using Revise, DaemonMode; runargs(', ...
                           num2str(manager.options.port)        , ...
                           ')" '                                , ...
                           manager.options.script_source        , ...
                           ' -load_options '                    , ...
                           opt_file];

            manager.DaemonCall(option_call);

            % Shut down server on cleanup
            manager.cleanup = onCleanup(@() delete(manager));

        end

        function updateJuliaPackages(manager)

            fprintf('Updating packages on Julia server (may take few seconds) ... ');
            instantiate_call = '"using Pkg; Pkg.update();"';
            manager.DaemonCall(instantiate_call);
            fprintf('done\n');
            
        end

        function load(manager, varargin)
        % Load data onto server (requires shared = true to make sense)
        % Add warning if shared is false
        %
        %  In varargin:
        %
        % 'inputType' : String, either 'Matlab' (default) or 'JSON'
        % 'kwargs'    : Struct, which is passed to the solver
        %
        % If 'inputType' is 'Matlab'
        %   - 'data' : Structure with fields 'model', 'schedule', 'initState'
        %
        % If 'inputType' is 'JSON'
        %   - 'inputFileName' : String, path to json file
        
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

            call_load = ['"using Revise, DaemonMode; runargs(', ...
                         num2str(manager.options.port)        , ...
                         ')" '                                , ...
                         manager.options.script_source        , ...
                         ' -load '                            , ...
                         loadingDataFilename];

            if manager.options.debug | manager.options.verbose
                fprintf("Loading data into Julia \n")
            end

            manager.DaemonCall(call_load);

        end

        function result = run(manager)

            outputFileName = [tempname,'.json'];

            %Call DaemonMode.runargs
            call_battery = ['"using DaemonMode; runargs(', ...
                            num2str(manager.options.port), ...
                            ')" '                        , ...
                            manager.options.script_source, ...
                            ' -run '                     , ...
                            outputFileName];

            if manager.options.debug | manager.options.verbose
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
                result = battMojsondecode(str);

                if manager.options.gc
                    delete(outputFileName);
                end

            else
                result = [];
            end

        end

        function result = runBatteryJson(manager, jsonstruct, varargin)

            opt = struct('useDirectJsonInput', false);
            opt = merge_options(opt, varargin{:});

            if opt.useDirectJsonInput
                
                assert(strcmp(jsonstruct.Geometry.case, '1D'), 'Direct Json input for battmo.jl is only available for 1D geometry');
                inputType = 'JSON';

                jsonstruct = jsonencode(jsonstruct);
                inputFileName  = [tempname, '.json'];
                fileID = fopen(inputFileName,'w');
                fprintf(fileID,'%s', jsonstruct);
                fclose(fileID);

                data = [];
                
            else
                
                inputType     = 'Matlab';
                data          = setupSimulationForJuliaBridge(jsonstruct);
                inputFileName = [];
                
            end

            manager.load('inputType'    , inputType, ...
                         'data'         , data     , ...
                         'inputFileName', inputFileName);

            result = manager.run();
            
        end

        function startup(manager)

            if ~ping_server(manager.options) %If server is already live, do nothing
                                             % Create DaemonMode.serve call

                if ispc

                    message = ['You are running on a Windows machine.'                                                           , ...
                               ' In this case Matlab cannot start a process in the background and you will have to '             , ...
                               'proceed manually as follows\n\n'                                                                 , ...
                               'First, julia should be installed, see https://julialang.org . '                                  , ...
                               'Then, open a terminal (NOT a powershell) and run:\n\n'                                           , ...
                               'julia --startup-file=no --project=/path/to/RunFromMatlab'                                        , ...
                               ' -e "using Revise, DaemonMode; serve(3000, true, call_stack=true, async=true)"\n\n'              , ...
                               'where /path/to/RunFromMatlab.m is the path to the directory '                                    , ...
                               '/Utilities/JuliaBridge/JuliaInterface/RunFromMatlab in you BattMo installation folder\n\n'       , ...
                               'Running this command will block the command prompt. '                                            , ...
                               'The server will remain active until the window is closed or it is deactivated in any other way. ', ...
                               'Calls to the server can now be made using the ServerManager\n\n'                                 , ...
                              ];

                    fprint(message);

                else

                    startup_call = ['"using Revise, DaemonMode; serve(' , ...
                                    num2str(manager.options.port)       , ...
                                    ', '                                , ...
                                    jl_bool(manager.options.shared)     , ...
                                    ', print_stack='                    , ...
                                    jl_bool(manager.options.print_stack), ...
                                    ')" &'];
                    % ', async='                          , ...
                    % jl_bool(manager.options.async)      , ...

                    if manager.options.debug | manager.options.verbose
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

        end

        function shutdown(manager)
        % Close server if active

            kill_call= ['"using Revise, DaemonMode; sendExitCode(', ...
                        num2str(manager.options.port), ...
                        ');"'];
            cmd = [manager.base_call, ' -e ', kill_call];
            [st, result] = system(cmd);
            assert(st == 0, "System call failed: \n %s \nSystem call returned:\n %s\n", cmd, result);

            if manager.options.debug | manager.options.verbose
                fprintf("Shutting down server \n");
            end

        end

        function restart(manager, varargin)
        % Close server and restart

            manager.shutdown();
            manager.startup()

        end

        function success = DaemonCall(manager, call)
        % Call the DaemonMode server

            cmd = [manager.base_call, ' -e ', call];

            if manager.options.debug
                fprintf("Call to julia: %s \n", cmd);
            elseif manager.options.verbose
                fprintf("Call to julia (It will take time the first time due to just in time compilation)\n");
                fprintf("  %s \n", cmd);
            end

            [st, result] = system(cmd);
            assert(st == 0, "System call failed: \n %s \nSystem call returned:\n %s\n", cmd, result);

            success = true;

            manager.call_history{end+1} = cmd;

        end

        function f = sweep(manager, experiment, values, name)
        % Perform a parameter sweep

            save_folder = fullfile(manager.use_folder, name);
            [st, result] = mkdir(save_folder);
            assert(st == 1, 'Unable to use mkdir: %s', result);

            options_file = fullfile(save_folder, [name, '.mat']);
            save(options_file, 'save_folder', 'experiment', 'values', manager.file_version);

            sweep_source = fullfile(fileparts(mfilename('fullpath')), 'RunFromMatlab', 'api', 'ParameterSweepControl.jl');
            sweep_call = ['"using Revise, DaemonMode; runargs(', num2str(manager.options.port), ')" ', sweep_source, ' ', options_file, ' &'];
            manager.DaemonCall(sweep_call);

            f = checkSweep(experiment, save_folder);

        end

        function result = collect_results(manager, f, index)
        % Collect results from a Futures object at given indices

            if manager.options.debug | manager.options.verbose
                disp("Waiting for futures object to finish...")
            end

            while isprop(f, 'State') && f.State == "running"
                pause(0.1);
            end

            if manager.options.debug | manager.options.versbose
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
                states = battMojsondecode(str);
                result(i).states = states;
            end
        end


        function deleteTempFolder(manager)

            if mrstPlatform('octave')
                confirm_recursive_rmdir(0);
            end
            [st, result] = rmdir(manager.use_folder, 's');
            assert(st == 1, 'Deleting temp folder failed with reason: %s', result);

        end

    end
end


function [result, count] = checkSweep(experiment, save_folder)

    output_file = fullfile(save_folder, [experiment, '_output.json']);

    while ~exist(output_file,'file')
        count = length(dir(fullfile(save_folder,'*.json')));
        pause(0.1);
    end

    % Read file
    fid = fopen(output_file);
    raw = fread(fid, inf);
    str = char(raw');
    fclose(fid);
    result = battMojsondecode(str);

end


function runtime = try_find_julia_runtime()
% Determine correct julia source

    % Default value
    runtime = 'julia';

    try
        if isunix
            [st, res] = system('which julia');
            if st == 1
                % Try juliaup directory
                home = getenv('HOME');
                runtime = [strtrim(home), '/.juliaup/bin/julia']
            end
        elseif ispc
            [st, res] = system('where julia');
        else
            return % default to 'julia'
        end
        if st == 0
            runtime = strtrim(res);
        end
    catch
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



%{
Copyright 2021-2024 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BattMo

BattMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BattMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}
