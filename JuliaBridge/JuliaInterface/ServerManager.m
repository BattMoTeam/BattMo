classdef ServerManager < handle 

    properties
        
        options
        base_call
        use_folder
        call_history = {}
        
    end
    
    properties (Hidden)
        
        cleanup
        
    end
    
    methods
        
        function manager = ServerManager(varargin)
        % Parse inputs

            p=inputParser;

            serverFolder = fileparts(mfilename('fullpath'));
            addParameter(p, 'julia', try_find_julia_runtime, @ischar);
            addParameter(p, 'project', fullfile(serverFolder, 'RunFromMatlab'), @ischar);
            addParameter(p, 'script_source', fullfile(serverFolder, 'DaemonHandler.jl'), @ischar);
            addParameter(p, 'startup_file', 'no', @ischar);
            addParameter(p, 'threads', 'auto', @validate_threads);
            addParameter(p, 'cwd', serverFolder, @ischar);
            addParameter(p, 'port', 3000, @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'positive'}));
            addParameter(p, 'shared', true, @(x) validateattributes(x, {'logical'}, {'scalar'}));
            addParameter(p, 'print_stack', true, @(x) validateattributes(x, {'logical'}, {'scalar'}));
            addParameter(p, 'async', true, @(x) validateattributes(x, {'logical'}, {'scalar'}));
            addParameter(p, 'gc', true, @(x) validateattributes(x, {'logical'}, {'scalar'}));
            addParameter(p, 'debug', false, @(x) validateattributes(x, {'logical'}, {'scalar'}));

            parse(p, varargin{:});
            manager.options=p.Results;
            
            manager.use_folder = [manager.options.cwd, '/Server_id=',num2str(manager.options.port)];
            [~,~] = mkdir([manager.use_folder,'/tmp']);

            % Save options in file
            op       = manager.options;
            opt_file = [manager.use_folder,'/options.mat'];
            save(opt_file, "op");
       
            % Build DaemonMode call:
            % Manage options
            manager.base_call = [manager.options.julia, ' '                              , ...
                                 '--startup-file='    , manager.options.startup_file, ' ', ...
                                 '--project='         , manager.options.project     , ' ', ...
                                 '--threads='         , manager.options.threads];
            
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
        
        function load(manager, data, kwargs)
            % Load data onto server (requires shared=true to make sense)
            % Add warning if shared is false

            write_to_file = [fullfile(manager.use_folder,tempname), '.mat'];
            call_load = ['"using Revise, DaemonMode; runargs(', num2str(manager.options.port), ')" ', manager.options.script_source, ...
                ' -load ', write_to_file];
            save(write_to_file, "data", "kwargs");
            if manager.options.debug
                fprintf("Loading data into Julia \n")
            end
            manager.DaemonCall(call_load);
            
        end

        function result = run_battery(manager, use_data, data, kwargs)
            
            % Run battery, either using already loaded data or after loading
            if ~use_data
                manager.load(data, kwargs);
            end
            
            load_from_file = [fullfile(manager.use_folder, tempname),'.json'];
            call_battery = ['"using Revise, DaemonMode; runargs(',num2str(manager.options.port) ,')" ', manager.options.script_source, ...
                ' -run_battery ', load_from_file];
            
            if manager.options.debug
                fprintf("Calling run battery \n")
            end
            
            st = manager.DaemonCall(call_battery);
            
            if st
                % Read julia output
                fid = fopen(load_from_file); 
                raw = fread(fid, inf); 
                str = char(raw'); 
                fclose(fid); 
                result = {jsondecode(str)};

                if manager.options.gc
                    delete(load_from_file);
                end
            else
                result=[];
            end
        end
        
        function iterate_values(manager, values)
            %Run parameter values as flags to julia script if shared is
            %true. Possible enable same procedure entirely done in Julia
        end

        function startup(manager)
        % Create serve call
            
            if ~ping_server(manager.options)
                startup_call = ['"using Revise, DaemonMode; serve(', ...
                                num2str(manager.options.port), ', ', jl_bool(manager.options.shared), ...
                                ', print_stack=', jl_bool(manager.options.print_stack), ...
                                ', async=', jl_bool(manager.options.async), ')" &'];
                
                if manager.options.debug
                    fprintf("Starting Julia server \n")
                end

                manager.DaemonCall(startup_call);
                
                while ~ping_server(manager.options)
                    pause(0.1);
                end
            end
            
        end

        function delete(manager)
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
    end
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

