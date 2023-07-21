clear jlcall_mod

casenames = {'p1d_40',
             'p2d_40'};

casenames = casenames(2);
% casenames = {'3d_demo_case'};

% casenames = {'4680_case'};

battmo_folder = fileparts(mfilename('fullpath'));
battmo_folder = fullfile(battmo_folder, '../../BattMo.jl');

jsonfolder = fullfile(battmo_folder, 'test/battery/data/jsonfiles/');

for icase = 1 : numel(casenames)
    casename = casenames{icase};
    [model,schedule,state0]=setupMatlabModel(casename, jsonfolder);
end

jlcall('', ...
     'project', '/data2/andreas/json_experiment/json_env', ... % activate a local Julia Project
     'setup', '/data2/andreas/json_experiment/BattMo/JuliaBridge/setup.jl',...
     'modules', {'json_env'}, ... % load a custom module and some modules from Base Julia
     'threads', 'auto', ... % use the default number of Julia threads
     'restart', true, ... % start a fresh Julia server environment
     'debug',true ...
     )

export=struct('model', model, ...
                'schedule',schedule,...
                'state0', state0,...
                'states',[]);
jlcall('json_env.setup_wrapper', {export})



