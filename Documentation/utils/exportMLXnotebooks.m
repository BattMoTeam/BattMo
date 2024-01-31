%% Export to Python notebooks require MATLAB versions >= 2023

notebooknames = {'tutorial_1_a_simple_p2d_model_live.mlx'              , ...
                 'tutorial_2_changing_control_protocol_live.mlx'       , ...
                 'tutorial_3_modify_structural_parameters_live.mlx'    , ...
                 'tutorial_4_modify_material_parameters_live.mlx'      , ...
                 'tutorial_5_simulate_CCCV_cycling_live.mlx'           , ...
                 'tutorial_6_simulate_thermal_performance_live.mlx'    , ...
                 'tutorial_7_a_simple_p4d_model_live.mlx'              , ...
                 'tutorial_8_simulate_a_multilayer_pouch_cell_live.mlx', ...
                 'tutorial_9_simulate_a_cylindrical_cell_live.mlx'};

inputdir  = fullfile(battmoDir(), 'Examples', 'NoteBooks');
outputdir = fullfile(battmoDir(), 'Documentation', '_static', 'notebooks');

%% setup options
% Most of the options are in fact the default options.
options = {'format', 'html', ...
           'run', false};

% To run and update the mlx notebook programmatically, it is possible to use: matlab.internal.liveeditor.executeAndSave('fullpathnameto.mlx')

for inote = 1 : numel(notebooknames)

    notebookname = notebooknames{inote};
    inputfile = fullfile(inputdir, notebookname);
    
    notebookname = strrep(notebookname, '.mlx', '.html');
    outputfile = fullfile(outputdir, notebookname);
    
    export(inputfile, outputfile, options{:});
    fprintf('Exported %s\n', notebookname);
    
end
