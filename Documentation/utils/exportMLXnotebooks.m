function exportMLXnotebooks()
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

    inputdir  = fullfile(battmoDir(), 'Examples', 'Notebooks');
    outputdir = fullfile(battmoDir(), 'Documentation', '_static', 'notebooks');

    %% Helper functions
    % exportMLXtoMfile
    % openMfileSaveToMlx
    % runMlxAndSave
    % exportMLXtoHtml
    
    function exportMLXtoMfile(notebookname)

        options = {'format', 'm', ...
                   'run', false};
        
        inputfile = fullfile(inputdir, notebookname);
        
        notebookname = strrep(notebookname, '.mlx', '.m');
        outputfile = fullfile(outputdir, notebookname);
        
        export(inputfile, outputfile, options{:});
        
    end

    function openMfileSaveToMlx(notebookname)

        notebookname = sprintf('%s.m', notebookname);
        inputfile = fullfile(inputdir, notebookname);
        
        notebookname = strrep(notebookname, '.m', '.mlx');
        outputfile = fullfile(outputdir, notebookname);
        
        matlab.internal.liveeditor.openAsLiveCode(fileread(inputfile)); % open as live script
        activeDoc = matlab.desktop.editor.getActive(); % get the active editor (the new file)
        activeDoc.saveAs(outputfile); % Save with same file name but .mlx
        activeDoc.close()
        
    end

    function runMlxAndSave(notebookname)

    % To run and update the mlx notebook programmatically, it is possible to use: matlab.internal.liveeditor.executeAndSave('fullpathnameto.mlx')
        
        inputfile = fullfile(inputdir, notebookname);
        matlab.internal.liveeditor.executeAndSave(inputfile);

    end

    function exportMLXtoHtml(notebookname)

        options = {'format', 'html', ...
                   'run', false};
        
        inputfile = fullfile(inputdir, notebookname);
        
        notebookname = strrep(notebookname, '.mlx', '.html');
        outputfile = fullfile(outputdir, notebookname);
        
        export(inputfile, outputfile, options{:});
        
    end

    %% Main loop
    
    for inote = 1 : numel(notebooknames)

        notebookname = notebooknames{inote};

        fprintf('Starting with %s ...', notebookname);

        exportMLXtoHtml(notebookname);

        fprintf('Done.\n', notebookname);
        
    end


end
