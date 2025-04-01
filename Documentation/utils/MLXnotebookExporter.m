classdef MLXnotebookExporter

    properties

        % list of registered notebooks
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

    end

    methods
        
        function exportMLXtoMfile(mne, notebookname, varargin)
            
            opt = struct('format', 'm', ...
                          'run', false);
            opt = merge_options(opt, varargin{:});
            
            inputfile = fullfile(mne.inputdir, notebookname);
            
            notebookname = strrep(notebookname, '.mlx', '.m');
            outputfile = fullfile(mne.outputdir, notebookname);
            
            optvals = cellfun(@(x) opt.(x), fieldnames(opt), 'uniformoutput', false);
            options = reshape(vertcat(fieldnames(opt)', optvals'), [], 1);
            
            export(inputfile, outputfile, options{:});
            
        end

        function openMfileSaveToMlx(mne, notebookname)

            notebookname = sprintf('%s.m', notebookname);
            inputfile = fullfile(mne.inputdir, notebookname);
            
            notebookname = strrep(notebookname, '.m', '.mlx');
            outputfile = fullfile(mne.outputdir, notebookname);
            
            matlab.internal.liveeditor.openAsLiveCode(fileread(inputfile)); % open as live script
            activeDoc = matlab.desktop.editor.getActive(); % get the active editor (the new file)
            activeDoc.saveAs(outputfile); % Save with same file name but .mlx
            activeDoc.close()
            
        end

        function runMlxAndSave(mne, notebookname)

        % To run and update the mlx notebook programmatically, it is possible to use: matlab.internal.liveeditor.executeAndSave('fullpathnameto.mlx')
            
            inputfile = fullfile(mne.inputdir, notebookname);
            matlab.internal.liveeditor.executeAndSave(inputfile);

        end

        function exportMLXtoHtml(mne, notebookname, varargin)

            opt = struct('format', 'html', ...
                          'run', false);
            opt = merge_options(opt, varargin{:});
            
            inputfile = fullfile(mne.inputdir, notebookname);
            
            notebookname = strrep(notebookname, '.mlx', '.html');
            outputfile = fullfile(mne.outputdir, notebookname);
            
            optvals = cellfun(@(x) opt.(x), fieldnames(opt), 'uniformoutput', false);
            options = reshape(vertcat(fieldnames(opt)', optvals'), [], 1);

            export(inputfile, outputfile, options{:});
            
        end


        function exportAll(mne)
            
            for inote = 1 : numel(mne.notebooknames)

                notebookname = mne.notebooknames{inote};

                fprintf('Starting with %s ...', notebookname);

                mne.exportMLXtoHtml(notebookname);

                fprintf('Done.\n', notebookname);
                
            end
        end


    end
    

end

