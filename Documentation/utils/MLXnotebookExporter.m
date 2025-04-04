classdef MLXnotebookExporter

    properties

        % list of registered notebooks
        notebooknames = {'part_1_battery_modeling_guide.mlx'                   , ...
                         'part_2_battery_modeling_guide.mlx'                   , ...
                         'tutorial_1_a_simple_p2d_model_live.mlx'              , ...
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

        function exportMLX(mne, notebookname, varargin)
            
            opt = struct('format', 'm', ...
                         'run', false);
            [opt, extra] = merge_options(opt, varargin{:});

            ext = ['.', mne.getExtension(opt.format)];
            
            inputfile = fullfile(mne.inputdir, notebookname);
            
            notebookname = strrep(notebookname, '.mlx', ext);
            outputfile = fullfile(mne.outputdir, notebookname);
            
            optvals = cellfun(@(x) opt.(x), fieldnames(opt), 'uniformoutput', false);
            options = reshape(vertcat(fieldnames(opt)', optvals'), [], 1);
            
            export(inputfile, outputfile, options{:}, extra{:});

        end
        
        function exportMLXtoMfile(mne, notebookname, varargin)
            
            mne.exportMLX(notebookname, 'format', 'm', varargin{:});
            
        end

        function exportMLXtoHtml(mne, notebookname, varargin)

            mne.exportMLX(notebookname, 'format', 'html', varargin{:});

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

        function exportAll(mne)
            
            for inote = 1 : numel(mne.notebooknames)

                notebookname = mne.notebooknames{inote};

                fprintf('Starting with %s ...', notebookname);

                mne.exportMLXtoHtml(notebookname);

                fprintf('Done.\n', notebookname);
                
            end
        end


    end

    methods (Static)

        function ext = getExtension(format)
            switch format
              case {'m', 'html', 'ipynb'}
                ext = format;
              otherwise
                error('Unknown format %s', format)

            end
            
        end

    end

end

