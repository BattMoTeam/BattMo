classdef MLXnotebookExporter

    properties

        % list of the registered mlx notebooks (for info and automatic iteration)
        notebooknames = {'part_1_battery_modeling_guide'                   , ...
                         'part_2_battery_modeling_guide'                   , ...
                         'tutorial_1_a_simple_p2d_model_live'              , ...
                         'tutorial_2_changing_control_protocol_live'       , ...
                         'tutorial_3_modify_structural_parameters_live'    , ...
                         'tutorial_4_modify_material_parameters_live'      , ...
                         'tutorial_5_simulate_CCCV_cycling_live'           , ...
                         'tutorial_6_simulate_thermal_performance_live'    , ...
                         'tutorial_7_a_simple_p4d_model_live'              , ...
                         'tutorial_8_simulate_a_multilayer_pouch_cell_live', ...
                         'tutorial_9_simulate_a_cylindrical_cell_live'};

        % list of registered m-scripts (obtained from the test suite)
        mscripts
        
        % input directories for the mlx notebooks (for info)
        inputdir  = fullfile(battmoDir(), 'Examples', 'Notebooks');
        outputdir = fullfile(battmoDir(), 'Documentation', '_static', 'notebooks');
        
    end

    methods

        function mne = MLXnotebookExporter()
            
            testrunexample = TestRunExamples();
            mne.mscripts = testrunexample.filename;
            
        end

        function updateIpynbs(mne)
        % Update all the ipynb in the documentation.

            run_note_book = false;
            
            inputdir  = mne.inputdir;
            outputdir = mne.outputdir;
            
            for inote = 1 : numel(mne.notebooknames)
                
                notebookname = mne.notebooknames{inote};
                
                inputfilename  = fullfile(inputdir, [notebookname, '.mlx']);
                outputfilename = fullfile(outputdir, [notebookname, '.ipynb']);
                
                export(inputfilename, outputfilename, 'format', 'ipynb', 'Run', run_note_book);
                
            end

        end

        function convertMlxToM(mne)
        % Convert all the registered mlx notebooks to m script
            
            inputdir  = mne.inputdir;

            for inote = 1 : numel(mne.notebooknames)
                
                notebookname = mne.notebooknames{inote};
                inputfilename  = fullfile(inputdir, [notebookname, '.mlx']);
                outputfilename = fullfile(inputdir, [notebookname, '.m']);
                export(inputfilename, outputfilename, 'format', 'm');
                
            end
        end
        
        function setupIpynbFromMlx(mne, filename, varargin)

            opt = struct('outputDirectory', []);
            opt = merge_options(opt, varargin{:});

            [inputfile, outputfile] = getIOfiles(mne, filename, ...
                                                 'outputDirectory', opt.outputDirectory, ...
                                                 'outputFormat', 'ipynb');

            export(inputfile, outputfile, 'format', 'ipynb');
            
        end

        
        function setupMlxFromM(mne, filename, varargin)

            opt = struct('run'            , false, ...
                         'outputDirectory', []);
            opt = merge_options(opt, varargin{:});

            [inputfile, outputfile] = getIOfiles(mne, filename, ...
                                                 'outputDirectory', opt.outputDirectory, ...
                                                 'outputFormat', 'mlx');
            
            matlab.internal.liveeditor.openAndSave(inputfile, outputfile);

            if opt.run
                matlab.internal.liveeditor.executeAndSave(outputfile);
            end
            
        end
        
        function setupIpynbFromM(mne, filename, varargin)

        % From a M-file, generate the notebooks (mlx and ipynb)
        % The default directory for the generated notebooks is a sub-directory where the M-file is located called 'notebooks'

            opt = struct('run'            , false, ...
                         'outputDirectory', []   , ...
                         'generateIpynb'  , true);
            opt = merge_options(opt, varargin{:});

            mne.setupMlxFromM(filename, ...
                              'run'            , opt.run, ...
                              'outputDirectory', opt.outputDirectory);

            if opt.generateIpynb
                mne.setupIpynbFromMlx(filename, 'outputDirectory', opt.outputDirectory);
            end
            
        end

        function [inputfile, outputfile] = getIOfiles(mne, filename, varargin)
            
            opt = struct('outputDirectory', [], ...
                         'outputFormat', 'mlx');
            
            opt = merge_options(opt, varargin{:});

            assert(exist(filename, 'file') == 2, 'File %s not found.', filename); 
            fullfilename = which(filename);
            [dirpath, filename, ext] = fileparts(fullfilename);

            if isempty(opt.outputDirectory)
                outputDirectory = fullfile(dirpath, 'notebooks');
                if ~exist(outputDirectory)
                    % Create the directory if it does not exist
                    mkdir(outputDirectory);
                end
            else
                outputDirectory = opt.outputDirectory
            end

            switch opt.outputFormat
              case 'mlx'
                inputfile  = fullfile(dirpath, [filename, '.m']);
                outputfile = fullfile(outputDirectory, [filename, '.mlx']);
              case 'ipynb'
                inputfile  = fullfile(outputDirectory, [filename, '.mlx']);
                outputfile = fullfile(outputDirectory, [filename, '.ipynb']);
              otherwise
                error('outputFormat not recognized');
            end
            
        end

        function runMlxAndSave(mne, filename)

        % To run and update the mlx notebook programmatically, it is possible to use:
        % matlab.internal.liveeditor.executeAndSave('fullpathnameto.mlx')
            
            matlab.internal.liveeditor.executeAndSave(filename);

        end

    end


end

