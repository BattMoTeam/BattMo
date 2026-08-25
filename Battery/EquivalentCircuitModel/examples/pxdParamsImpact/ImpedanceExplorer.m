classdef ImpedanceExplorer < handle
    properties (SetAccess = immutable)
        jsonstruct % Holds the base P2D model configuration structure
        model      % Initialized BattMo model object
    end
    
    properties
        soc = 1 % Default State of Charge (0 to 1)
        parameters % Cell array storing the hierarchical JSON paths of observed parameters
    end
    
    methods
        function impexp = ImpedanceExplorer(jsonstruct)
            % Constructor: auto-loads Chen2020 P2D model if no structure is provided
            if nargin < 1 || isempty(jsonstruct)
                fprintf('Initializing default Chen2020 P2D model...\n');
                try
                    jsonstruct_material = parseBattmoJson(fullfile('ParameterData', 'ParameterSets', 'Chen2020', 'chen2020_lithium_ion_battery.json'));
                    jsonstruct_geometry = parseBattmoJson(fullfile('Examples', 'JsonDataFiles', 'geometryChen.json'));
                    jsonstruct = mergeStructs({jsonstruct_material, jsonstruct_geometry});
                    
                    % Safe default fallback for capacitance
                    ne = 'NegativeElectrode'; co = 'Coating'; am = 'ActiveMaterial'; itf = 'Interface';
                    if ~isfield(jsonstruct.(ne).(co).(am).(itf), 'doubleLayerCapacitance')
                        jsonstruct.(ne).(co).(am).(itf).useDoubleLayerCapacity = true;
                        jsonstruct.(ne).(co).(am).(itf).doubleLayerCapacitance = 0.2;
                    end
                catch ME
                    error('ImpedanceExplorer:AutoLoadFailed', 'Failed to auto-load default files: %s', ME.message);
                end
            end
            
            impexp.jsonstruct = jsonstruct;
            impexp.model      = setupModelFromJson(jsonstruct);
            impexp.setupDefaultParameters();
        end
        
        function viewOptions(impexp)
            % Helper to print all available scalar properties in the model JSON
            fjv = flattenStruct(impexp.jsonstruct, 'doprint', false);
            fprintf('\n\n*** List of all the scalar parameters in the model ***\n');
            fjv.print('filter', {'value', @(val) (isnumeric(val) && isscalar(val))});
        end
            
        function setupDefaultParameters(impexp)
            % Defines the default nested JSON paths for the 3 target sliders
            ne  = 'NegativeElectrode';
            co  = 'Coating';
            am  = 'ActiveMaterial';
            itf = 'Interface';
            sd  = 'SolidDiffusion';
 
            impexp.parameters = {{ne, co, am, sd,  'referenceDiffusionCoefficient'}, ...
                                 {ne, co, am, itf, 'reactionRateConstant'}, ...
                                 {ne, co, am, itf, 'doubleLayerCapacitance'}};
        end
        
        function pnames = setupParameterLegendNames(impexp)
            % Joins path cells into unified text labels for the GUI
            params = impexp.parameters;
            pnames = cell(1, numel(params));
            for iparam = 1 : numel(params)
                pnames{iparam} = params{iparam}{end};
            end
        end
        
        function sld = setupSlider(impexp, fig, ypos, txt, limits, majorTicks, initVal)
            % Generates a single slider with its associated label
            if nargin < 5, limits = [-5, 5]; end
            if nargin < 6, majorTicks = -5:0.5:5; end
            if nargin < 7, initVal = 0; end
            
            uilabel(fig, ...
                    'Position', [50, ypos + 10, 500, 40], ...
                    'WordWrap', 'on', ...
                    'Text', txt);
            
            sld = uislider(fig, ...
                           'Position', [50, ypos, 500, 3], ...
                           'Limits', limits, ...
                           'MajorTicks', majorTicks, ...
                           'Value', initVal);
        end
        
        function vals = getParameterValues(impexp, jsonstruct)
            % Reads parameter values from JSON structure. Error out if path is invalid.
            if nargin < 2
                jsonstruct = impexp.jsonstruct;
            end
            
            parameters = impexp.parameters;
            vals = zeros(1, numel(parameters)); 
            for iparam = 1 : numel(parameters)
                val = getStructField(jsonstruct, parameters{iparam});
                if ~isnumeric(val) || ~isscalar(val)
                    pathStr = strjoin(parameters{iparam}, '.');
                    error('ImpedanceExplorer:ParameterNotFound', ...
                        ['\n\n=== PARAMETER NOT FOUND ===\n' ...
                         'Invalid path in BattMo JSON:\n -> "jsonstruct.%s"\n' ...
                         'Check spelling or use app.viewOptions() to debug.\n\n'], pathStr);
                end
                vals(iparam) = val;
            end
        end
        
        function jsonstruct = setParameterValues(impexp, vals, jsonstruct)
            % Writes updated parameter values back into the JSON structure
            if nargin < 3
                jsonstruct = impexp.jsonstruct;
            end
            
            parameters = impexp.parameters;
            for iparam = 1 : numel(parameters)
                jsonstruct = setStructField(jsonstruct, parameters{iparam}, vals(iparam), 'handleMisMatch', 'quiet');
            end
        end
        
        function start(impexp)
            
            % Computes sizes and runs the App Designer UI Window
            pnames = impexp.setupParameterLegendNames();
            vals   = impexp.getParameterValues();
            
            numParams = numel(pnames);
            numSliders = numParams + 1; % Params + SoC slider
            
            % Layout geometry constants
            sliderSpacing = 75;
            bottomOffset  = 40;
            
            lblY = bottomOffset + numSliders * sliderSpacing + 10;
            axY  = lblY + 35;
            figHeight = axY + 250 + 40; 
            
            fig = uifigure('Name', 'Impedance Explorer', 'Position', [100, 50, 600, figHeight]);
            
            % Nyquist plot setup
            ax = uiaxes(fig, 'Position', [50, axY, 500, 250]);
            title(ax, 'Nyquist Diagram');
            xlabel(ax, 'Re(Z)');
            ylabel(ax, '-Im(Z)');
            grid(ax, 'on');
            
            % Status bar setup
            lbl = uilabel(fig, 'Position', [50, lblY, 500, 22], 'Text', 'Ready, adjust a slider to run simulation.');
            lbl.FontWeight = 'bold';
            
            % Create parameter sliders
            sliders = cell(1, numParams);
            for iparam = 1 : numParams
                txt = sprintf('%d) %s * 10^x (initial: %g)', iparam, pnames{iparam}, vals(iparam));
                ypos = bottomOffset + (numParams - iparam + 1)*sliderSpacing; 
                sliders{iparam} = impexp.setupSlider(fig, ypos, txt);
            end
            
            % Create SoC slider (at the absolute bottom)
            socTxt = sprintf('State of Charge (SoC)');
            socSlider = impexp.setupSlider(fig, bottomOffset, socTxt, [0, 1], 0:0.1:1, impexp.soc);

            impexp.updatePlot(ax, lbl, sliders, socSlider);
            
            % Wire slider events to the updatePlot method
            callback_fcn = @(src, event) impexp.updatePlot(ax, lbl, sliders, socSlider);
            for islider = 1 : numParams
                sliders{islider}.ValueChangedFcn = callback_fcn;
            end
            socSlider.ValueChangedFcn = callback_fcn;
            
        end
        
        function updatePlot(impexp, ax, lbl, sliders, socSlider)
            % Event handler to recalculate physics and update the plot
            jsonstruct = impexp.jsonstruct;
            lbl.Text = 'Computing EIS... Please wait.';
            lbl.FontColor = '#D95319';
            drawnow; 
            
            try
                % Read updated GUI values
                impexp.soc = socSlider.Value;
                baseValues = impexp.getParameterValues();
                vals = zeros(size(baseValues));
                for ival = 1 : numel(baseValues)
                    vals(ival) = baseValues(ival) * (10^(sliders{ival}.Value));
                end
                jsonstruct = impexp.setParameterValues(vals);
                
                % Setup BattMo Model and state initialization
                [~, inputparams, ~] = setupModelFromJson(jsonstruct);
                options = [];
                options.stateInitialization.initializationSetup = 'soc';
                options.stateInitialization.soc = impexp.soc;
                
                % Execute Impedance Solver
                impsolv = ImpedanceSolver(inputparams, options);
                frequencies = logspace(-2, 3, 100); 
                Z = impsolv.computeImpedance(frequencies);
                
                % Draw new Nyquist trace
                plot(ax, real(Z), -imag(Z), '-', 'LineWidth', 2);
                
                lbl.Text = 'Done';
                lbl.FontColor = '#77AC30';
                
            catch ME
                
                lbl.Text = 'Calculation failed.';
                lbl.FontColor = '#A2142F';
                fprintf('\n--- BATTMO CALCULATION ERROR ---\n%s\n', ME.message);
                
            end
        end
    end
end
