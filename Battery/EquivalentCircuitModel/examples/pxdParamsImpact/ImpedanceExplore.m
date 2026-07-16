classdef ImpedanceExplore < handle
    properties (SetAccess = immutable)
        jsonstruct
        model
        
    end
    
    properties
        
        
        % done : add option to change soc in slider
        soc = 1
        
        parameters
        
    end
    
    methods
        %%
        % Constructor method that initializes the application, auto-loading the physical P2D Chen2020 battery dataset from BattMo files if no input JSON structure is provided.
        function impexp = ImpedanceExplore(jsonstruct)
            if nargin < 1 
                fprintf('Automatic initialization of Chen2020 configuration (BattMo)...\n');
                try
                    jsonstruct_material = parseBattmoJson(fullfile('ParameterData', 'ParameterSets', 'Chen2020', 'chen2020_lithium_ion_battery.json'));
                    jsonstruct_geometry = parseBattmoJson(fullfile('Examples', 'JsonDataFiles', 'geometryChen.json'));
                    
                    jsonstruct = mergeJsonStructs({jsonstruct_material, jsonstruct_geometry});
                    if ~isfield(jsonstruct.NegativeElectrode.Coating.ActiveMaterial.Interface, 'doubleLayerCapacitance')
                        jsonstruct.NegativeElectrode.Coating.ActiveMaterial.Interface.doubleLayerCapacitance = 0.2; 
                        fprintf('-> default doubleLayerCapacitance, value = 0.2 F/m^2\n');
                    end
                    
                catch ME
                    error('ImpedanceExplore:AutoLoadFailed', ...
                        ME.message);
                end
            end
            model = setupModelFromJson(jsonstruct);
            impexp.jsonstruct = jsonstruct;
            impexp.model      = model;
            
            impexp.setupDefaultParameters()
            
        end
        
        %%
        % Utility method to flatten the active JSON structure and print all available scalar parameters in the model to the command window.
        function viewOptions(impexp)
            fjv = flattenJsonStruct(impexp.jsonstruct, 'doprint', false);
            fprintf('\n\n*** List of all the scalar parameters in the model ***\n');
            fjv.print('filter', {'value', @(val) (isnumeric(val) && isscalar(val))});
            
        end
        
        %%
        % Configures the default internal hierarchical paths of the active material parameters to be observed and modified.
        function setupDefaultParameters(impexp)
            elyte   = 'Electrolyte';
            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            co      = 'Coating';
            am      = 'ActiveMaterial';
            itf     = 'Interface';
            sd      = 'SolidDiffusion';
            cc      = 'CurrentCollector';
            ctrl    = 'Control';
            thermal = 'ThermalModel';
 
            impexp.parameters = {{ne, co, am, sd,  'referenceDiffusionCoefficient'}, ...
                                 {ne, co, am, itf, 'reactionRateConstant'}, ...
                                 {ne, co, am, itf, 'doubleLayerCapacitance'}};
        end
        
        %%
        % Helper method that generates human-readable parameter names by joining the individual sub-keys of their nested JSON paths.
        function pnames = setupParameterLegendNames(impexp)
            params = impexp.parameters;
            for iparam = 1 : numel(params)
                pnames{iparam} = strjoin(params{iparam});
            end
            
        end
        
        %%
        % Dynamically instantiates a slider control alongside its descriptive label at a specific vertical coordinate.
        function sld = setupSlider(impexp, fig, ypos, txt, limits, majorTicks, initVal)
            if nargin < 5, limits = [-5, 5]; end
            if nargin < 6, majorTicks = -5:0.5:5; end
            if nargin < 7, initVal = 0; end
            
            uilabel(fig, ...
                    'Position', [50, ypos + 20, 500, 40], ...
                    'WordWrap', 'on', ...
                    'Text', txt);
            
            sld = uislider(fig, ...
                           'Position', [50, ypos, 500, 3], ...
                           'Limits', limits, ...
                           'MajorTicks', majorTicks, ...
                           'Value', initVal);
        end
        
        %%
        % Retrieves the current numerical values of targeted active parameters from the application''s JSON model.
        function vals = getParameterValues(impexp, jsonstruct)
            if nargin < 2
                jsonstruct = impexp.jsonstruct;
            end
            
            parameters = impexp.parameters;
            vals = zeros(1, numel(parameters)); 
            for iparam = 1 : numel(parameters)
                val = getJsonStructField(jsonstruct, parameters{iparam});
                if ~isnumeric(val) || ~isscalar(val)
                pathStr = strjoin(parameters{iparam}, '.');
                error('ImpedanceExplore:ParameterNotFound', ...
                    ['\n\n=== PARAMETER NOT FOUND ===\n' ...
                     'The following path is incorrect in the BattMo JSON:\n' ...
                     ' -> "jsonstruct.%s"\n' ...
                     'Received type: %s (instead of a double).\n\n' ...
                     'To find the correct path, use the command:\n' ...
                     ' >> app = ImpedanceExplore(); app.viewOptions();\n\n'], ...
                    pathStr, class(val));
                end
            
                vals(iparam) = val;
            end
            
        end
        
        %%
        % Updates specific nested fields within the model''s JSON structure using values calculated from the current slider states.
        function jsonstruct = setParameterValues(impexp, vals, jsonstruct)
            if nargin < 3
                jsonstruct = impexp.jsonstruct;
            end
            
            parameters = impexp.parameters;
            for iparam = 1 : numel(parameters)
                jsonstruct = setJsonStructField(jsonstruct, parameters{iparam}, vals(iparam), 'handleMisMatch', 'quiet');
            end
            
        end
        
        %%
        % Dynamically calculates element positions, builds the graphical user interface, spawns the sliders, and registers event callbacks.
        function start(impexp)
            pnames = impexp.setupParameterLegendNames();
            vals   = impexp.getParameterValues();
            
            numParams = numel(pnames);
            numSliders = numParams + 1; % active parameters + SoC slider
            
            % --- DYNAMIC GEOMETRIC CONFIGURATION ---
            sliderSpacing = 75;    % Height allocated to each block (Label + Slider)
            bottomOffset  = 40;    % Margin at the bottom of the window
            
            % Calculate Y positions for status label and axes
            lblY = bottomOffset + numSliders * sliderSpacing + 10;
            axY  = lblY + 35;
            
            % Dynamic calculation of the window height
            % (axY + 250 pixels for the plot + 40 pixels margin at the top)
            figHeight = axY + 250 + 40; 
            
            % 1. Create the user interface at the correct size
            fig = uifigure('Name', 'Impedance Explorer', 'Position', [100, 50, 600, figHeight]);
            
            % 2. Create automatically positioned axes
            ax = uiaxes(fig, 'Position', [50, axY, 500, 250]);
            title(ax, 'Nyquist Diagram');
            xlabel(ax, 'Re(Z)');
            ylabel(ax, '-Im(Z)');
            grid(ax, 'on');
            
            % 3. Create automatically positioned status label
            lbl = uilabel(fig, 'Position', [50, lblY, 500, 22], 'Text', 'Ready, use slider to start.');
            lbl.FontWeight = 'bold';
            
            sliders = cell(1, numParams);
            
            % --- CREATE PARAMETER SLIDERS ---
            % Stack them from bottom to top. The first parameter will be placed right above the SoC slider.
            for iparam = 1 : numParams
                txt = sprintf('%d. %s * 10^x (initial value: %g)', iparam, pnames{iparam}, vals(iparam));
                
                % Calculate Y position: offsets by "iparam" heights
                ypos = bottomOffset + iparam * sliderSpacing; 
                
                sliders{iparam} = impexp.setupSlider(fig, ypos, txt);
            end
            
            % --- CREATE SOC SLIDER (At the very bottom) ---
            socTxt = sprintf('State of Charge (SoC) - Current: %g', impexp.soc);
            ypos_soc = bottomOffset; % Starting position at the very bottom
            socSlider = impexp.setupSlider(fig, ypos_soc, socTxt, [0, 1], 0:0.1:1, impexp.soc);
            
            % --- EVENT CONNECTIONS ---
            callback_fcn = @(src, event) impexp.updatePlot(ax, lbl, sliders, socSlider);
            
            % Bind callback to all parameter sliders
            for islider = 1 : numParams
                sliders{islider}.ValueChangedFcn = callback_fcn;
            end
            % Bind callback to the SoC slider
            socSlider.ValueChangedFcn = callback_fcn;
        end
        
        %%
        % Core callback method that retrieves current slider values, updates the physical state, runs the BattMo impedance solver, and updates the Nyquist plot.
        function updatePlot(impexp, ax, lbl, sliders, socSlider)
            jsonstruct = impexp.jsonstruct;
            % Freeze display
            lbl.Text = 'Processing... ';
            lbl.FontColor = '#D95319';
            drawnow; 
            try
                impexp.soc = socSlider.Value;
                baseValues = impexp.getParameterValues();
                for ival = 1 : numel(baseValues)
                    vals(ival) = baseValues(ival)*(10^(sliders{ival}.Value));
                end
                jsonstruct = impexp.setParameterValues(vals);
                
                % 4. Execute BattMo
                [model, inputparams, ~] = setupModelFromJson(jsonstruct);
                %%
                % We use a dedicated function to compute the initial state for the given electrode concentrations.
                options = [];
                options.stateInitialization.initializationSetup = 'soc';
                options.stateInitialization.soc = impexp.soc;  % now the soc is updated
                %%
                % Setup impedance solver and compute the Impedance for a given set of frequencies
                impsolv = ImpedanceSolver(inputparams, options);
                frequences = logspace(-2, 3, 30); 
                Z = impsolv.computeImpedance(frequences);
                
                % 5. Plotting
                plot(ax, real(Z), -imag(Z), '-o', 'LineWidth', 2, 'Color', '#0072BD');
                
                % Return to normal
                lbl.Text = 'Done';
                lbl.FontColor = '#77AC30';
            catch ME
                lbl.Text = 'Calculation error';
                lbl.FontColor = '#A2142F';
                disp('--- Calculation error ---');
                disp(ME.message);
                for k = 1:length(ME.stack)
                    disp(['Line ', num2str(ME.stack(k).line), ' in ', ME.stack(k).name]);
                end
            end
        end
    end
end