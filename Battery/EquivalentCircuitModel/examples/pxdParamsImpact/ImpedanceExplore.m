classdef ImpedanceExplore < handle

    properties (SetAccess = private)

        jsonstruct
        
    end
    
    properties
        
        %%
        % TODO: add option to change soc in slider
        soc = 1
        
        parameters
        parameterLegendNames
        
    end
    
    methods

        function impexp = ImpedanceExplore(jsonstruct)

            impexp.jsonstruct = jsonstruct;
            impexp.setupDefaultParameters()
            
        end

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
            
            %%
            % TODO : setup default way to get legend
            impexp.parameterLegendNames = {'Negative Electrode Solid Diffusion Coefficient', ...
                                           'Negative Electrode Reaction Rate Constant', ...
                                           'Negative Electrode Double Layer Capacitance'};

        end

        function setupDefaultParameterLegendNames(impexp)

            params = impexp.parameters;
            for iparam = 1 : numel(params)
                impexp.parameterLegendNames{iparam} = strjoin(params{iparam});
            end
            
        end
        
        function updateImpedance(impexp)
            
        end

        function sld = setupSlider(imexp, fig, ypos, txt)
            
            uilabel(fig                                 , ...
                    'Position', [50, ypos + 20, 500, 40], ...
                    'WordWrap', 'on'                    , ...
                    'Text', txt);
            

            sld = uislider(fig, ...
                           'Position', [50, ypos, 500, 3], ...
                           'Limits', [-5, 5]             , ... % Définit le minimum et le maximum
                           'MajorTicks', -5:0.5:5        , ... % Ajoute une graduation tous les 1
                           'Value', 0);
                
        end

        function vals = getParameterValues(impexp, jsonstruct)

            if nargin < 2
                jsonstruct = impexp.jsonstruct;
            end
            
            parameters = impexp.parameters;

            for iparam = 1 : numel(parameters)
                vals(iparam) = getJsonStructField(jsonstruct, parameters{iparam});
            end
            
        end
        
        function jsonstruct = setParameterValues(impexp, vals, jsonstruct)

            if nargin < 3
                jsonstruct = impexp.jsonstruct;
            end
            
            parameters = impexp.parameters;

            for iparam = 1 : numel(parameters)
                jsonstruct = setJsonStructField(jsonstruct, parameters{iparam}, vals(iparam), 'handleMisMatch', 'quiet');
            end
            
        end

        function start(impexp)
            
        % 2. Création de l'interface graphique (Fenêtre plus haute pour accueillir les curseurs)
            fig = uifigure('Name', 'Impedance Explorer', 'Position', [100, 50, 600, 650]);
            
            % Graphique
            ax = uiaxes(fig, 'Position', [50, 350, 500, 250]);
            title(ax, 'Nyquist Diagram');
            xlabel(ax, 'Re(Z)');
            ylabel(ax, '-Im(Z)');
            grid(ax, 'on');
            
            % Label de statut
            lbl = uilabel(fig, 'Position', [50, 310, 500, 22], 'Text', 'Ready, use slider to start.');
            lbl.FontWeight = 'bold';

            sliders = {};

            pnames = impexp.parameterLegendNames;
            vals   = impexp.getParameterValues();
            % TODO : setup logic to get y positions based on number of parameters
            yposs  = [240, 150, 60];
            
            for iparam = 1 : numel(pnames)
                txt = sprintf('%d. %s * 10^x (initial value: %g)', iparam, pnames{iparam}, vals(iparam));
                sliders{iparam} = impexp.setupSlider(fig, yposs(iparam), txt);
            end

            % 4. Connexion aux événements
            % On crée une fonction anonyme qui passe TOUS les curseurs
            callback_fcn = @(src, event) impexp.updatePlot(ax, lbl, sliders);
            for islider = 1 : numel(sliders)
                sliders{islider}.ValueChangedFcn = callback_fcn;
            end
            % Optionnel : Tracer l'état initial
            % updatePlot(ax, lbl, jsonstruct, c_ne, c_pe, sld1, sld2, sld3);
            
        end

        function updatePlot(impexp, ax, lbl, sliders)

            jsonstruct = impexp.jsonstruct;
            
            % Figer l'affichage
            lbl.Text = 'Processing... ';
            lbl.FontColor = '#D95319';
            drawnow; 

            try

                baseValues = impexp.getParameterValues();
                
                for ival = 1 : numel(baseValues)
                    vals(ival) = baseValues(ival)*(10^(sliders{ival}.Value));
                end

                jsonstruct = impexp.setParameterValues(vals);
                
                % 4. Exécution de BattMo

                [model, inputparams, ~] = setupModelFromJson(jsonstruct);

                %%
                % We use a dedicated function to compute the initial state for the given electrode concentrations.
                options = [];
                options.stateInitialization.initializationSetup = 'soc';
                options.stateInitialization.soc = impexp.soc;

                %%
                % setup impedance solver and compute the Impedance for a given set of frequencies
                impsolv = ImpedanceSolver(inputparams, options);

                frequences = logspace(-2, 3, 30); 
                Z = impsolv.computeImpedance(frequences);
                
                % 5. Tracé
                plot(ax, real(Z), -imag(Z), '-o', 'LineWidth', 2, 'Color', '#0072BD');
                
                % Retour à la normale
                lbl.Text = 'Done';
                lbl.FontColor = '#77AC30';
                
            catch ME
                
                lbl.Text = 'Erreur lors du calcul (vérifiez la console).';
                lbl.FontColor = '#A2142F';
                disp('--- ERREUR DANS LE CALCUL ---');
                disp(ME.message);
                for k = 1:length(ME.stack)
                    disp(['Ligne ', num2str(ME.stack(k).line), ' dans ', ME.stack(k).name]);
                end
            end
        end
        
    end


end
