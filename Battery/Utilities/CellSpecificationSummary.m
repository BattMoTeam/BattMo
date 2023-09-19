classdef CellSpecificationSummary


    properties (SetAccess = private)
        % We want the computed values to remain synchronized with the model. To do so we set the "setAccess" property to
        % private. See function updateModel to change model and updatePackingMass to change the packingMass.
        model
        
        packingMass
        
        mass
        masses % structure with mass of each of the cell components

        volume
        volumes % structure with volume of each of the cell components
        
        energy
        specificEnergy
        energyDensity

        NPratio

        capacity
        capacities % structure with capacity of negative and positive electrode
        
        initialVoltage

        dischargeFunction % Maximum energy discharge function (voltage versus state of charge)

        temperature

        dischargeSimulations % cell array with struct elements with field
                             % - CRate
                             % - energy
                             % - specificEnergy
                             % - energyDensity
                             % - dischargeCurve
        
    end

    methods

        function css = CellSpecificationSummary(model, varargin)

            opt = struct('packingMass', 0, ...
                         'temperature', 298);
            opt = merge_options(opt, varargin{:});

            css.packingMass          = opt.packingMass;
            css.temperature          = opt.temperature;
            css.dischargeSimulations = {};
            css.model                = model;

            css = css.computeSpecs();
            
        end

        function css = updateModel(css, model)

            css.model = model;
            css = css.computeSpecs();
        end

        function css = updatePackingMass(css, packingMass)

            css.packingMass = packingMass;
            css = css.computeSpecs();
            
        end
            
        
        function css = computeSpecs(css)

            % Reset the simulations 
            css.dischargeSimulations = {};

            model = css.model;
            
            temperature = css.temperature;
            packingMass = css.packingMass;
            
            [mass, masses, volumes]     = computeCellMass(model, 'packingMass', packingMass);
            [capacity, capacities]      = computeCellCapacity(model);
            [energy, dischargeFunction] = computeCellEnergy(model, 'temperature', temperature);
            
            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';
            am  = 'ActiveMaterial';
            itf = 'Interface';
            sd  = 'SolidDiffusion';

            % Compute specific energy and energy density

            volume = volumes.val;

            specificEnergy = energy/mass;
            energyDensity  = energy/volume;
            
            % Compute NP ratio
            
            NPratio = capacities.(ne)/capacities.(pe);

            % Compute initial voltage
            
            itfmodel = model.(ne).(am).(itf);
            cmax  = itfmodel.cmax;
            cinit = itfmodel.theta100*cmax;
            
            U = -itfmodel.computeOCPFunc(cinit, temperature, cmax);
            
            itfmodel = model.(pe).(am).(itf);
            cmax  = itfmodel.cmax;
            cinit = itfmodel.theta100*cmax;
            
            U = U + itfmodel.computeOCPFunc(cinit, temperature, cmax);

            %  Assign values
            
            css.mass              = mass;
            css.masses            = masses;
            css.volume            = volume;
            css.volumes           = volumes;
            css.energy            = energy;
            css.specificEnergy    = specificEnergy;
            css.energyDensity     = energyDensity;
            css.NPratio           = NPratio;
            css.capacity          = capacity;
            css.capacities        = capacities;
            css.initialVoltage    = U;
            css.dischargeFunction = dischargeFunction;
        
        end

        function printSpecifications(css)

            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';

            function lines = addLine(lines, description, unit, value)

                line.description = description;
                line.unit        = unit;
                line.value       = value;

                lines{end + 1} = line;
                
            end
            
            lines = {};

            lines = addLine(lines, 'Packing mass'               , 'kg'   , css.packingMass);
            lines = addLine(lines, 'Temperature'                , 'C'    , css.temperature -  273.15);
            lines = addLine(lines, 'Mass'                       , 'kg'   , css.mass);
            lines = addLine(lines, 'Volume'                     , 'L'    , css.volume/litre);
            lines = addLine(lines, 'Total Capacity'             , 'Ah'   , css.capacity/hour);
            lines = addLine(lines, 'Negative Electrode Capacity', 'Ah'   , css.capacities.(ne)/hour);
            lines = addLine(lines, 'Positive Electrode Capacity', 'Ah'   , css.capacities.(pe)/hour);
            lines = addLine(lines, 'N/P ratio'                  , '-'    , css.NPratio);
            lines = addLine(lines, 'Energy'                     , 'Wh'   , css.energy/hour);
            lines = addLine(lines, 'Specific Energy'            , 'Wh/kg', css.specificEnergy/hour);
            lines = addLine(lines, 'Energy Density'             , 'Wh/L' , (css.energyDensity/hour)*litre);
            lines = addLine(lines, 'Initial Voltage'            , 'V'    , css.initialVoltage);

            function str = appendCrate(str, crate)

                str = sprintf('%s (CRate = %g)', str, crate);
                
            end
            
            for isim = 1 : numel(css.dischargeSimulations)
                
                simres = css.dischargeSimulations{isim};
                ac = @(str) appendCrate(str, simres.CRate); 
                lines = addLine(lines, ac('Energy'), 'Wh', simres.energy/hour);
                lines = addLine(lines, ac('Specific Energy'), 'Wh/kg', simres.specificEnergy/hour);
                lines = addLine(lines, ac('Energy Density') , 'Wh/L' , (simres.energyDensity/hour)*litre);
                
            end
            
            function printLines(lines)
            %% print lines
                
            % Setup format for fprintf
            % Get maximum string length for the description field
                descriptions = cellfun(@(line) line.description, lines, 'uniformoutput', false);
                lgths        = cellfun(@(str) length(str), descriptions);
                s            = max(lgths);

                fmt = sprintf('%%%ds : %%g [%%s]\n', s);

                for iline = 1 : numel(lines)
                    line = lines{iline};
                    fprintf(fmt             , ...
                            line.description, ...
                            line.value      , ...
                            line.unit);
                end
                
            end
            
            printLines(lines);
            
        end

        function css = addCrates(css, CRates, varargin)

            for icrate = 1 : numel(CRates)

                CRate = CRates(icrate);
                css = css.addCrate(CRate, varargin{:});

            end
            
        end
        
        function css = addCrate(css, CRate, varargin)

            % the extras options are passed to computeCellEnergy
            extras = varargin;
            
            model = css.model;
            
            [energy, dischargeFunction] = computeCellEnergy(model, 'CRate', CRate, extras{:});

            specificEnergy = energy/css.mass;
            energyDensity  = energy/css.volume;

            simresult = struct('CRate'            , CRate            , ...
                               'energy'           , energy           , ...
                               'specificEnergy'   , specificEnergy   , ...
                               'energyDensity'    , energyDensity    , ...
                               'dischargeFunction', dischargeFunction);

            css.dischargeSimulations{end + 1} = simresult;
            
        end
            
        
    end

    methods (Static)


        
    end
    
end
