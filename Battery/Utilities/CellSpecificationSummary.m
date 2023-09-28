classdef CellSpecificationSummary

    properties
        % Note that all values are in SI (standard units)
        
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
        
    end

    methods

        function css = CellSpecificationSummary(model, varargin)

            opt = struct('packingMass', 0, ...
                         'temperature', 298);
            opt = merge_options(opt, varargin{:});

            css.packingMass = opt.packingMass;
            css.temperature = opt.temperature;
            
            css = css.update(model);
            
        end
        
        function css = update(css, model)

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

            lines = {};
            
            lines = css.addLine(lines, 'Packing mass'               , 'kg'   , css.packingMass);
            lines = css.addLine(lines, 'Temperature'                , 'C'    , css.temperature -  273.15);
            lines = css.addLine(lines, 'Mass'                       , 'kg'   , css.mass);
            lines = css.addLine(lines, 'Volume'                     , 'L'    , css.volume/litre);
            lines = css.addLine(lines, 'Total Capacity'             , 'Ah'   , css.capacity/hour);
            lines = css.addLine(lines, 'Negative Electrode Capacity', 'Ah'   , css.capacities.(ne)/hour);
            lines = css.addLine(lines, 'Positive Electrode Capacity', 'Ah'   , css.capacities.(pe)/hour);
            lines = css.addLine(lines, 'N/P ratio'                  , '-'    , css.NPratio);
            lines = css.addLine(lines, 'Energy'                     , 'Wh'   , css.energy/hour);
            lines = css.addLine(lines, 'Specific Energy'            , 'Wh/kg', css.specificEnergy/hour);
            lines = css.addLine(lines, 'Energy Density'             , 'Wh/L' , (css.energyDensity/hour)*litre);
            lines = css.addLine(lines, 'Initial Voltage'            , 'V'    , css.initialVoltage);
            
            css.printLines(lines);
            
        end
        
    end

    methods (Static)

        function lines = addLine(lines, description, unit, value)

            line.description = description;
            line.unit        = unit;
            line.value       = value;

            lines{end + 1} = line;
            
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
        
        
    end
    
end