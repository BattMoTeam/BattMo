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
            
            U = itfmodel.computeOCPFunc(cinit, temperature, cmax);
            
            itfmodel = model.(pe).(am).(itf);
            cmax  = itfmodel.cmax;
            cinit = itfmodel.theta100*cmax;
            
            U = U - itfmodel.computeOCPFunc(cinit, temperature, cmax);

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
            
        end
    end
end
