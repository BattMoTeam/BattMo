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

        Temperature
        
    end

    methods

        function css = CellSpecificationSummary(varargin)

            opt = struct('packingMass', [],
                         'Temperature', 298);
            opt = merge_options(opt, varargin{:});

            css.packingMass = opt.packingMass;
            css.Temperature = opt.Temperature;
            
            css = css.update(model);
            
        end
        
        function css = update(css, model)

            T           = css.Temperature;
            packingMass = css.packingMass;
            
            [mass, masses, volumes] = computeCellMass(model, packingMass);
            [capacity, capacities]  = computeCellCapacity(model);
            energy                  = computeCellEnergy(model, 'Temperature', T);
            
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
            
            U = itfmodel.computeOCPFunc(cinit, T, cmax);
            
            itfmodel = model.(pe).(am).(itf);
            cmax  = itfmodel.cmax;
            cinit = itfmodel.theta100*cmax;
            
            U = U - itfmodel.computeOCPFunc(cinit, T, cmax);

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
            css.initialVoltage    = initialVoltage;
            css.dischargeFunction = dischargeFunction;
            css.Temperature       = Temperature;
        
        end

        function printSpecifications(css)
            
        end
    end
end
