classdef CellSpecificationSummary
% Utility class to compute standard cell specifications (see list in properties below) using the model
%
% Energy computation for given DRate can be added using the method addDrate. The results will be stored in the property
% dischargeSimulations.
    
    properties (SetAccess = private)
        % We want the computed values to remain synchronized with the model. To do so we set the "setAccess" property to
        % private. See function updateModel to change model and updatePackingMass to change the packingMass.
        model
        
        packingMass
        thicknesses % Used to compute mass loading
        mass
        masses % structure with mass of each of the cell components

        massLoadings % Computed if thicknesses are available
        
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
                             % - DRate
                             % - energy
                             % - specificEnergy
                             % - energyDensity
                             % - dischargeCurve
                             % - E    % Voltage output 
                             % - I    % Current output
                             % - time % time output
        
    end

    methods

        function css = CellSpecificationSummary(model, varargin)

            opt = struct('packingMass', 0 , ...
                         'thicknesses', [], ...
                         'jsonstruct', [] , ...
                         'temperature', 298);
            opt = merge_options(opt, varargin{:});

            css.packingMass          = opt.packingMass;
            css.temperature          = opt.temperature;
            css.dischargeSimulations = {};
            css.model                = model;

            thicknesses = css.extractThicknessFromModel(model); % can be done in 1D model
            
            if ~isempty(opt.jsonstruct)
                % We fetch the thicknesses from the jsonstruct
                thicknesses.NegativeElectrode = opt.jsonstruct.NegativeElectrode.Coating.thickness;
                thicknesses.PositiveElectrode = opt.jsonstruct.PositiveElectrode.Coating.thickness;
            end

            if ~isempty(opt.thicknesses)
                % We fetch the thicknesses from the jsonstruct
                thicknesses = opt.thicknesses;
            end

            css.thicknesses = thicknesses;
            
            css = css.computeSpecs();
            
        end

        function css = updatePackingMass(css, packingMass)

            css.packingMass = packingMass;
            css = css.computeSpecs();
            
        end

        function thicknesses = extractThicknessFromModel(css, model)

            if model.grid.griddim > 1
                thicknesses = [];
                return
            end

            ne = 'NegativeElectrode';
            pe = 'PositiveElectrode';
            co = 'Coating';
            
            eldes = {ne, pe};

            for ielde = 1 : numel(eldes)

                elde = eldes{ielde};

                G = model.(elde).(co).grid;
                xmax = max(G.faces.centroids(:, 1));
                xmin = min(G.faces.centroids(:, 1));

                thicknesses.(elde) = xmax - xmin;
                
            end
            
        end

        
        function css = updateNegativeElectrodeThickness(css, thickness)
            
            css.thicknesses.NegativeElectrode = thickness;
            css = css.computeSpecs();
            
        end
        
        function css = updatePositiveElectrodeThickness(css, thickness)
            
            css.thicknesses.PositiveElectrode = thickness;
            css = css.computeSpecs();
            
        end

        function css = updateThicknesses(css, thicknesses)
            
            css.thicknesses = thicknesses;
            css = css.computeSpecs();
            
        end
        

        function css = computeSpecs(css)

            % Reset the simulations 

            css.dischargeSimulations = {};

            model = css.model;
            
            temperature = css.temperature;
            packingMass = css.packingMass;
            
            [mass, masses, volumes] = computeCellMass(model, 'packingMass', packingMass);
            [capacity, capacities]  = computeCellCapacity(model);
            [energy, output]        = computeCellEnergy(model, 'temperature', temperature);
            
            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';
            co  = 'Coating';
            am  = 'ActiveMaterial';
            itf = 'Interface';
            sd  = 'SolidDiffusion';

            % Compute Mass Loadings

            if ~isempty(css.thicknesses)
                massLoadings.(ne) = css.thicknesses.(ne)*model.(ne).(co).effectiveDensity;
                massLoadings.(pe) = css.thicknesses.(pe)*model.(pe).(co).effectiveDensity;
            else
                massLoadings.(ne) = [];
                massLoadings.(pe) = [];
            end
            
            % Compute specific energy and energy density

            volume = volumes.val;

            specificEnergy = energy/mass;
            energyDensity  = energy/volume;
            
            % Compute NP ratio
            
            NPratio = capacities.(ne)/capacities.(pe);

            % Compute initial voltage
            
            itfmodel = model.(ne).(co).(am).(itf);
            cmax  = itfmodel.saturationConcentration;
            cinit = itfmodel.guestStoichiometry100*cmax;
            
            U = -itfmodel.computeOCPFunc(cinit, temperature, cmax);
            
            itfmodel = model.(pe).(co).(am).(itf);
            cmax  = itfmodel.saturationConcentration;
            cinit = itfmodel.guestStoichiometry100*cmax;
            
            U = U + itfmodel.computeOCPFunc(cinit, temperature, cmax);

            %  Assign values
            
            css.mass              = mass;
            css.masses            = masses;
            css.volume            = volume;
            css.volumes           = volumes;
            css.massLoadings      = massLoadings;
            css.energy            = energy;
            css.specificEnergy    = specificEnergy;
            css.energyDensity     = energyDensity;
            css.NPratio           = NPratio;
            css.capacity          = capacity;
            css.capacities        = capacities;
            css.initialVoltage    = U;
            css.dischargeFunction = output.dischargeFunction;
        
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

            lines = addLine(lines, 'Packing mass', 'kg', css.packingMass);
            lines = addLine(lines, 'Temperature' , 'C' , css.temperature -  273.15);
            lines = addLine(lines, 'Mass'        , 'kg', css.mass);
            lines = addLine(lines, 'Volume'      , 'L' , css.volume/litre);
            if ~isempty(css.thicknesses)
                lines = addLine(lines, 'Negative Electrode Coating Thickness', 'µm', css.thicknesses.(ne)/(micro*meter));
                lines = addLine(lines, 'Positive Electrode Coating Thickness', 'µm', css.thicknesses.(pe)/(micro*meter));
                lines = addLine(lines, 'Negative Electrode Mass Loading', 'mg/cm^2', css.massLoadings.(ne)/(milli*gram/((centi*meter)^2)));
                lines = addLine(lines, 'Positive Electrode Mass Loading', 'mg/cm^2', css.massLoadings.(pe)/(milli*gram/((centi*meter)^2)));
            end
            lines = addLine(lines, 'Total Capacity'             , 'Ah'   , css.capacity/hour);
            lines = addLine(lines, 'Negative Electrode Capacity', 'Ah'   , css.capacities.(ne)/hour);
            lines = addLine(lines, 'Positive Electrode Capacity', 'Ah'   , css.capacities.(pe)/hour);
            lines = addLine(lines, 'N/P ratio'                  , '-'    , css.NPratio);
            lines = addLine(lines, 'Energy'                     , 'Wh'   , css.energy/hour);
            lines = addLine(lines, 'Specific Energy'            , 'Wh/kg', css.specificEnergy/hour);
            lines = addLine(lines, 'Energy Density'             , 'Wh/L' , (css.energyDensity/hour)*litre);

            function str = appendDrate(str, crate)

                str = sprintf('%s (DRate = %g)', str, crate);
                
            end
            
            for isim = 1 : numel(css.dischargeSimulations)
                
                simres = css.dischargeSimulations{isim};
                ac = @(str) appendDrate(str, simres.DRate); 
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

                fmt = sprintf('%%%ds : %%-8g %%s\n', s);

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

        function css = addDrates(css, DRates, varargin)

            for icrate = 1 : numel(DRates)

                DRate = DRates(icrate);
                css = css.addDrate(DRate, varargin{:});

            end
            
        end
        
        function css = addDrate(css, DRate, varargin)

            % the extras options are passed to computeCellEnergy
            extras = varargin;
            
            model = css.model;
            
            [energy, output] = computeCellEnergy(model, 'DRate', DRate, extras{:});
            
            specificEnergy = energy/css.mass;
            energyDensity  = energy/css.volume;

            simresult = struct('DRate'            , DRate                   , ...
                               'energy'           , energy                  , ...
                               'specificEnergy'   , specificEnergy          , ...
                               'energyDensity'    , energyDensity           , ...
                               'dischargeFunction', output.dischargeFunction, ...
                               'E'                , output.E                , ...
                               'I'                , output.I                , ...
                               'time'             , output.time);

            css.dischargeSimulations{end + 1} = simresult;
            
        end
            
        
    end

    
end



%{
Copyright 2021-2024 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BattMo

BattMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BattMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BattMo.  If not, see <http://www.gnu.org/licenses/>.
%}
