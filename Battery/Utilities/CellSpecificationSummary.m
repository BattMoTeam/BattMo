classdef CellSpecificationSummary
% Utility class to compute standard cell specifications (see list in properties below) using the model
%
% Energy computation for given DRate can be added using the method addDrate. The results will be stored in the property
% dischargeSimulations.

    properties (SetAccess = immutable)

        % We want the computed values to remain synchronized with the model. To do so we set the "setAccess" property to
        % immutable.

        model

        gridGenerator

        jsonstruct

    end

    properties (SetAccess = private)

        has_packing

        packing_mass
        packing_volume

        thicknesses % Used to compute mass loading

        mass   % total mass (includes packing mass if packing is present)
        masses % structure with mass of each of the cell components

        massLoadings % Computed if thicknesses are available

        volume  % total volume (includes packing volume if packing is present)
        volumes % structure with volume of each of the cell components

        energy
        specificEnergy
        energyDensity

        NPratio

        capacity
        capacities % structure with capacity of negative and positive electrode

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

        function css = CellSpecificationSummary(S, varargin)

            % S may be either a model or a jsonstruct

            opt = struct('packing_mass'  , []   , ...
                         'packing_volume', []   , ...
                         'total_volume'  , []   , ...
                         'thicknesses'   , []   , ...
                         'jsonstruct'    , []   , ...
                         'gridGenerator' , []   , ...
                         'temperature'   , 298);
            opt = merge_options(opt, varargin{:});

            if isstruct(S)
                assert(isempty(opt.jsonstruct), 'If input is jsonstruct, the opt.jsonstruct should not be given');
                opt.jsonstruct = S;
                css.jsonstruct = S;
                [css.model, ~, ~, opt.gridGenerator] = setupModelFromJson(S);
            else
                css.model = S;
            end

            css.packing_mass         = opt.packing_mass;
            css.temperature          = opt.temperature;
            css.gridGenerator        = opt.gridGenerator;
            css.dischargeSimulations = {};

            thicknesses = css.extractThicknessFromModel(); % can be done in 1D model or when gridGenerator is given (if possible)

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

        function css = setupPacking(css, varargin)

            opt = struct('packing_mass'  , []  , ...
                         'packing_volume', []  , ...
                         'total_volume'  , []);

            opt = merge_options(opt, varargin{:});

            packing_mass   = opt.packing_mass;
            packing_volume = opt.packing_volume;
            total_volume   = opt.total_volume;

            if isempty(packing_volume) & isempty(total_volume) & isempty(packing_mass)
                css.has_packing = false;
                return
            end

            css.has_packing = true;

            if isempty(packing_mass)
                fprintf('Packing mass is not given. We set the packing mass to zero.\n');
                packing_mass = 0;
            end
            css.packing_mass = packing_mass;

            if ~isempty(total_volume)

                css.volume = total_volume;
                assert(css.volume > css.volumes.val, 'Total volume given is smaller that the sum of the components in the cell');
                if ~isempty(packing_volume)
                    fprintf('Both total volume and packing volume are given. We use the total volume\n')
                end
                if ~isempty(css.volumes)
                    packing_volume = css.volume - css.volumes.val;
                else
                    packing_volume = [];
                end
                css.packing_volume = packing_volume;

            else

                if ~isempty(packing_volume)
                    css.packing_volume = packing_volume;
                else
                    fprintf('Packing volume or total volume are not given. We set the packing volume to zero\n');
                    css.packing_volume = 0;
                end

                if ~isempty(css.volumes)
                    total_volume = css.packing_volume + css.volumes.val;
                else
                    total_volume = [];
                end

                css.volume = total_volume;

            end

            css = css.computeSpecs();

        end

        function thicknesses = extractThicknessFromModel(css)

            gridgen = css.gridGenerator;
            model   = css.model;

            if isempty(gridgen) & model.grid.griddim > 1
                thicknesses = [];
                return
            end

            ne = 'NegativeElectrode';
            pe = 'PositiveElectrode';
            co = 'Coating';

            if model.grid.griddim == 1

                % we recover the lengths from the model directly

                eldes = {ne, pe};

                for ielde = 1 : numel(eldes)

                    elde = eldes{ielde};

                    G = model.(elde).(co).grid;
                    xmax = max(G.faces.centroids(:, 1));
                    xmin = min(G.faces.centroids(:, 1));

                    thicknesses.(elde) = xmax - xmin;

                end

                return
            end

            switch class(gridgen)

              case {'BatteryGeneratorMultilayerPouch'}

                thicknesses.(ne) = gridgen.unit_cell_thickness(2);
                thicknesses.(pe) = gridgen.unit_cell_thickness(4);

              case {'BatteryGeneratorP4D'}

                thicknesses.(ne) = gridgen.zlength(2);
                thicknesses.(pe) = gridgen.zlength(4);

              otherwise

                error('grid generator class not recognized');

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

            [mass, masses, volumes] = computeCellMass(model);

            if css.has_packing
                mass = mass + css.packing_mass;
            end

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

            if css.has_packing
                if isempty(css.volume)
                    volume = ccs.packing_volume + volumes.val;
                else
                    volume = css.volume;
                end
            else
                volume = volumes.val;
            end

            specificEnergy = energy/mass;
            energyDensity  = energy/volume;

            % Compute NP ratio

            NPratio = capacities.(ne)/capacities.(pe);

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
            css.dischargeFunction = output.dischargeFunction;

        end

        function output = printSpecifications(css, varargin)

            opt = struct('returnStruct', false);
            opt = merge_options(opt, varargin{:});

            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';

            function lines = addLine(lines, description, unit, value)

                line.description = description;
                line.unit        = unit;
                line.value       = value;

                lines{end + 1} = line;

            end

            lines = {};

            if css.has_packing
                lines = addLine(lines, 'Packing mass', 'kg', css.packing_mass);
                lines = addLine(lines, 'Packing volume', 'L', css.packing_volume/litre);
            else
                lines = addLine(lines, 'No packing', [], []);
            end
            lines = addLine(lines, 'Temperature' , 'C' , css.temperature -  273.15);
            if css.has_packing
                lines = addLine(lines, 'Mass (including packing)'  , 'kg', css.mass);
                lines = addLine(lines, 'Volume (including packing)', 'L' , css.volume/litre);
            else
                lines = addLine(lines, 'Mass'        , 'kg', css.mass);
                lines = addLine(lines, 'Volume'      , 'L' , css.volume/litre);
            end
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
            lines = addLine(lines, 'Nominal Energy'             , 'Wh'   , css.energy/hour);
            lines = addLine(lines, 'Nominal Specific Energy'    , 'Wh/kg', css.specificEnergy/hour);
            lines = addLine(lines, 'Nominal Energy Density'     , 'Wh/L' , (css.energyDensity/hour)*litre);

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
                fmt2 = sprintf('%%%ds\n', s);

                for iline = 1 : numel(lines)
                    line = lines{iline};
                    if ~isempty(line.value)
                        fprintf(fmt             , ...
                                line.description, ...
                                line.value      , ...
                                line.unit);
                    else
                        fprintf(fmt2            , ...
                                line.description);
                    end
                end

            end

            printLines(lines);

            if opt.returnStruct

                s = [lines{:}];
                output = struct();

                for i = 1:numel(s)
                    fieldName = matlab.lang.makeValidName(s(i).description);
                    output.(fieldName) = struct( ...
                        'Unit',  s(i).unit, ...
                        'Value', s(i).value);

                end

            else

                output = lines;

            end

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
