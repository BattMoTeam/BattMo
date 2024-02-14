function [E, energyDensity, energy] = runJsonFunction(jsonfiles, varargin)
%
%
% SYNOPSIS:
%   function [E, energyDensity, energy] = runJsonFunction(jsonfiles, varargin)
%
% DESCRIPTION: Run a battery simulation as specified by one or several json files.
%
% PARAMETERS:
%   jsonfiles - Specification of material and geometry as one or several json files.
%               Either:
%                1. A single json file as a character array (single quotes)
%                2. A cell array ({}) with one or several json files.
%
% OPTIONAL PARAMETERS:
%   do_plot - True/false (default is false).
%
% RETURNS:
%   E             -
%   energyDensity -
%   energy        -
%
% EXAMPLE:
%
%   [E, energyDensity, energy] = runJsonFunction('a_json_file.json', 'do_plot', true);
%
%   E = runJsonFunction({'lithium_ion_battery_nmc_graphite.json', 'geometry1d.json});
%
% SEE ALSO:
%

    mrstModule add ad-core mrst-gui mpfa agmg linearsolvers

    opt = struct('do_plot', false);
    opt = merge_options(opt, varargin{:});

    if ~exist('jsonfiles', 'var')
        jsonfiles = {
            fullfile('ParameterData', 'BatteryCellParameters', 'LithiumIonBatteryCell',...
                     'lithium_ion_battery_nmc_graphite.json'),...
            fullfile('Examples', 'JsonDataFiles', 'geometry1d.json'), ...
            fullfile('Examples', 'JsonDataFiles', 'extra_output.json'), ...
                    };
    end

    % Parse
    if ischar(jsonfiles)

        % Single json file
        jsonstruct = parseBattmoJson(jsonfiles);

    elseif iscell(jsonfiles)

        jsonstructs = cell(numel(jsonfiles), 1);
        for k = 1:numel(jsonfiles)
            jsonstructs{k} = parseBattmoJson(jsonfiles{k});
        end

        % Merge
        jsonstruct = mergeJsonStructs(jsonstructs);

    else
        error('Cannot parse input variable jsonfiles. It should either be a character array ('''') or a cell array of one or more character arrays {'''', '''', ...}.');
    end

    % Run battery simulation with function that takes json input
    % NB: the output is a struct with many fields.
    output = runBatteryJson(jsonstruct);

    %%

    E = output.E;
    energyDensity = output.energyDensity;
    energy = output.energy;

    if opt.do_plot
        figure
        plot(energyDensity, E)
        xlabel('Energy Density [Wh/L]');
        ylabel('Voltage [V]');

        figure
        plot(energy, E)
        xlabel('Energy [Wh]');
        ylabel('Voltage [V]');
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
