classdef TestJsonFiles < matlab.unittest.TestCase

    properties (TestParameter)

        jsonfile = {fullfile('ParameterData', 'BatteryCellParameters', 'LithiumIonBatteryCell', 'lithium_ion_battery_nmc_graphite.json'), ...
                    fullfile('ParameterData', 'ParameterSets', 'Chen2020', 'chen2020_lithium_ion_battery.json'), ...
                    fullfile('ParameterData', 'ParameterSets', 'Lin2015', 'lfp.json')};

        jsonlintfile = arrayfun(@(s) fullfile(s.folder, s.name), dir(fullfile(battmoDir(), '**', '*.json')), 'uniformoutput', false);

    end

    properties

        exclude = { };

    end

    methods (Test)

        function testJsonLint(test, jsonlintfile)

            if ~contains(jsonlintfile, test.exclude)
                % fprintf('Linting %s\n', jsonlintfile);
                loadModule('checkLint');
                is_valid = py.checkLint.check(jsonlintfile);
                assert(is_valid, jsonlintfile);
            end

        end

        function testJson(test, jsonfile)

            if ~contains(jsonfile, test.exclude)
                % fprintf('Validating %s\n', jsonfile);
                loadModule('validationJsonScript')
                is_valid = py.validationJsonScript.validate(jsonfile);
                assert(is_valid);
            end

        end

    end

end

%{
Copyright 2021-2023 SINTEF Industry, Sustainable Energy Technology
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
