classdef TestJsonFiles < matlab.unittest.TestCase

    properties (TestParameter)

        jsonfile = {fullfile('ParameterData', 'BatteryCellParameters', 'LithiumIonBatteryCell', 'lithium_ion_battery_nmc_graphite.json'),
                    fullfile('ParameterData', 'ParameterSets', 'Xu2015', 'lfp.json'),
                    fullfile('ParameterData', 'ParameterSets', 'Chen2020', 'chen2020_lithium_ion_battery.json')};
    end
    
    methods (Test)

        function testJson(test, jsonfile)

            loadModule('validationJsonScript')
            is_valid = py.validationJsonScript.validate(jsonfile);
            assert(is_valid);
            
        end

    end

end

%{
Copyright 2009-2022 SINTEF Digital, Mathematics & Cybernetics.

This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).

MRST is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MRST is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MRST.  If not, see <http://www.gnu.org/licenses/>.
%}
