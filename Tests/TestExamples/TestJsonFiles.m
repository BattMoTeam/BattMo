classdef TestJsonFiles < matlab.unittest.TestCase

    properties (TestParameter)

        jsonLintFile = arrayfun(@(s) fullfile(s.folder, s.name), dir(fullfile(battmoDir(), '**', '*.json')), 'uniformoutput', false);

        jsonDataSet = {
            fullfile('ParameterData', 'BatteryCellParameters', 'LithiumIonBatteryCell'      , 'lithium_ion_battery_nmc_graphite.json'), ...
            fullfile('ParameterData', 'ParameterSets'        , 'Chen2020'                   , 'chen2020_lithium_ion_battery.json')    , ...
            fullfile('ParameterData', 'ParameterSets'        , 'Xu2015'                     , 'lfp.json')                             , ...
            fullfile('Examples'     , 'JsonDataFiles'        , 'p2d_40_jl.json')            , ...
            fullfile('Examples'     , 'JsonDataFiles'        , 'p2d_40_jl_ud.json')         , ...
            fullfile('Examples'     , 'JsonDataFiles'        , 'p2d_40.json')               , ...
            fullfile('Examples'     , 'JsonDataFiles'        , 'sample_input.json')         , ...
                      };

        jsonDataFile = {
            {'4680-geometry.json'          , 'Geometry'}    , ...
            {'geometry1d.json'             , 'Geometry'}    , ...
            {'geometry3d.json'             , 'Geometry'}    , ...
            {'geometryChen.json'           , 'Geometry'}    , ...
            {'geometryMultiLayerPouch.json', 'Geometry'}    , ...
            {'cccv_control.json'           , 'ControlModel'}, ...
            {'cc_discharge_control.json'   , 'ControlModel'}, ...
            {'extra_output.json'           , 'Output'}      , ...
            {'linear_solver_setup.json'    , 'Solver'}      , ...
            {'silicongraphite.json'        , 'Coating'}     , ...
            {'simulation_parameters.json'  , 'TimeStepping'}, ...
                       };

    end

    properties

        excludeJsonLintFile = {};
        excludeJsonDataSet  = {'p2d_40_jl_ud'};
        excludeJsonDataFile = {};

        lintModule = 'checkLint';
        validateModule = 'validateJsonFiles';
        isPySetup;

    end

    methods

        function test = TestJsonFiles()

            try
                setupPythonExecutable();
                setupPythonPath();

                modulenames = {test.lintModule, test.validateModule};
                loadModule(modulenames, 'setupPython', false);
                test.isPySetup = true;

            catch
                test.isPySetup = false;
            end

        end

    end

    methods (Test)

        function testJsonLint(test, jsonLintFile)

            ok = false;

            if test.isPySetup
                dispif(mrstVerbose, 'Linting %s\n', jsonLintFile);
                test.assumeFalse(contains(jsonLintFile, test.excludeJsonLintFile));
                ok = py.(test.lintModule).check(jsonLintFile);
            end

            assert(ok);

        end

        function testJsonDataSet(test, jsonDataSet)

            ok = false;

            if test.isPySetup
                dispif(mrstVerbose, 'Validating %s\n', jsonDataSet);
                test.assumeFalse(contains(jsonDataSet, test.excludeJsonDataSet));
                %ok = py.(test.validateModule).validate(battmoDir(), jsonDataSet);
                ok = validateJsonFiles(jsonDataSet);
            end

            assert(ok);

        end

        function testJsonDataFile(test, jsonDataFile)

            ok = false;

            if test.isPySetup
                jsonfile = fullfile(battmoDir(), 'Examples', 'JsonDataFiles', jsonDataFile{1});
                schemafile = [jsonDataFile{2}, '.schema.json'];

                dispif(mrstVerbose, 'Validating %s against %s\n', jsonfile, schemafile);
                test.assumeFalse(contains(jsonfile, test.excludeJsonDataFile));
                ok = py.(test.validateModule).validate(battmoDir(), jsonfile, schemafile);
            end

            assert(ok);

        end

        function testJsonStructs(test, jsonDataFile)

            function pyDict = struct_to_pydict(matStruct)
                if isstruct(matStruct)
                    pyDict = py.dict();
                    fields = fieldnames(matStruct);
                    for i = 1:numel(fields)
                        field = fields{i};
                        value = matStruct.(field);
                        pyDict{field} = struct_to_pydict(value); % Recursively convert nested structs
                    end
                else
                    pyDict = matStruct; % Base case: non-struct values
                end
            end

            ok = false;

            if test.isPySetup
                jsonfile = fullfile(battmoDir(), 'Examples', 'JsonDataFiles', jsonDataFile{1});
                jsonstruct = parseBattMoJson(jsonfile);
                pyDict = struct_to_pydict(jsonstruct);
                schemafile = [jsonDataFile{2}, '.schema.json'];

                dispif(mrstVerbose, 'Validating %s against %s\n', jsonfile, schemafile);
                test.assumeFalse(contains(jsonfile, test.excludeJsonDataFile));
                ok = py.(test.validateModule).validateStructs(battmoDir(), pyDict, schemafile);
            end

            assert(ok);

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
