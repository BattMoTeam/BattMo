classdef TestJsonFiles < matlab.unittest.TestCase

    properties (TestParameter)

        jsonLintFile = TestJsonFiles.findJsonFiles('*.json');

        jsonSchemaFile = TestJsonFiles.findJsonFiles('*.schema.json');

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
        resolveModule = 'resolveFileInputJson';
        validateModule = 'validateJsonFiles';
        isPySetup;

    end

    methods

        function test = TestJsonFiles()

            try
                setupPythonExecutable();
                setupPythonPath();

                modulenames = {test.lintModule, test.resolveModule, test.validateModule};
                loadModule(modulenames, 'setupPython', false);
                for imod = 1 : numel(modulenames)
                    pymodule = py.importlib.import_module(modulenames{imod});
                    py.importlib.reload(pymodule);
                end
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

            assert(ok, 'JSON linting failed for %s', jsonLintFile);

        end

        function testJsonSchema(test, jsonSchemaFile)

            ok = false;

            if test.isPySetup
                dispif(mrstVerbose, 'Validating JSON schema %s\n', jsonSchemaFile);
                ok = py.(test.validateModule).checkSchema(jsonSchemaFile);
            end

            assert(ok, 'JSON schema validation failed for %s', jsonSchemaFile);

        end

        function testJsonDataSet(test, jsonDataSet)

            ok = false;

            if test.isPySetup
                dispif(mrstVerbose, 'Validating %s\n', jsonDataSet);
                test.assumeFalse(contains(jsonDataSet, test.excludeJsonDataSet));
                ok = py.(test.validateModule).validate(battmoDir(), jsonDataSet);
            end

            assert(ok, 'JSON dataset validation failed for %s', jsonDataSet);

        end

        function testJsonDataFile(test, jsonDataFile)

            ok = false;
            jsonfile = fullfile(battmoDir(), 'Examples', 'JsonDataFiles', jsonDataFile{1});
            schemafile = [jsonDataFile{2}, '.schema.json'];

            if test.isPySetup
                dispif(mrstVerbose, 'Validating %s against %s\n', jsonfile, schemafile);
                test.assumeFalse(contains(jsonfile, test.excludeJsonDataFile));
                ok = py.(test.validateModule).validate(battmoDir(), jsonfile, schemafile);
            end

            assert(ok, 'JSON validation failed for %s against %s', jsonfile, schemafile);

        end

    end

    methods (Static)

        function filenames = findJsonFiles(pattern)

            rootdir = getCanonicalPath(battmoDir());
            files = dir(fullfile(rootdir, '**', pattern));
            folders = {files.folder};

            excludedRoots = {fullfile(rootdir, 'Externals'), ...
                             fullfile(rootdir, '.git'), ...
                             fullfile(rootdir, 'Documentation', '.venv')};

            include = true(1, numel(files));
            for iroot = 1 : numel(excludedRoots)
                root = [excludedRoots{iroot}, filesep];
                include = include & ~startsWith(folders, root);
            end

            files = files(include);
            filenames = arrayfun(@(s) fullfile(s.folder, s.name), files, 'uniformoutput', false);

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
