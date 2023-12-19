classdef TestRunExamples < matlab.unittest.TestCase

    properties (TestParameter)

        % List of all the supported scripts

        filename = {'battMoTutorial'              , ...
                    'runBatteryMLP'               , ...
                    'runBatteryP2D'               , ...
                    'runBatteryP3D'               , ...
                    'runBatteryP4D'               , ...
                    'runBatteryPrecondTestP2D.m'  , ...
                    'runChen2020.m'               , ...
                    'runCR.m'                     , ...
                    'runGittTest.m'               , ...
                    'runJellyRoll.m'              , ...
                    'runSector.m'                 , ...
                    'runActiveMaterial.m'         , ...
                    'runSEIActiveMaterial.m'      , ...
                    'runSiliconGraphiteBattery'   , ...
                    'runBatteryLinearSolve.m'     , ...
                    'runJsonFunction.m'           , ...
                    'runJsonScript.m'             , ...
                    'runSetupModel.m'             , ...
                    'runBatteryP3DMech.m'         , ...
                    'runBattery1DOptimize.m'      , ...
                    'runParameterIdentification.m', ...
                    'runElectrolyser.m'           , ...
                    'magnesium_script'};

    end

    properties

        % Scripts that are not yet fully supported 
        
        exclude = {'runJellyRollLinearSolve.m', ...
                   'runSingleParticleSEI.m'   , ...
                   'runParameterSweep'        , ...
                   'runDissElectrolyser.m'};
        
    end
    
    methods (Test)

        function testRunExample(test, filename)

            run(fullfile(battmoDir, 'startupBattMo.m'));

            if ~contains(filename, exclude)
                % FIXME Disable plotting
                set(0, 'defaultFigureVisible', 'off');
                fprintf('\n\nRunning %s...\n\n', filename);
                run(filename);
                close all
                set(0, 'defaultFigureVisible', 'on');
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
