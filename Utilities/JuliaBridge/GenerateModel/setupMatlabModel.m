function export = setupMatlabModel(casename, jsonfolder, generate_reference_solution,varargin)
%% Script for BattMo.m to produce reference solution
% - casename : is used to identify json file name for the input and saved data output
% - jsonfolder : folder where the json file is fetched
% - datafolder : folder where the computed data is saved
    
    % load json setup file

    json_filename = sprintf('%s.json', casename);
    json_filename = fullfile(jsonfolder, json_filename);

    jsonstruct = parseBattmoJson(json_filename);
    
    %% To run the simulation, you need to install matlab battmo

    export = runBatteryJson(jsonstruct, 'runSimulation', generate_reference_solution);

    %% Save solution as a matlab struct that can be imported in Julia

    export.model    = class2data(export.model);
    export.schedule = class2data(export.schedule);
    
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
