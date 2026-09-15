clear all
close all

printer = @(s) disp(jsonencode(s, 'PrettyPrint', true));

experiment = {
    % 'Discharge at 1C for 0.5 hours';
    % 'Discharge at C/20 for 0.5 hours';
    % 'Charge at 0.5 C for 45 minutes';
    'Discharge at 1 A for 0.5 hours';
    'Charge at 200 mA for 45 minutes';
    %'Discharge at 1W for 0.5 hours';
    %'Charge at 200mW for 45 minutes';
    'Rest for 10 minutes';
    'Hold at 1V for 20 seconds';
    %'Charge at 1 C until 4.1V';
    'Hold at 4.1 V until 50mA';
    % 'Hold at 3V until C/50';
    % 'Discharge at C/3 for 2 hours or until 2.5 V';
             };

for itest = 1:numel(experiment)
    jsonstruct = convertPyBaMMtoJson(experiment{itest});
    disp(experiment{itest});
    printer(jsonstruct);

    schema = fullfile(battmoDir, 'Utilities', 'JsonSchemas', 'GenericControl.schema.json');
    validateStruct(jsonstruct, 'useTmpFile', false, 'schema', schema);

end

%{
Copyright 2021-2026 SINTEF Industry, Sustainable Energy Technology
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
