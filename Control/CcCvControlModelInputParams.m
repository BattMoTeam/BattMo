classdef CcCvControlModelInputParams < ControlModelInputParams
%
% Constant-Current-Constant-Voltage Control. We switch between charge and discharge scenarii. The switch occurs when the
% derivative of the non-controlled variable gets below a given threshold
%
    
    properties
        %
        % When voltage control, we wait for the derivative of the current to be below the value of dIdtLimit between we switch to the following constant current control.
        %
        dIdtLimit
        
        %
        % Not used for the moment
        %
        dEdtLimit         
        
    end
    
    methods

        function paramobj = CcCvControlModelInputParams(jsonstruct);
            
            paramobj = paramobj@ControlModelInputParams(jsonstruct);
            paramobj.controlPolicy = 'CCCV';
            
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
