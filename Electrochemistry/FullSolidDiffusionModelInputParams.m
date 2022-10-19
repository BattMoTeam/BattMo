classdef FullSolidDiffusionModelInputParams < SolidDiffusionModelInputParams

    properties

        N  % Number of discretization intervals in the diffusion model [-]
        np % Number of computational grid cells (typically set by parent model :class:`ActiveMaterial <Electrochemistry.ActiveMaterial>`)

        volumeFraction
        activeMaterialFraction = 1
        
        D    % Function to update D value given as a struct with fields
             % D.type is in {'function', 'constant'}. If 'constant' is chosen the value of D0 is used as in parent class
             % D.functionname :  matlab function name (should be available in path)
             % D.argumentlist = ["c", "cmin", "cmax"]

        % needed if function is used
        cmax     % maximum concentration [mol/m^3]
        theta0   % Minimum lithiation, 0% SOC    [-]
        theta100 % Maximum lithiation, 100% SOC  [-]

    end
    
    methods
        
        function paramobj = FullSolidDiffusionModelInputParams(jsonstruct)
            paramobj = paramobj@SolidDiffusionModelInputParams(jsonstruct);
        end

        function paramobj = validateInputParams(paramobj)

            D0 = paramobj.D0;
            D  = paramobj.D;
            
            assert(~isempty(D0) || ~isempty(D), 'Either D0 or D should be provided');

            if ~isempty(D) & strcmp(D.type, 'constant')
                assert(~isempty(D0), 'D0 should be provided when D is set to be constant')
            end
            
        end
        
    end
    
    
end



%{
Copyright 2021-2022 SINTEF Industry, Sustainable Energy Technology
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
