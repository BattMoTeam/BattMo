classdef FullSolidDiffusionModelInputParams < SolidDiffusionModelInputParams
%
% Full diffusion model (standard PXD)
%

    
    properties
        %
        % Number of discretization intervals in the diffusion model [-]
        %
        N

        %
        %  Number of computational grid cells (typically set by parent model :class:`ActiveMaterial <Electrochemistry.ActiveMaterialInputParams>`)
        %
        np 

        %
        % Volume Fraction (typically set by parent model :class:`ActiveMaterial <Electrochemistry.ActiveMaterialInputParams>`)
        %
        volumeFraction

        %
        %  Active Material Fraction (typically set by parent model :class:`ActiveMaterial <Electrochemistry.ActiveMaterialInputParams>`)
        %
        activeMaterialFraction

        %
        % Function to update D value given as a struct with fields
        % - D.type is in {'function', 'constant'}. If 'constant' is chosen the value of D0 defined in parent class
        % - D.functionname :  matlab function name (should be available in path)
        % - D.argumentlist = ["c", "cmin", "cmax"]
        % 
        D    

        %
        % maximum concentration [mol/m^3] (only needed if D is a function)
        %
        cmax
        %
        % Minimum lithiation, 0% SOC [-] (only needed if D is a function)
        %
        theta0
        %
        % Maximum lithiation, 100% SOC [-] (only needed if D is a function)
        %
        theta100 

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
                paramobj.D0 = paramobj.D.value;
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
