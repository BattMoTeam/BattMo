classdef FullSolidDiffusionModelInputParams < SolidDiffusionModelInputParams
%
% Full diffusion model (standard PXD)
%
    
    properties


        %% Standard Parameters


        % Function to update diffusion coefficient value, given as a struct with fields
        % - type         : element in {'function', 'constant'}. If 'constant' is chosen the value of referenceDiffusionCoefficient defined in parent class is used
        % - functionName : matlab function name (should be available in path)
        % - argumentList : should be  ["c", "cmin", "cmax"]
        diffusionCoefficient    

        % The following parameters are only need in the case the diffusionCoefficient is given as a function (see argument list)
        saturationConcentration % the saturation concentration of the guest molecule in the host material (symbol: cmax)
        guestStoichiometry100   % the ratio of the concentration of the guest molecule to the saturation concentration
                                % of the guest molecule in a phase at a cell voltage that is defined as 100% SOC(symbol: theta100)
        guestStoichiometry0     % the ratio of the concentration of the guest molecule to the saturation concentration
                                % of the guest molecule in a phase at a cell voltage that is defined as 0% SOC (symbol: theta0)


        %% Discretization parameters
        
        % Number of discretization intervals in the diffusion model [-]
        N
        % Number of computational grid cells (typically set by parent model :class:`ActiveMaterial <Electrochemistry.ActiveMaterialInputParams>`)
        np 

    end
    
    methods
        
        function inputparams = FullSolidDiffusionModelInputParams(jsonstruct)

            D0 = getJsonStructField(jsonstruct, 'referenceDiffusionCoefficient');
            D  = getJsonStructField(jsonstruct, 'diffusionCoefficient');
            
            assert(isAssigned(D0) || isAssigned(D), 'Either D0 or D should be provided');

            inputparams = inputparams@SolidDiffusionModelInputParams(jsonstruct);
            
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
