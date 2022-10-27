classdef SEIActiveMaterialInputParams < ActiveMaterialInputParams
%
% Input parameter class for :class:`ActiveMaterial <Electrochemistry.ActiveMaterial>`
% 
    properties

        SolidElectrodeInterface
        SideReaction        
    end

    methods

        function paramobj = SEIActiveMaterialInputParams(jsonstruct)

            paramobj = paramobj@ActiveMaterialInputParams(jsonstruct);

            pick = @(fd) pickField(jsonstruct, fd);
            % For SEI, we always use full diffusion model
            paramobj.SolidDiffusion          = FullSolidDiffusionModelInputParams(pick('SolidDiffusion'));
            paramobj.SideReaction            = SideReactionInputParams(pick('SideReaction'));
            paramobj.SolidElectrodeInterface = SolidElectrodeInterfaceInputParams(pick('SolidElectrodeInterface'));

            paramobj = paramobj.validateInputParams();
            
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
