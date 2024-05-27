classdef SEIActiveMaterialInputParams < ActiveMaterialInputParams
%
% Input parameter class for :code:`SEIActiveMaterial` model
% 
    properties

        SolidElectrodeInterface
        SideReaction
        
    end

    methods

        function inputparams = SEIActiveMaterialInputParams(jsonstruct)
            
            errorMessage = 'We should have the Safari SEI model as input (otherwise this model should not be called)';
            jsonstruct = setJsonStructField(jsonstruct, 'SEImodel', 'Safari', 'errorMessage', errorMessage);
            
            errorMessage = 'For the Safari SEI model, the use of the full diffusion model is required';
            jsonstruct = setJsonStructField(jsonstruct, 'diffusionModelType', 'full', 'errorMessage', errorMessage);
            
            isRootSimulationModel = getJsonStructField(jsonstruct, 'isRootSimulationModel');
            
            if isAssigned(isRootSimulationModel) && isRootSimulationModel
                % only one particle in the stand-alone model
                jsonstruct = setJsonStructField(jsonstruct, {'SolidElectrodeInterface', 'np'}, 1, 'handleMisMatch', 'quiet');
            end

            inputparams = inputparams@ActiveMaterialInputParams(jsonstruct);
            
            pick = @(fd) pickField(jsonstruct, fd);

            inputparams.SideReaction            = SideReactionInputParams(pick('SideReaction'));
            inputparams.SolidElectrodeInterface = SolidElectrodeInterfaceInputParams(pick('SolidElectrodeInterface'));

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
