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


            inputparams = inputparams@ActiveMaterialInputParams(jsonstruct);
            
            pick = @(fd) pickField(jsonstruct, fd);

            inputparams.SideReaction            = SideReactionInputParams(pick('SideReaction'));
            inputparams.SolidElectrodeInterface = SolidElectrodeInterfaceInputParams(pick('SolidElectrodeInterface'));

        end

        function inputparams = setupSEImodel(inputparams, jsonstruct)

            if isfield(jsonstruct, 'SEImodel') && ~strcmp(jsonstruct.SEImodel, 'Safari')
                error('Inconsitent model, check SEI model values')
            end

            inputparams.SEImodel = 'Safari';
            
        end
                
        function inputparams = setupSolidDiffusion(inputparams, jsonstruct)

            sd = 'SolidDiffusion';
            input.diffusionModelType = inputparams.diffusionModelType;
            inputparams.(sd) = FullSolidDiffusionModelInputParams(pickField(jsonstruct, sd));

        end
        
        function inputparams = validateInputParams(inputparams)

            inputparams = validateInputParams@ActiveMaterialInputParams(inputparams);

            assert(strcmp(inputparams.SEImodel, 'Safari'), 'The SEI model should be set to Safari');
            
            if inputparams.isRootSimulationModel
                % only one particle in the stand-alone model
                inputparams.SolidElectrodeInterface.np = 1;
            end

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
