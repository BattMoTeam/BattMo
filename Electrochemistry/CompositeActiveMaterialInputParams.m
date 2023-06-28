classdef CompositeActiveMaterialInputParams < ElectronicComponentInputParams
%
% Input parameter class for :class:`ActiveMaterial <Electrochemistry.ActiveMaterial>`
% 
    properties

        FirstMaterial

        SecondMaterial

        externalCouplingTerm % structure to describe external coupling (used in absence of current collector)
        
        volumeFraction % Volume fraction of the whole material (binder and so on included)

        use_particle_diffusion

    end

    methods

        function paramobj = CompositeActiveMaterialInputParams(jsonstruct)

            paramobj = paramobj@ElectronicComponentInputParams(jsonstruct);

            pick = @(fd) pickField(jsonstruct, fd);

            paramobj.FirstMaterial = ActiveMaterialInputParams(pick('FirstMaterial'));
            paramobj.SecondMaterial = ActiveMaterialInputParams(pick('SecondMaterial'));

            paramobj = paramobj.validateInputParams();
            
        end

        function paramobj = validateInputParams(paramobj)

            gr = 'FirstMaterial';
            si = 'SecondMaterial';


            %paramobj = mergeParameters(paramobj,  {{gr, 'volumeFraction'}, {'volumeFraction'}});
            %paramobj = mergeParameters(paramobj,  {{si, 'volumeFraction'}, {'volumeFraction'}});

            paramobj = validateInputParams@ElectronicComponentInputParams(paramobj);
            
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
