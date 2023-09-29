classdef CoatingInputParams < ElectronicComponentInputParams
%
% Input parameter class for :code:`Coating` model
% 
    properties

        %% Sub-Models
        
        ActiveMaterial
        Binder
        ConductingAdditive
        
        %% Standard parameters

        density              % the mass density of the material (symbol: rho)
        bruggemanCoefficient % the Bruggeman coefficient for effective transport in porous media (symbol: beta)
        
        %% Coupling parameters
        
        externalCouplingTerm % structure to describe external coupling (used in absence of current collector)

    end

    methods

        function paramobj = CoatingInputParams(jsonstruct)

            paramobj = paramobj@ElectronicComponentInputParams(jsonstruct);
            
            pick = @(fd) pickField(jsonstruct, fd);

            paramobj.ActiveMaterial     = ActiveMaterialInputParams(pick('ActiveMaterial'));
            paramobj.Binder             = BinderInputParams(pick('Binder'));
            paramobj.ConductingAdditive = ConductingAdditiveInputParams(pick('ConductingAdditive'));
            
            paramobj = paramobj.validateInputParams();
            
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
