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

        density              % the mass density of the material (symbol: rho). Important : the density is computed with respect to total volume (including the empty pores)
        bruggemanCoefficient % the Bruggeman coefficient for effective transport in porous media (symbol: beta)

        %% Advanced parameters

        volumeFractions
        volumeFraction
        thermalConductivity             % (if not given computed from the subcomponents)
        specificHeatCapacity            % (if not given computed from the subcomponents)
        effectiveThermalConductivity    % (account for volume fraction)
        effectiveVolumetricHeatCapacity % (account for volume fraction and density)
        
        %% Coupling parameters
        
        externalCouplingTerm % structure to describe external coupling (used in absence of current collector)

        %% particle case
        
        activematerial_type % 'default' (only one particle type) or 'composite' (two different particles)
        
    end

    methods

        function paramobj = CoatingInputParams(jsonstruct)

            paramobj = paramobj@ElectronicComponentInputParams(jsonstruct);
            
            pick = @(fd) pickField(jsonstruct, fd);

            if isempty(paramobj.activematerial_type)            
                paramobj.activematerial_type = 'default';
            end

            am = 'ActiveMaterial';
            switch paramobj.activematerial_type
              case 'default'
                paramobj.(am) = ActiveMaterialInputParams(jsonstruct.(am));
              case 'composite'
                paramobj.(am) = CompositeActiveMaterialInputParams(jsonstruct.(am));
              otherwise
                error('activematerial_type not recognized');
            end
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
