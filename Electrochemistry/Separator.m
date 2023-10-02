classdef Separator < BaseModel
    
    properties
        
        porosity             % the ratio of the volume free space to the total volume (symbol: varepsilon)
        density              % the mass density of the material (symbol: rho)
        bruggemanCoefficient % coefficient to determine effective transport parameters in porous media (symbol: beta)
        
        thermalConductivity  % Intrinsic Thermal conductivity of the electrolyte
        specificHeatCapacity % Specific Heat capacity of the electrolyte
        EffectiveThermalConductivity
        EffectiveVolumetricHeatCapacity

        use_thermal

    end
    
    methods

        function model = Separator(paramobj)
            
            model = model@BaseModel();

            % OBS : All the models should have same backend (this is not assigned automaticallly for the moment)
            % in the case of the separator, probably this does not matter as no computation is actually done on this grid
            model.AutoDiffBackend = SparseAutoDiffBackend('useBlocks', false);
            
            fdnames = {'G'                   , ...
                       'porosity'            , ...
                       'density'             , ...
                       'bruggemanCoefficient', ...
                       'thermalConductivity' , ...
                       'specificHeatCapacity', ...
                       'use_thermal'};
            model = dispatchParams(model, paramobj, fdnames);
            model.porosity = model.porosity*ones(model.G.cells.num,1);
            model = model.setupDependentProperties();

        end
        
        function model = setupDependentProperties(model)
            model.volumeFraction = 1 - model.porosity;

            if model.use_thermal
                model.EffectiveThermalConductivity = model.thermalConductivity.*(model.volumeFraction).^model.BruggemanCoefficient;
                model.EffectiveVolumetricHeatCapacity = model.specificHeatCapacity.*model.volumeFraction*model.density;
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
