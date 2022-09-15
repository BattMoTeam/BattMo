classdef Separator < BaseModel
    
    properties
        
        porosity            % Porosity [-]
        thermalConductivity % intrinsic thermal conductivity value
        heatCapacity        % intrinsic heat capacity value
        density             % Density [kg m^-3]        
        
        volumeFraction      % Volume fraction [-]
        EffectiveThermalConductivity
        EffectiveHeatCapacity

        use_thermal

        BruggemanCoefficient
    end
    
    methods

        function model = Separator(paramobj)
            
            model = model@BaseModel();

            % OBS : All the models should have same backend (this is not assigned automaticallly for the moment)
            % in the case of the separator, probably this does not matter as no computation is actually done on this grid
            model.AutoDiffBackend = SparseAutoDiffBackend('useBlocks', false);
            
            fdnames = {'G'                  , ...
                       'porosity'           , ...
                       'thermalConductivity', ...
                       'heatCapacity'       , ...
                       'density'            , ...
                       'use_thermal'        , ...
                       'BruggemanCoefficient'};
            model = dispatchParams(model, paramobj, fdnames);
            model.porosity = model.porosity*ones(model.G.cells.num,1);
            model = model.setupDependentProperties();

        end
        
        function model = setupDependentProperties(model)
            model.volumeFraction = 1 - model.porosity;

            if model.use_thermal
                model.EffectiveThermalConductivity = model.thermalConductivity.*(model.volumeFraction).^model.BruggemanCoefficient;
                model.EffectiveHeatCapacity = model.heatCapacity.*model.volumeFraction;
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
