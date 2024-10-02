classdef Binder < BaseModel

    properties

        %% Input Parameters

        %  Standard parameters
        
        electronicConductivity % the electronic conductivity of the material (symbol: sigma)
        density                % the mass density of the material (symbol: rho)
        massFraction           % the ratio of the mass of the material to the total mass of the phase or mixture (symbol: gamma)
        thermalConductivity    % Thermal conductivity of current collector
        specificHeatCapacity   % Heat capacity of current collector
        
    end

    methods
        
        function model = Binder(inputparams)
            
            model = model@BaseModel();
            
            fdnames = {'electronicConductivity', ... 
                       'density'               , ...                
                       'massFraction'          , ...
                       'thermalConductivity'   , ...
                       'specificHeatCapacity'};

            model = dispatchParams(model, inputparams, fdnames);
        end

        function jsonstruct = exportParams(model)

            jsonstruct = exportParams@BaseModel(model);

            fdnames = {'electronicConductivity', ... 
                       'density'               , ...                
                       'massFraction'          , ...           
                       'thermalConductivity'   , ...    
                       'specificHeatCapacity'};
            
            for ifd = 1 : numel(fdnames)
                fdname = fdnames{ifd};
                jsonstruct.(fdname) = model.(fdname);
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
