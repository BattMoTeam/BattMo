classdef Separator < BaseModel
    
    properties
        
        % Physicochemical properties
        porosity       % Porosity,         [-]
        rp             % Pore radius,      [m]
        Gurley         % Gurley number,    [s]
        volumeFraction % Volume fraction,  [-]
        
        thermalConductivity % intrinsic thermal conductivity value
        heatCapacity % intrinsic heat capacity value

        EffectiveThermalConductivity
        EffectiveHeatCapacity
        density % [kg m^-3]        
    end
    
    methods
        function model = Separator(paramobj)
            
            model = model@BaseModel();

            % OBS : All the models should have same backend (this is not assigned automaticallly for the moment)
            % in the case of the separator, probably this does not matter as no computation is actually done on this grid
            model.AutoDiffBackend = SparseAutoDiffBackend('useBlocks', false);
            
            fdnames = {'G'                  , ...
                       'porosity'           , ...
                       'rp'                 , ...
                       'Gurley'             , ...
                       'thermalConductivity', ...
                       'heatCapacity'       , ...
                       'density'};
            model = dispatchParams(model, paramobj, fdnames);
            
            model.volumeFraction = 1 - model.porosity;
            model.EffectiveThermalConductivity = model.thermalConductivity.*(model.volumeFraction).^1.5;
            model.EffectiveHeatCapacity = model.heatCapacity.*model.volumeFraction;


            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)
            
            % no dynamical variable in the separator
            
        end
    end
    
end




%{
Copyright 2009-2021 SINTEF Industry, Sustainable Energy Technology
and SINTEF Digital, Mathematics & Cybernetics.

This file is part of The Battery Modeling Toolbox BatMo

BatMo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

BatMo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with BatMo.  If not, see <http://www.gnu.org/licenses/>.
%}
