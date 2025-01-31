classdef OxygenElectrode < ZincElectrode

    
    methods
        
        function model = registerVarAndPropfuncNames(model)
            
            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)
            
            model = registerVarAndPropfuncNames@ZincElectrode(model);
            
            model = model.registerVarName('E');
            
            fn = @() ZincElectrode.updatePotential;
            inputnames = {'E'};
            model = model.registerPropFunction({'phi', fn, inputnames});
            
       end

       
        function state = updatePotential(model, state)
        % this is a special setup : We impose given potential and assume infinite conductivity in hydrogen electrode
        % (likely to change later)
            nc = model.G.getNumberOfCells();
            
            E = state.E;
            
            state.phi = E.*ones(nc, 1);
        end
        
        function state = updateVolumeFraction(model, state)
            nc = model.G.getNumberOfCells();
            
            state.volumeFraction = model.volumeFraction*ones(nc, 1);
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
