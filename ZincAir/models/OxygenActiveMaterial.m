classdef HydrogenActiveMaterial < SeaWaterActiveMaterial

    properties
        Asurf % 
    end
    
    methods
        
        function model = HydrogenActiveMaterial(inputparams)
            
            model = model@SeaWaterActiveMaterial(inputparams);
            fdnames = {'Asurf'};
            model = dispatchParams(model, inputparams, fdnames);
            
        end
        
        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@SeaWaterActiveMaterial(model);
            
            fn = @() HydrogenActiveMaterial.updateENernst;
            inputnames = {'T', 'cElectrolyte'};
            model = model.registerPropFunction({'ENernst', fn, inputnames});
            
            fn = @() HydrogenActiveMaterial.updateReactionRate;
            inputnames = {'eta', 'T'};
            model = model.registerPropFunction({'R', fn, inputnames});
            
        end
        
        function state = updateENernst(model, state)
            
            R = model.con.R;
            F = model.con.F;
            
            T = state.T;
            c = state.cElectrolyte;
            
            ENernst = 0 - R.*T./(2*F).*log(1000^2./c.^2);
            
            state.ENernst = ENernst;
        end
                
        function state = updateReactionRate(model, state)

            F     = model.con.F;
            Asurf = model.Asurf;
            
            eta = state.eta;
            T = state.T;
            
            R = Asurf./(2.*F).*SmoothButlerVolmerEquation(5e-4, 0.5, 2, eta, T, inf);
            
            state.R = R;
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
