classdef ZincActiveMaterial < ZincAirActiveMaterial

    
    properties
        Asurf
    end
    
    methods

        function model = ZincActiveMaterial(inputparams)
            
            model = model@ZincAirActiveMaterial(inputparams);
            fdnames = {'Asurf', ...
                       'etaMax'};
            model = dispatchParams(model, inputparams, fdnames);
            
        end
        
        function model = registerVarAndPropfuncNames(model)
        
            model = registerVarAndPropfuncNames@ZincAirActiveMaterial(model);
        
            varnames = {};
            
            varnames{end + 1} = 'specificSurfaceArea';
            
            model = model.registerVarNames(varnames);
            
            fn = @ZincActiveMaterial.updateENernst;
            inputnames = {'T', 'cElectrolyte'};
            model = model.registerPropFunction({'ENernst', fn, inputnames});
            
            fn = @ZincActiveMaterial.updateReactionRate;
            inputnames = {'eta', 'T', 'specificSurfaceArea'};
            model = model.registerPropFunction({'R', fn, inputnames}); 
            
        end
        
       
        function state = updateENernst(model, state)
            
            R = model.con.R;
            F = model.con.F;
            
            T = state.T;
            c = state.cElectrolyte;
            
            ENernst = -2.7 - R.*T./(2.*F).*log(1000./c);
            
            state.ENernst = ENernst;

        end
        
        function state = updateReactionRate(model, state)
            
            R = model.con.R;
            F = model.con.F;
            etaMax = model.etaMax;
            
            eta   = state.eta;
            T     = state.T;
            Asurf = state.specificSurfaceArea;
            
            R = Asurf./(2*F).*SmoothButlerVolmerEquation(1e-6, 0.5, 2, eta, T, etaMax); % mol m^-3 s^-1
            
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
