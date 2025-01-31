classdef ZincAirActiveMaterial < BaseModel

    properties
        
        con = PhysicalConstants();
        chargeCarrierName
        stochElectron
        etaMax
    
    end
    
    methods
        
        function model = ZincAirActiveMaterial(inputparams)

            model = model@BaseModel();
            fdnames = {'chargeCarrierName', ...
                       'stochElectron', ...
                       'etaMax'};
            model = dispatchParams(model, inputparams, fdnames);
            
        end

        
        function model = registerVarAndPropfuncNames(model)
            
            varnames = {};
            % Temperature
            varnames{end + 1} = 'T'; 
            % Reaction rate (in mol/s/m^3)
            varnames{end + 1} = 'R'; 
            % electric potential from electrode 
            varnames{end + 1} = 'phiElectrode'; 
            % electric potential from electrolyte 
            varnames{end + 1} = 'phiElectrolyte'; 
            % concentration from electrolyte
            varnames{end + 1} = 'cElectrolyte'; 
            % value for eta (overpotential)
            varnames{end + 1} = 'eta'; 
            % value of ENernst coefficient
            varnames{end + 1} = 'ENernst'; 

            model = model.registerVarNames(varnames);

            fn = @() ZincAirActiveMaterial.updateEta;
            inputnames = {'phiElectrolyte', 'phiElectrode', 'ENernst'};
            model = model.registerPropFunction({'eta', fn, inputnames});
            
        end
        
        
        function state = updateEta(model, state)          

            phiElectrolyte = state.phiElectrolyte;
            phiElectrode   = state.phiElectrode;
            ENernst        = state.ENernst;

            eta = phiElectrode - phiElectrolyte - ENernst;
            
            state.eta = eta;
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
