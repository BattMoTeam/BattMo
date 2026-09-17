classdef ProtonicMembraneStokesGasSupply < ProtonicMembraneGasSupply
    
    properties


    end
    
    methods
        
        function model = ProtonicMembraneStokesGasSupply(inputparams)
            
            model = model@ProtonicMembraneGasSupply(inputparams);

        end


        function model = registerVarAndPropfuncNames(model)
            
            model = registerVarAndPropfuncNames@ProtonicMembraneGasSupply(model);
            
            nGas = model.nGas;

            varnames = {};
            % Velocities / m/s (the 2D spatial indexing follows from the stokes solver setup)
            varnames = {'veloticies'};
            % Stokes equation
            varnames = {'stokesEquation'};
            % Total fluxes at the grid faces
            varnames = {'totalFlux'};
            
            model = model.registerVarNames(varnames);

            fn = @ProtonicMembraneStokesGasSupply.setupStokesEquation;
            inputvarnames = {'velocities', 'pressure'};
            model = model.registerPropFunction({'stokesEquation', fn, inputvarnames});
            
            fn = @ProtonicMembraneStokesGasSupply.updateTotalFlux;
            inputvarnames = {'velocities', 'pressure'};
            model = model.registerPropFunction({'totalFlux', fn, inputvarnames});
            
            for igas = 1 : nGas

                fn = @ProtonicMembraneStokesGasSupply.setupMassFluxes;
                inputvarnames = {VarName({}, 'massfractions', nGas, igas), ...
                                 VarName({}, 'densities', nGas, igas), ...
                                 'density', 'totalFlux'};
                outputvarname = VarName({}, 'massFluxes', nGas, igas);
                model = model.registerPropFunction({outputvarname, fn, inputvarnames});
                
            end
            
        end


        
    end
    
end

%{
Copyright 2021-2026 SINTEF Industry, Sustainable Energy Technology
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
