classdef NormedCurrentCollector < CurrentCollector


    properties

        potentialDifferenceScaling
        
    end
    
    methods

        function model = NormedCurrentCollector(inputparams)

            model = model@CurrentCollector(inputparams);
            
            fdnames = {'potentialDifferenceScaling'};
            model = dispatchParams(model, inputparams, fdnames);
            
        end

        function model = registerVarAndPropfuncNames(model)

            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@CurrentCollector(model);

            varnames = {'phiRef', ...
                        'scaledDeltaPhi'};
            model = model.registerVarNames(varnames);


            fn = @NormedCurrentCollector.updateCurrent;
            inputnames = {'scaledDeltaPhi', 'conductivity'};
            model = model.registerPropFunction({'j', fn, inputnames});
            
            fn = @NormedCurrentCollector.updatePhi;
            inputnames = {'scaledDeltaPhi', 'phiRef'};
            model = model.registerPropFunction({'phi', fn, inputnames});

        end

        function state = updateCurrent(model, state)
        % Assemble electrical current which is stored in :code:`state.j`
            
            sigma          = state.conductivity;
            scaledDeltaPhi = state.scaledDeltaPhi;

            scaling = model.potentialDifferenceScaling;
            
            j = scaling*assembleHomogeneousFlux(model, scaledDeltaPhi, sigma);

            state.j = j;

        end

        function state = updatePhi(model, state)
        % Assemble electrical current which is stored in :code:`state.j`
            
            scaling = model.potentialDifferenceScaling;

            phiRef         = state.phiRef;
            scaledDeltaPhi = state.scaledDeltaPhi;

            state.phi = phiRef + scaling*scaledDeltaPhi;

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
