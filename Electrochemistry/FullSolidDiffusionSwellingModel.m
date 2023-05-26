classdef FullSolidDiffusionSwellingModel < FullSolidDiffusionModel

    properties

    end

    methods

        function model = FullSolidDiffusionSwellingModel(paramobj)
            model = model@FullSolidDiffusionModel(paramobj);
        end

        function model = registerVarAndPropfuncNames(model)

            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@FullSolidDiffusionModel(model);
                       
            fn = @FullSolidDiffusionSwellingModel.updateFlux;
            inputnames = {'c', 'D'};
            model = model.registerPropFunction({'flux', fn, inputnames});
            
            
            fn = @FullSolidDiffusionSwellingModel.updateMassSource;
            model = model.registerPropFunction({'massSource', fn, {'radius', 'volumeFraction', 'Rvol'}});
            
            fn = @FullSolidDiffusionSwellingModel.updateMassAccum;
            fn = {fn, @(propfunction) PropFunction.accumFuncCallSetupFn(propfunction)};
            model = model.registerPropFunction({'massAccum', fn, {'c'}});
            
            fn = @FullSolidDiffusionSwellingModel.assembleSolidDiffusionEquation;
            model = model.registerPropFunction({'solidDiffusionEq', fn, {'c', 'cSurface', 'massSource', 'D'}});

            
        end
        

        function state = updateMassSource(model, state)

        %% Modification of mass source
           
            op  = model.operators;
            rp0 = model.rp;
            amf = model.activeMaterialFraction;
            
            rp   = state.radius;
            vf   = state.volumeFraction;
            Rvol = state.Rvol;
            
            % Rvol and radius are discretized in  celltbl
            massSource = - Rvol.*((4*pi*rp.^3)./(3*amf.*vf)).*(rp0./rp).^3;
            massSource = op.mapFromBc*massSource;

            state.massSource = massSource;
            
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
    
