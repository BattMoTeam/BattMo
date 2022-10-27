classdef CurrentCollector < ElectronicComponent

    properties
        
        externalCouplingTerm

        thermalConductivity
        heatCapacity
        density
        
    end
    
    methods
        
        function model = CurrentCollector(paramobj)
            
            model = model@ElectronicComponent(paramobj);

            fdnames = {'externalCouplingTerm', ...
                       'thermalConductivity' , ...
                       'heatCapacity'        , ...
                       'density'};
            
            model = dispatchParams(model, paramobj, fdnames);
            
            % The parameter EffectiveElectricalConductivity in CurrentCollectorInputParams is given as scalar
            model.EffectiveElectricalConductivity = model.EffectiveElectricalConductivity*ones(model.G.cells.num, 1);
            
        end
        
        function model = registerVarAndPropfuncNames(model)

            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)
            
            model = registerVarAndPropfuncNames@ElectronicComponent(model);
            
            varnames = {'jCoupling', ...
                        'jExternal', ...
                        'jFaceCoupling', ...
                        'jFaceExternal'};
            model = model.registerVarNames(varnames);
            
            fn = @CurrentCollector.updatejBcSource;
            model = model.registerPropFunction({'jBcSource', fn, {'jCoupling', 'jExternal'}});

            if model.use_thermal
                fn = @CurrentCollector.updatejFaceBc;
                model = model.registerPropFunction({'jFaceBc', fn, {'jFaceCoupling', 'jFaceExternal'}});
            end
            
        end
        
        function state = updatejBcSource(model, state)
            
            state.jBcSource = state.jCoupling + state.jExternal;
            
        end
        
        function state = updatejFaceBc(model, state)
            
            state.jFaceBc = state.jFaceCoupling + state.jFaceExternal;
            
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
