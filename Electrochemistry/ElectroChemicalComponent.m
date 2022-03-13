classdef ElectroChemicalComponent < ElectronicComponent
    
    properties

        EffectiveDiffusionCoefficient % Effective diffusion coefficient
        
    end

    methods
        
        function model = ElectroChemicalComponent(paramobj)
        % Here, :code:`paramobj` is instance of :class:`ElectronicChemicalInputParams <Electrochemistry.ElectronicChemicalInputParams>`
        
            model = model@ElectronicComponent(paramobj);

        end

        
        function model = registerVarAndPropfuncNames(model)
        
            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)
            
            model = registerVarAndPropfuncNames@ElectronicComponent(model);

            varnames = {'c'         , ...
                        'massSource', ...
                        'massFlux'  , ...
                        'massAccum' , ...
                        'massCons'};
            model = model.registerVarNames(varnames);
            
            fn = @ElectroChemicalComponent.updateMassFlux;
            model = model.registerPropFunction({'massFlux', fn, {'c'}});
            
            fn = @ElectroChemicalComponent.updateMassConservation;
            model = model.registerPropFunction({'massCons', fn, {'massFlux', 'massSource', 'massAccum'}});
            
            
        end
        
        function state = updateMassFlux(model, state)
        % Assemble diffusion flux which is stored in :code:`state.Flux`

            D = model.EffectiveDiffusionCoefficient;
            
            c = state.c;

            massflux = assembleFlux(model, c, D);
            
            state.massFlux = massflux;
            
        end
        
        function state = updateMassConservation(model, state)
        % Assemble residual of the mass conservation equation which is stored in :code:`state.massCons`
            
            flux   = state.massFlux;
            source = state.massSource;
            accum  = state.massAccum;
            bcsource = 0;
            
            masscons = assembleConservationEquation(model, flux, bcsource, source, accum);
            
            state.massCons = masscons;
            
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
