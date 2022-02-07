classdef ElectroChemicalComponent < ElectronicComponent
    
    properties

        EffectiveDiffusionCoefficient % Effective diffusion coefficient
        
        
    end

    methods
        
        function model = ElectroChemicalComponent(paramobj)
        % Here, :code:`paramobj` is instance of :class:`ElectronicChemicalInputParams <Electrochemistry.ElectronicChemicalInputParams>`
        
            model = model@ElectronicComponent(paramobj);
            
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
