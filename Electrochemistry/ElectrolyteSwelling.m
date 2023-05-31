classdef ElectrolyteSwelling < Electrolyte

    properties

    end

    methods

        function model = ElectrolyteSwelling(paramobj)
        % paramobj is instance of ElectrolyteInputParams or a derived class
            model = model@Electrolyte(paramobj);
        end

%% Declaration of the Dynamical Variables and Function of the model (setup of varnameList and propertyFunctionList)
        function model = registerVarAndPropfuncNames(model)
       
            model = registerVarAndPropfuncNames@Electrolyte(model);


            varnames = {'volumeFraction'  ,...             % volume fraction of the electrolyte (= porosity of the SwellingMaterial)
                        'convFlux'};                       % Convective flux
            model = model.registerVarNames(varnames);


            fn = @ElectrolyteSwelling.updateAccumTerm;
            model = model.registerPropFunction({'massAccum', fn, {'c','volumeFraction'}});

            fn = @ElectrolyteSwelling.updateDiffusionCoefficient;
            model = model.registerPropFunction({'D', fn, {'c','T','vf'}});

            fn = @ElectrolyteSwelling.updateCurrent;
            model = model.registerPropFunction({'j', fn, {'dmudcs','phi','T','c','conductivity','vf'}});

            fn = @ElectrolyteSwelling.updateMassFlux;
            model = model.registerPropFunction({'massFlux', fn, {'c','j','D'}});
            
            fn = @ElectrolyteSwelling.updateConvFlux;
            model = model.registerPropFunction({'convFlux', fn, {}});

            fn = @ElectrolyteSwelling.updateVolumeFraction;
            model = model.registerPropFunction({'volumeFraction', fn, {}});
        end

    %% Definition of the accumulation term (dc/dt)
        function state = updateAccumTerm(model, state, state0, dt)

            c = state.c;
            vf = state.volumeFraction;

            c0 = state0.c;

            cdotcc  = (c - c0)/dt;
            effectiveVolumes = vf.*model.G.cells.volumes;

            state.massAccum  = effectiveVolumes.*cdotcc;
            
        end

    %% Update the diffusion coefficient which vary at each step as the volumeFraction is no more constant
        function state = updateDiffusionCoefficient(model, state)

            brcoef = model.BruggemanCoefficient; 
            computeD = model.computeDiffusionCoefficientFunc;

            c = state.c;
            T = state.T;
            vf = state.volumeFraction;
            
            D = computeD(c, T);
            
            % set effective coefficient
            state.D = D .* vf .^ brcoef;

        end

    %% Update the current in the electrolyte (due to the movement of the ions)
        function state  = updateCurrent(model, state)

            ncomp  = model.ncomp;
            sp     = model.sp;
            R      = model.constants.R;
            F      = model.constants.F;
            brcoef = model.BruggemanCoefficient;
            
            dmudcs       = state.dmudcs;
            phi          = state.phi;
            T            = state.T;
            c            = state.c;
            conductivity = state.conductivity;
            vf = state.volumeFraction;

            % Compute effective ionic conductivity in porous media
            conductivityeff = conductivity.*vf.^brcoef;

            state.conductivityeff = conductivityeff;
            j = assembleFlux(model, phi, conductivityeff);

            sum_dmudc = dmudcs{1} + dmudcs{2};
            coef = (1/F)*(1 - sp.t(1))*conductivityeff.*sum_dmudc;
            jchem = assembleFlux(model, c, coef);

            j = j - jchem;

            state.j = j;

        end

    %% Update the mass flux which is the sum of the three fluxes defines above.
        function state = updateMassFlux(model, state)


            sp = model.sp;
            
            % We assume that the current and the diffusion coefficient D has been updated when this function is called
            c = state.c;
            j = state.j;
            D = state.D;

            %% 1. Flux from diffusion
            diffFlux = assembleFlux(model, c, D);
            state.diffFlux = diffFlux;

            %% 2. Flux from electrical forces
            F = model.constants.F;
            fluxE = sp.t ./ (sp.z .* F) .* j;

            %% 3. Flux from the convective term to be taken into account for swelling materials(cf eq 6, paper Analysis of the Lithium-ion Insertion
            %Silicon composite electrode/separator/lithium foil cell of Rajeswari Chandrasekaran and Thomas F. Fuller
            convFlux = state.convFlux;

            %% 3. Sum the two flux contributions
            flux = diffFlux + fluxE + convFlux;

            state.massFlux = flux;

        end

        %%By default definitions. The real values are set in the BatterySwelling class (respectively in updateElectrolyteVolumeFraction and
        %%updateConvFlux functions)

        function state = updateVolumeFraction(model, state)
            state.volumeFraction = model.volumeFraction;
        end

        function state = updateConvFlux(model, state)
            state.convFlux = 0;
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
