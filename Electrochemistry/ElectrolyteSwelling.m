classdef ElectrolyteSwelling < Electrolyte

    properties

        
    end

    methods

        function model = ElectrolyteSwelling(inputparams)
        % inputparams is instance of ElectrolyteInputParams or a derived class
            model = model@Electrolyte(inputparams);
        end

        %% Declaration of the Dynamical Variables and Function of the model (setup of varnameList and propertyFunctionList)
        function model = registerVarAndPropfuncNames(model)
            
            model = registerVarAndPropfuncNames@Electrolyte(model);

            varnames = {'volumeFraction'  ,...             % volume fraction of the electrolyte (= porosity of the SwellingMaterial)
                        'convFlux'};                       % Convective flux
            model = model.registerVarNames(varnames);

            
            fn = @ElectrolyteSwelling.updateDiffusionCoefficient;
            model = model.registerPropFunction({'D', fn, {'c','T','volumeFraction'}});

            fn = @ElectrolyteSwelling.updateConductivity;
            inputnames = {'c', 'volumeFraction'};
            model = model.registerPropFunction({'conductivity', fn, inputnames});

            fn = @ElectrolyteSwelling.updateMassFlux;
            model = model.registerPropFunction({'massFlux', fn, {'c', 'j', 'D', 'convFlux'}});
            
            fn = @ElectrolyteSwelling.updateConvFlux;
            model = model.registerPropFunction({'convFlux', fn, {}});

            fn = @ElectrolyteSwelling.updateVolumeFraction;
            model = model.registerPropFunction({'volumeFraction', fn, {}});
            
        end

        function state = updateConductivity(model, state)
            
            brcoef = model.bruggemanCoefficient;
            
            c  = state.c;
            T  = state.T;
            vf = state.volumeFraction;
            
            kappa = model.computeConductivity(c, T);
            state.conductivity = kappa.*vf.^brcoef;
            
        end

        %% Update the diffusion coefficient which vary at each step as the volumeFraction is no more constant
        function state = updateDiffusionCoefficient(model, state)

            brcoef = model.bruggemanCoefficient; 
            computeD = model.computeDiffusionCoefficient;

            c  = state.c;
            T  = state.T;
            vf = state.volumeFraction;
            
            D = computeD(c, T);
            
            % set effective coefficient
            state.D = D .* vf .^ brcoef;

        end

        %% Update the mass flux which is the sum of the three fluxes defines above.
        function state = updateMassFlux(model, state)

            state = updateMassFlux@Electrolyte(model, state);
            
            convFlux = state.convFlux;

            state.massFlux = state.massFlux + convFlux;

        end

        %%By default definitions. The real values are set in the BatterySwelling class (respectively in updateElectrolyteVolumeFraction and
        %% updateConvFlux functions)

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
