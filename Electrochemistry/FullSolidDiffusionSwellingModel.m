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
 % Modification of mass source according to eq 6 in  Analysis of Lithium Insertion/Deinsertion in a Silicon
 % Electrode Particle at Room Temperature Rajeswari Chandrasekaran, Alexandre Magasinski, Gleb Yushin, and
 % Thomas F. Fuller ([ref 3])

            op   = model.operators;
            rp0  = model.rp;
            amf = model.activeMaterialFraction;

            radius   = state.radius;
            vf   = state.volumeFraction;
            Rvol = state.Rvol;
            
            % One can notice that Rvol.*((4*pi*rp0.^2 .* rp)./(3.*vf)) is strictly equal to the reaction 
            % rate for a non swelling particle (using a = 3*vf/rp), which is the rate i used in the 
            % equation 6. We then multiply it by the scalingCoeff defined in [ref3]
            %scalingCoeff = (rp0./radius).^3;
            scalingCoeff = (rp0./radius).^3;

            massSource = - Rvol.*((4*pi* rp0.^2 .* radius)./(3*vf.*amf)).* scalingCoeff;
            massSource = op.mapFromBc*massSource;

            state.massSource = massSource;
        end


   %% The other variables are directly updated in methods defined in the class SwellingMaterial because they need parameters unavailable in this class
   
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
    
