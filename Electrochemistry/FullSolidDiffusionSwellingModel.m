classdef FullSolidDiffusionSwellingModel < FullSolidDiffusionModel

    properties

        molarVolumeSi = 1.2e-05;
        molarVolumeLi = 8.8e-06;
        
    end

    methods

        function model = FullSolidDiffusionSwellingModel(inputparams)

            model = model@FullSolidDiffusionModel(inputparams);
            
        end

        function model = registerVarAndPropfuncNames(model)

            sd = 'SolidDiffusion';

            model = registerVarAndPropfuncNames@FullSolidDiffusionModel(model);

            varnames = {'radiusElongation', ...
                        'radius'          , ...
                        'intercalationFlux'};

            model = model.registerVarNames(varnames);

            model = model.unsetAsExtraVarName('cAverage');

            model = model.setAsExtraVarName('radius'); % the radius is not needed directly. We use the radius elongation value instead.
            
            fn = @FullSolidDiffusionSwellingModel.updateMassSource;
            model = model.registerPropFunction({'massSource', fn, {'intercalationFlux', 'radiusElongation'}});

            fn = @FullSolidDiffusionSwellingModel.updateRadius;
            model = model.registerPropFunction({'radius', fn, {'radiusElongation'}});

            fn = @FullSolidDiffusionSwellingModel.updateRadiusElongation;
            model = model.registerPropFunction({'radiusElongation', fn, {'cAverage'}});

            fn = @FullSolidDiffusionSwellingModel.updateMassAccum;
            fn = {fn, @(propfunction) PropFunction.accumFuncCallSetupFn(propfunction)};            
            model = model.registerPropFunction({'massAccum', fn, {'c', 'radiusElongation'}});
            
        end

        function state = updateMassAccum(model, state, state0, dt)

            op = model.operators;

            c      = state.c;
            delta  = op.mapToParticle*state.radiusElongation;
            c0     = state0.c;
            delta0 = op.mapToParticle*state0.radiusElongation;

            state.massAccum = 1./(delta.*dt).*op.vols.*(delta.^3.*c - delta0.^3.*c0);
            
        end
        
        function state = updateRadius(model, state)

            state.radius = model.particleRadius*state.radiusElongation;
            
        end
        
        function state = updateRadiusElongation(model, state)
        % eq 4 in Modelling capacity fade in silicon-graphite composite electrodes for lithium-ion batteries
        % Shweta Dhillon, Guiomar Hernández, Nils P. Wagner, Ann Mari Svensson, Daniel Brandell ([ref 1])

            cmax = model.saturationConcentration;

            cAverage = state.cAverage;

            Q = (3.75.*model.molarVolumeLi)./(model.molarVolumeSi);

            state.radiusElongation = (1 + Q.*cAverage./cmax).^(1/3);
            
        end

        function state = updateMassSource(model, state)
        % Modification of mass source according to eq 6 in  Analysis of Lithium Insertion/Deinsertion in a Silicon
        % Electrode Particle at Room Temperature Rajeswari Chandrasekaran, Alexandre Magasinski, Gleb Yushin, and
        % Thomas F. Fuller ([ref 3])

            op  = model.operators;
            
            rp0  = model.particleRadius;
            vf0  = model.volumeFraction;
            vsa0 = model.volumetricSurfaceArea;
            
            delta = state.radiusElongation;
            R     = state.intercalationFlux;

            massSource = - (R./(delta.^2)).*((4/3)*pi*rp0^3.*vsa0/vf0);
            massSource = op.mapFromBc*massSource;

            state.massSource = massSource;
            
        end

        function cleanState = addStaticVariables(model, cleanState, state)

            cleanState.radiusElongation = state.radiusElongation;
            
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
    
