classdef FullSolidDiffusionSwellingModel < FullSolidDiffusionModel

    properties

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
            
            fn = @FullSolidDiffusionSwellingModel.updateMassSource;
            model = model.registerPropFunction({'massSource', fn, {'intercalationFlux'}});

            fn = @FullSolidDiffusionSwellingModel.updateRadius;
            model = model.registerPropFunction({'radius', fn, {'radiusElongation'}});

            fn = @FullSolidDiffusionSwellingModel.updateElongationRadius;
            model = model.registerPropFunction({'radiusElongation', fn, {'cAverage'}});
            
        end
        
        function state = updateRadius(model, state)
        % eq 4 in Modelling capacity fade in silicon-graphite composite electrodes for lithium-ion batteries
        % Shweta Dhillon, Guiomar Hernández, Nils P. Wagner, Ann Mari Svensson, Daniel Brandell ([ref 1])

            R_delith = model.particleRadius;
            cmax     = model.saturationConcentration;

            cAverage = state.cAverage;

            R_lith   = computeRadius(cAverage, cmax, R_delith);

            %vf = state.volumeFraction;
            %vf0 = 0.07;
            %R_lith = R_delith .* (vf./vf0).^(1/3);

            state.radius = R_lith;
            
        end

        function state = updateFlux(model, state)
            
            useDFunc = model.useDFunc;
            op = model.operators;
            r0 = model.particleRadius;
            
            c = state.c;
            D = state.D;
            radius = state.radius;

            beta = (radius/r0);

            if useDFunc
                state.flux = op.flux(((1./beta).^2).* D, c);
            else
                D = op.mapToParticle*D;
                beta = op.mapToParticle*beta;
                state.flux =  op.flux(((1./beta).^2).*D, c);
            end
            
            
        end

 
        function state = updateMassSource(model, state)
        % Modification of mass source according to eq 6 in  Analysis of Lithium Insertion/Deinsertion in a Silicon
        % Electrode Particle at Room Temperature Rajeswari Chandrasekaran, Alexandre Magasinski, Gleb Yushin, and
        % Thomas F. Fuller ([ref 3])

            op  = model.operators;
            rp0 = model.particleRadius;

            vf = state.volumeFraction;

            radius = state.radius;
            Rvol   = state.Rvol;

            % One can notice that Rvol.*((4*pi*rp0.^2 .* rp)./(3.*vf)) is strictly equal to the reaction 
            % rate for a non swelling particle (using a = 3*vf/rp), which is the rate i used in the 
            % equation 6. We then multiply it by the scalingCoeff defined in [ref3]

            massSource = - Rvol.*((4*pi* radius.^2 .* rp0)./(3*vf));
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
    
