classdef SimplifiedSolidDiffusionModel < SolidDiffusionModel
%
% See simplified model presented in *Comparison of approximate solution methods for the solid phase diffusion equation in a porous electrode model* by Zhang, Qi and White, Ralph E (
% Journal of power sources, 2007) :cite:p:`zhang2007comparison`
%
    methods

        function model = SimplifiedSolidDiffusionModel(inputparams)

            model = model@SolidDiffusionModel(inputparams);

        end

        function model = registerVarAndPropfuncNames(model)

            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@SolidDiffusionModel(model);

            varnames = {};
            % concentration at surface
            varnames{end + 1} = 'cSurface';
            % concentration (over the whole particle)
            varnames{end + 1} = 'cAverage';
            % Solid diffusion equation
            varnames{end + 1} = 'solidDiffusionEq';

            model = model.registerVarNames(varnames);

            fn = @SimplifiedSolidDiffusionModel.updateDiffusionCoefficient;
            inputnames = {'T'};
            model = model.registerPropFunction({'D', fn, inputnames});

            fn = @SimplifiedSolidDiffusionModel.assembleSolidDiffusionEquation;
            inputnames = {'cSurface', 'cAverage', 'Rvol', 'D'};
            model = model.registerPropFunction({'solidDiffusionEq', fn, inputnames});

            fn = @SimplifiedSolidDiffusionModel.assembleAccumTerm;
            fn = {fn, @(propfunction) PropFunction.accumFuncCallSetupFn(propfunction)};
            inputnames = {'cAverage'};
            model = model.registerPropFunction({'massAccum', fn, inputnames});

            fn = @SimplifiedSolidDiffusionModel.updateMassConservation;
            inputnames = {'massAccum', 'massSource'};
            model = model.registerPropFunction({'massCons', fn, inputnames});

            fn = @SimplifiedSolidDiffusionModel.updateMassSource;
            inputnames = {'Rvol'};
            model = model.registerPropFunction({'massSource', fn, inputnames});

        end


        function state = updateMassConservation(model, state)
        % Used when diffusionModelType == 'simple' or no particle diffusion

            source = state.massSource;
            accum = state.massAccum;

            state.massCons = accum - source;

        end


        function state = assembleAccumTerm(model, state, state0, dt)

            vf = model.volumeFraction;

            c  = state.cAverage;
            c0 = state0.cAverage;

            state.massAccum = vf.*(c - c0)/dt;

        end

        function state = updateMassSource(model, state)
        % used when diffusionModelType == simple

            Rvol = state.Rvol;

            state.massSource = - Rvol;

        end


        function state = assembleSolidDiffusionEquation(model, state)
        % We update the surface concentration of the charge carrier in the active material.
        % The surface concentration value is computed following polynomial method, as described in ref1 (see below)

            rp  = model.particleRadius;
            vsa = model.volumetricSurfaceArea;
            
            csurf = state.cSurface;
            caver = state.cAverage;
            D     = state.D;
            Rvol  = state.Rvol;

            state.solidDiffusionEq = csurf - caver + (rp.*Rvol)./(5*vsa*D);

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
