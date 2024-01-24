classdef SeaWaterElectrolyte < SeaWaterElectrolyteNoPrecipitation

    %% Model with one precipitate

    properties

        dischargeProductMolarVolume
        superOversaturationRatio
        solidPrecipitatePorosity
        characteristicPoreDiameter
        nucleationMaximum
        nucleationRate
        nucleationActivation

    end

    methods

        function model = SeaWaterElectrolyte(inputparams)
            model = model@SeaWaterElectrolyteNoPrecipitation(inputparams);

            fdnames = {'dischargeProductMolarVolume', ...
                       'superOversaturationRatio'   , ...
                       'solidPrecipitatePorosity'   , ...
                       'nucleationMaximum'          , ...
                       'nucleationRate'             , ...
                       'nucleationActivation'       , ...
                       'characteristicPoreDiameter'};
            model = dispatchParams(model, inputparams, fdnames);

            model = model.setupReactionRates(inputparams);

            model = model.setupMainIonIndex();
        end

        function model = setupMainIonIndex(model)
            error('Virtual function. It depends on the chosen electrolyte');
        end


        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@SeaWaterElectrolyteNoPrecipitation(model);

            varnames = {};

            % Solid volume fraction (composed of precipitate)
            varnames{end + 1} = 'solidVolumeFraction';
            % Solid volume fraction (composed of precipitate)
            varnames{end + 1} = 'availableVolumeFraction';
            % Reaction rate for the precipitation equation
            varnames{end + 1} = 'Rprecipitation';
            % Nucleation variable (when equal to 1 the nucle is formed)
            varnames{end + 1} = 'nucleation';
            % Saturation constant for 'Mg+2'
            varnames{end + 1} = 'cSat';
            % residual of the evolution equation for the 'nucleation' variable
            varnames{end + 1} = 'nucleationEquation';
            % accumulation term for the discharge
            varnames{end + 1} = 'dischargeAccum';
            % source term for the discharge equation
            varnames{end + 1} = 'dischargeSource';
            % residual for the mass conservation of discharge.
            varnames{end + 1} = 'dischargeMassCons';
            % driving force in computation of precipitation rate
            varnames{end + 1} = 'k';
            % Indicator use in computation of time stepping
            varnames{end + 1} = 'indicator';
            %
            varnames{end + 1} = 'volumeFractionEquation';

            model = model.registerVarNames(varnames);

            spdict    = model.spdict;
            nsp       = model.nsp;
            nqp       = model.nqp;
            indmainsp = model.indmainsp;

            %% update Atomic Mass Conservation equation residual
            fn = @() SeaWaterElectrolyte.updateAtomicMassCons;
            for inqp = 1 : nqp
                inputnames = {VarName({}, 'qpepscs', nqp, inqp), VarName({}, 'cs', nsp), 'volumeFraction'};
                model = model.registerPropFunction({VarName({}, 'atomicMassCons', nqp, inqp), fn, inputnames});
            end

            %% update discharge source value
            fn = @() SeaWaterElectrolyte.updateDischargeSource;
            inputnames = {'Rprecipitation'};
            model = model.registerPropFunction({{'dischargeSource'}, fn, inputnames});

            %% update discharge mass conservation equation residual
            fn = @() SeaWaterElectrolyte.updateDischargeMassCons;
            inputnames = {'dischargeAccum', ...
                          'dischargeSource'};
            model = model.registerPropFunction({{'dischargeMassCons'}, fn, inputnames});

            %% update solid volume fraction variables
            fn = @() SeaWaterElectrolyte.updateSolidVolumeFraction;
            indsolid = model.indsolidsp(1); %
            inputnames = {VarName({}, 'cs', nsp, indsolid)};
            model = model.registerPropFunction({{'solidVolumeFraction'}, fn, inputnames});

            %% update nucleation equation residual
            % note : this function uses state from previous time step and time step value in its input
            indi = model.mainIonIndex;
            fn = @() SeaWaterElectrolyte.assembleNucleationEquation;
            fn = {fn, @(propfunc) PropFunction.accumFuncCallSetupFn(propfunc)};
            inputnames = {'cSat', ...
                          VarName({}, 'cs', nsp, indi), ...
                          'availableVolumeFraction'
                          };
            model = model.registerPropFunction({{'nucleationEquation'}, fn, inputnames});


            fn = @() SeaWaterElectrolyte.updateAvailableVolumeFraction;
            inputnames = {'volumeFraction', ...
                          'solidVolumeFraction'
                          };
            model = model.registerPropFunction({'availableVolumeFraction', fn, inputnames});

            %% update precipitation reaction rate
            fn = @() SeaWaterElectrolyte.updateRprecipitation;
            inputnames = {'solidVolumeFraction', ...
                          'volumeFraction', ...
                          'nucleation', ...
                          'k'};
            model = model.registerPropFunction({{'Rprecipitation'}, fn, inputnames});

            fn = @() SeaWaterElectrolyte.updateDischargeAccumTerm;
            fn = {fn, @(propfunc) PropFunction.accumFuncCallSetupFn(propfunc)};
            indsolid = model.indsolidsp(1);
            nsp = model.nsp;
            inputnames = {VarName({}, 'cs', nsp, indsolid)};
            model = model.registerPropFunction({'dischargeAccum', fn, inputnames});

            %% following commented functions are not used in assembly (not part of the residual assembly equations)
            %% update an indicator that is (or can be) used in computation of time-step size
            fn = @() SeaWaterElectrolyte.computeTimeStepIndicator;
            inputnames = {'volumeFraction'};
            model = model.registerPropFunction({{'indicator'}, fn, inputnames});

            model = model.setAsExtraVarName('indicator');

        end

        function state = updateDischargeAccumTerm(model, state, state0, dt)

            indsolid = model.indsolidsp(1);

            vols = model.G.getVolumes();
            state.dischargeAccum = 1/dt*vols.*(state.cs{indsolid} - state0.cs{indsolid});

        end

        function state = updateRprecipitation(model, state)

            spp      = model.solidPrecipitatePorosity;
            porediam = model.characteristicPoreDiameter;

            svf = state.solidVolumeFraction;
            vf  = state.volumeFraction;
            nuc = state.nucleation;
            k   = state.k;
            cs  = state.cs;

            freeVolumeFraction = vf - spp*svf;
            fvfSmoother = exp(-1e-1./freeVolumeFraction);
            fvfSmoother(value(freeVolumeFraction) < 0) = 0;

            precipitationSurfaceArea = (6./porediam).*nuc.*(freeVolumeFraction.*fvfSmoother);

            R = precipitationSurfaceArea.*k;

            indsolid = model.indsolidsp;
            csolid = cs{indsolid};

            ind = value(R) < 0;
            if any(ind)
                R(ind) = R(ind).*csolid(ind);
            end

            state.precipitationSurfaceArea = precipitationSurfaceArea;
            state.Rprecipitation = R;
            % debugging : we turn off precipitation
            % state.Rprecipitation = 0;

        end

        function state = assembleNucleationEquation(model, state, state0, dt)

            indi     = model.mainIonIndex;
            op       = model.operators;
            vols     = model.G.getVolumes();
            osr      = model.superOversaturationRatio;
            nucMax   = model.nucleationMaximum;
            nucRate  = model.nucleationRate;
            nucAct   = model.nucleationActivation;

            cs         = state.cs;
            cSat       = state.cSat;
            nucleation = state.nucleation;
            avf        = state.availableVolumeFraction;

            c = cs{indi};

            % nucleation regime
            ind1 = (value(c) >= osr*value(cSat));
            % inactive regime
            ind2 = (value(c) < osr*value(cSat)) & (value(c) >= value(cSat));
            % denucleation regime
            ind3 = (value(c) < value(cSat));

            nucMax = avf.*nucMax;

            rate = 0*nucleation; % initialization;

            if any(ind1)
                rate(ind1) = nucRate./nucMax(ind1).*exp(-nucAct./(log(c(ind1)./(osr*cSat(ind1))))).*(nucMax(ind1) - nucleation(ind1));
            end

            if any(ind3)
                rate(ind3) = - nucRate./nucMax(ind3).*exp(-nucAct./(log(cSat(ind3)./c(ind3)))).*nucleation(ind3);
            end

            eq = 1/dt*vols.*(nucleation - state0.nucleation)  - vols.*rate;

            state.nucleationEquation = eq;

        end

        function state = updateAvailableVolumeFraction(model, state)

            state.availableVolumeFraction = state.volumeFraction + state.solidVolumeFraction;

        end


        function state = updateSolidVolumeFraction(model, state)

            molarVolume = model.dischargeProductMolarVolume;
            indsolid = model.indsolidsp(1); %

            state.solidVolumeFraction = state.cs{indsolid}*molarVolume;

        end


        function state = updateAtomicMassCons(model, state)
        % We add solid part in the quasi particles

            state = updateAtomicMassCons@SeaWaterElectrolyteNoPrecipitation(model, state);

            nqp         = model.nqp;
            C           = model.qpCompositionMatrix;
            indsolidsp  = model.indsolidsp;

            atomicMassCons = state.atomicMassCons;
            cs = state.cs;

            for indqp = 1 : nqp
                for isolid = 1 : model.nsolid
                    indsp = indsolidsp(isolid);
                    if C(indqp, indsp) ~= 0
                        % NOTE : solid concentrations are given with respect to total volume (and not with respect to
                        % phase volume, as the liquid concentrations)
                        atomicMassCons{indqp} = atomicMassCons{indqp} - C(indqp, indsp)*cs{indsp};
                    end
                end
            end

            state.atomicMassCons = atomicMassCons;

        end


        function state = updateDischargeSource(model, state)

            vols = model.G.getVolumes();

            Rprec = state.Rprecipitation;

            state.dischargeSource = vols.*Rprec;

        end

        function state = updateDischargeMassCons(model, state)


            accum    = state.dischargeAccum;
            source   = state.dischargeSource;


            % TODO : should not be hard-coded but moved outside
            includeDiffusion = false;

            if includeDiffusion

                indsolid = model.indsolidsp(1);
                op       = model.operators;

                c        = state.cs{indsolid};

                % TODO : Hard coded constant should be moved outside of function
                D = 1e-9;

                useEffectivePorosity = false;

                if useEffectivePorosity
                    ssp      = model.solidPrecipitatePorosity;

                    elyte_vf = state.volumeFraction;
                    solid_vf = state.solidVolumeFraction;
                    % The solid concentration is a total concentration, we convert it into a concentration that exclude
                    % electrode part
                    poro = ssp.*solid_vf + elyte_vf;
                    c = c./poro;

                    Deff = (poro.^1.5).*D;
                    flux = assembleHomogeneousFlux(model, c, Deff);

                else
                    D = D*ones(size(value(c), 1), 1);
                    flux = assembleHomogeneousFlux(model, c, D);

                end

                divTerm = op.Div(flux);

            else

                divTerm = 0;

            end

            state.dischargeMassCons = accum + divTerm - source;

        end

        function state = computeTimeStepIndicator(model, state)

            spp = model.solidPrecipitatePorosity;

            svf = state.solidVolumeFraction;
            vf  = state.volumeFraction;

            th = vf - spp*svf;

            state.indicator = 5.*th.^(1/10);

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
