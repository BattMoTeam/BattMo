classdef IonomerMembrane < ElectronicComponent

    properties

        volumeFraction

        H2O % with fields
        % H2O.c0 : Reference concentration (value is overwritten at initialization in current implementation)
        % H2O.D  : diffusion coefficient for water
        % H2O.V0 : Molar mass (needed for function groupHydration which is only needed in setup of initial condition and not for assembly)

        OH % with fields
        % OH.xi : OH occupation
        % OH.z  : charge number
        % OH.t  : transference number

        cT % Total concentration of charged groups (one value per cell)

        V % molar volume (needed for function groupHydration which is only needed in setup of initial condition and not
          % for assembly, and also for activity computation)

        tortuosity

    end

    methods

        function model = IonomerMembrane(inputparams)

            inputparams.use_thermal = false;
            model = model@ElectronicComponent(inputparams);

            fdnames = {'volumeFraction', ...
                       'H2O'           , ...
                       'OH'            , ...
                       'V'             , ...
                       'tortuosity'};
            model = dispatchParams(model, inputparams, fdnames);

            cT = inputparams.cT;
            nc = model.G.getNumberOfCells();
            model.cT = cT*ones(nc, 1);

            model.constants = PhysicalConstants();

        end


        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@ElectronicComponent(model);

            varnames = {};

            % Water concentration (per total volume)
            varnames{end + 1} = 'H2Oceps';
            % Water concentration
            varnames{end + 1} = 'H2Oc';
            % OH concentration
            varnames{end + 1} = 'cOH';

            %% Component properties

            % H2O activity
            varnames{end + 1} = 'H2Oa';
            % Chemical potential derivative
            varnames{end + 1} = 'H2Odmudc';

            %% Flux variables

            % Chemical flux
            % TODO : find better name
            varnames{end + 1} = 'jchem';
            % H2O diffusion flux
            varnames{end + 1} = 'H2OdiffFlux';
            % H2O migration flux
            varnames{end + 1} = 'H2OmigFlux';

            %% Source terms

            % H2O source [mol/s] (one value per cell)
            varnames{end + 1} = 'H2OSource';
            % OH source [mol/s] (one value per cell)
            varnames{end + 1} = 'OHsource';

            %% Accumulation term

            % H2O accumulation term
            varnames{end + 1} = 'H2Oaccum';

            %% Mass conservation equation

            % H2O mass conservation equation
            varnames{end + 1} = 'H2OmassCons';

            %% Activity equation
            varnames{end + 1} = 'activityEquation';

            model = model.registerVarNames(varnames);

            % update water concentration
            fn = @() IonomerMembrane.updateH2Oc;
            inputnames = {'H2Oceps'};
            model = model.registerPropFunction({'H2Oc', fn, inputnames});

            % update water activity equation
            fn = @() IonomerMembrane.updateActivityEquation;
            inputnames = {'H2Oc', 'H2Oa', 'T'};
            model = model.registerPropFunction({'activityEquation', fn, inputnames});

            % update effective conductivity
            fn = @() IonomerMembrane.updateConductivity;
            inputnames = {'T', 'H2Oa'};
            model = model.registerPropFunction({'conductivity', fn, inputnames});

            % update chemical potential derivative
            fn = @() IonomerMembrane.updatedmudc;
            inputnames = {'T', 'H2Oc'};
            model = model.registerPropFunction({'H2Odmudc', fn, inputnames});

            % update chemical flux
            fn = @() IonomerMembrane.updatejchem;
            inputnames = {'conductivity', 'H2Odmudc', 'H2Oc'};
            model = model.registerPropFunction({'jchem', fn, inputnames});

            % update electronic flux
            fn = @() IonomerMembrane.updateCurrent;
            inputnames = {'conductivity', 'phi', 'jchem'};
            model = model.registerPropFunction({'j', fn, inputnames});

            % update H2O diffusion flux
            fn = @() IonomerMembrane.updateH2OdiffFlux;
            inputnames = {'H2Oc'};
            model = model.registerPropFunction({'H2OdiffFlux', fn, inputnames});

            % update H2O migration flux
            fn = @() IonomerMembrane.updateH2OmigFlux;
            inputnames = {'j'};
            model = model.registerPropFunction({'H2OmigFlux', fn, inputnames});

            % update electronic source
            fn = @() IonomerMembrane.updateEsource;
            inputnames = {'OHsource'};
            model = model.registerPropFunction({'eSource', fn, inputnames});

            % update electronic source
            fn = @() IonomerMembrane.assembleH2Oaccum;
            functionCallSetupFn = @(propfunction) PropFunction.accumFuncCallSetupFn(propfunction);
            fn = {fn, functionCallSetupFn};
            inputnames = {'H2Oceps'};
            model = model.registerPropFunction({'H2Oaccum', fn, inputnames});

            % assemble H2O mass conservation equation
            fn = @() IonomerMembrane.assembleH2OmassCons;
            inputnames = {'H2OdiffFlux', 'H2OmigFlux', 'H2Oaccum', 'H2OSource'};
            model = model.registerPropFunction({'H2OmassCons', fn, inputnames});

            % the OH concentration is constant
            fn = @() IonomerMembrane.setupOHconcentration;
            inputnames = {};
            model = model.registerPropFunction({'cOH', fn, inputnames});

        end

        function state = updateH2Oc(model, state)

            H2Oceps = state.H2Oceps;

            vf = model.volumeFraction;

            state.H2Oc = H2Oceps./vf;

        end


        function state = setupOHconcentration(model, state)

            state.cOH = model.cT;

        end

        function state = updateActivityEquation(model, state)

            V  = model.V;
            V0 = model.H2O.V0;

            c  = state.H2Oc;
            aw = state.H2Oa;
            T  = state.T;

            eq = c*V./(1 - c*V0) - ((-0.6.*aw.^3 + 0.85.*aw.^2 - 0.2.*aw + 0.153) .* (T - 313) ...
                                    + 39.*aw.^3 - 47.7.*aw.^2 + 23.4.*aw + 0.117).*29./18.877;

            state.activityEquation = eq;

        end


        function state = updateConductivity(model, state)

            vf  = model.volumeFraction;
            tau = model.tortuosity;

            T = state.T;
            a = state.H2Oa;

            dosimplified = true;

            if ~dosimplified

                a(a > 1) = 1;

                kappa = (0.1334                       ...
                         - 3.882e-4.*T                ...
                         + (0.01148.*T - 3.909).*a    ...
                         - (0.06690.*T - 23.01).*a.^2 ...
                         + (0.1227.*T - 42.61).*a.^3  ...
                         - (0.06021.*T - 21.80).*a.^4) .* 20;

            else

                kappa = 10 + 0*a;

            end

            state.conductivity = kappa.*vf.^tau;

        end


        function state = updatedmudc(model, state)

            T = state.T;
            c = state.H2Oc;

            R = model.constants.R;

            state.H2Odmudc = R.*T./c;

        end


        function state = updatejchem(model, state)

            xi  = model.OH.xi;
            z   = model.OH.z;
            F   = model.constants.F;

            kappaeff = state.conductivity;
            dmudc    = state.H2Odmudc;
            cH2O     = state.H2Oc;

            fluxcoef = dmudc.*kappaeff.*xi./(z.*F);

            state.jchem = assembleFlux(model, cH2O, fluxcoef);

        end


        function state = updateCurrent(model, state)

            kappaeff = state.conductivity;
            phi      = state.phi;
            jchem    = state.jchem;

            jelec = assembleFlux(model, phi, kappaeff);

            state.j = jelec + jchem;

        end

        function state = updateH2OdiffFlux(model, state)

            D   = model.H2O.D;
            vf  = model.volumeFraction;
            tau = model.tortuosity;

            c = state.H2Oc;

            Deff = D.*vf.^tau;

            state.H2OdiffFlux = assembleFlux(model, c, Deff);

        end


        function state = updateH2OmigFlux(model, state)

            j = state.j;

            xi = model.OH.xi;
            z  = model.OH.z;
            t  = model.OH.t;
            F  = model.constants.F;

            %% todo : check expression, not same as in paper
            state.H2OmigFlux = xi.*t/(z.*F).*j;

        end

        function state = assembleH2OmassCons(model, state)

            diffFlux = state.H2OdiffFlux;
            migFlux  = state.H2OmigFlux;
            accum    = state.H2Oaccum;
            source   = state.H2OSource;

            flux = diffFlux + migFlux;
            bcsource  = 0;

            state.H2OmassCons = assembleConservationEquation(model, flux, bcsource, source, accum);

        end

        function state = updateEsource(model, state)

            F = model.constants.F;
            z = model.OH.z;

            state.eSource = z*F*state.OHsource;

        end

        function state = assembleH2Oaccum(model, state, state0, dt)

            ceps  = state.H2Oceps;
            ceps0 = state0.H2Oceps;

            vols = model.G.getVolumes();

            state.H2Oaccum = 1/dt*vols.*(ceps - ceps0);

        end



    end

    methods(Static)

        function cH2O = groupHydration(model, aw, T)

            V = model.V;
            V0 = model.H2O.V0;

            aw(aw > 1) = 1;
            %   lambda corresponds to the number of water molecules per charged
            %   group, in the membrane and ionomer phase, given the membrane water
            %   activity and the temperature. The default polynomial relationship is based on
            %   Jiao et al., Int J Hydrogen Energy 2014, doi:10.1016/j.ijhydene.2014.01.180
            %   using only the water activity < 1 portion of their equation.
            lambda = ((-0.6.*aw.^3 + 0.85.*aw.^2 - 0.2.*aw + 0.153) .* (T - 313) ...
                      + 39.*aw.^3 - 47.7.*aw.^2 + 23.4.*aw + 0.117).*29./18.877;

            cH2O = lambda./(V0.*lambda + V);

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
