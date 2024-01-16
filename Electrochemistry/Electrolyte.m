classdef Electrolyte < BaseModel

    properties


        %% Input parameters

        % Standard parameters
        sp % Structure with following fields
           % - z : charge number
           % - t : transference number

        density              % the mass density of the material (symbol: rho)
        ionicConductivity    % a function to determine the ionic conductivity of the electrolyte under given conditions (symbol: kappa)
        diffusionCoefficient % a function to determine the diffusion coefficient of a molecule in the electrolyte under given conditions (symbol: D)
        bruggemanCoefficient % the coefficient for determining effective transport parameters in porous media (symbol: beta)
        thermalConductivity  % Intrinsic Thermal conductivity of the electrolyte
        specificHeatCapacity % Specific Heat capacity of the electrolyte

        % Advanced parameters

        volumeFraction
        effectiveThermalConductivity    % (account for volume fraction)
        effectiveVolumetricHeatCapacity % (account for volume fraction and density)

        %%  helper properties

        constants
        compnames
        ncomp
        computeConductivityFunc
        computeDiffusionCoefficientFunc
        use_thermal

    end

    methods

        function model = Electrolyte(inputparams)
        % inputparams is instance of ElectrolyteInputParams or a derived class

            model = model@BaseModel();

            fdnames = {'G'                              , ...
                       'sp'                             , ...
                       'compnames'                      , ...
                       'density'                        , ...
                       'ionicConductivity'              , ...
                       'diffusionCoefficient'           , ...
                       'bruggemanCoefficient'           , ...
                       'thermalConductivity'            , ...
                       'specificHeatCapacity'           , ...
                       'volumeFraction'                 , ...
                       'effectiveThermalConductivity'   , ...
                       'effectiveVolumetricHeatCapacity', ...
                       'use_thermal'};

            model = dispatchParams(model, inputparams, fdnames);

            model.computeConductivityFunc         = str2func(inputparams.ionicConductivity.functionname);
            model.computeDiffusionCoefficientFunc = str2func(inputparams.diffusionCoefficient.functionname);

            model.ncomp = numel(model.compnames);

            model.constants = PhysicalConstants();

            if model.use_thermal

                if isempty(model.effectiveThermalConductivity)

                    bg = model.bruggemanCoefficient;

                    model.effectiveThermalConductivity = (model.volumeFraction).^bg.*model.thermalConductivity;

                end

                if isempty(model.effectiveVolumetricHeatCapacity)

                    model.effectiveVolumetricHeatCapacity = (model.volumeFraction).*model.density.*model.specificHeatCapacity;

                end

            end

        end

        function model = registerVarAndPropfuncNames(model)
            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@BaseModel(model);

            varnames = {};

            % Temperature [K]
            varnames{end + 1} = 'T';
            % Electrical potential [V]
            varnames{end + 1} = 'phi';
            % Current Source [A] - one value per cell
            varnames{end + 1} = 'eSource';
            % Current Source at boundary [A] - one value per cell (will be typically set to zero for non-boundary cells)
            varnames{end + 1} = 'jBcSource';
            % Conductivity [S m^-1]
            varnames{end + 1} = 'conductivity';
            % Current density flux [A] - one value per face
            varnames{end + 1} = 'j';
            % Residual for the charge conservation equation
            varnames{end + 1} = 'chargeCons';
            % Concentration
            varnames{end + 1} = 'c';
            % Source term in the mass conservation equation
            varnames{end + 1} = 'massSource';
            % Flux term in the mass conservation equation
            varnames{end + 1} = 'massFlux';
            % Accumulation term in the mass conservation equation
            varnames{end + 1} = 'massAccum';
            % Residual of the mass conservation equation
            varnames{end + 1} = 'massCons';
            % Diffusion coefficient
            varnames{end + 1} = 'D';
            % Derivative of the chemical potential with respect to the concentrations
            varnames{end + 1} = VarName({}, 'dmudcs', 2);
            % ionic conductivity
            varnames{end + 1} = 'conductivity';
            % diffusion fluzes
            varnames{end + 1} = 'diffFlux';

            model = model.registerVarNames(varnames);

            if model.use_thermal
                varnames = {'jFace', ...
                            'jFaceBc'};
                model = model.registerVarNames(varnames);
            end

            fn = @Electrolyte.updateChargeConservation;
            inputnames = {'j', 'jBcSource', 'eSource'};
            model = model.registerPropFunction({'chargeCons', fn, inputnames});

            fn = @Electrolyte.updateConductivity;
            model = model.registerPropFunction({'conductivity', fn, {'c', 'T'}});

            fn = @Electrolyte.updateDmuDcs;
            model = model.registerPropFunction({VarName({}, 'dmudcs', 2), fn, {'c', 'T'}});

            fn = @Electrolyte.updateDiffusionCoefficient;
            model = model.registerPropFunction({'D', fn, {'c', 'T'}});

            fn = @Electrolyte.updateCurrent;
            inputnames = {'phi'                   , ...
                          VarName({}, 'dmudcs', 2), ...
                          'conductivity'};
            model = model.registerPropFunction({'j', fn, inputnames});

            fn = @Electrolyte.updateMassFlux;
            model = model.registerPropFunction({'massFlux', fn, {'c', 'j', 'D'}});
            model = model.registerPropFunction({'diffFlux', fn, {'c', 'j', 'D'}});

            fn = @Electrolyte.updateAccumTerm;
            fn = {fn, @(propfunction) PropFunction.accumFuncCallSetupFn(propfunction)};
            model = model.registerPropFunction({'massAccum', fn, {'c'}});

            fn = @Electrolyte.updateCurrentBcSource;
            model = model.registerPropFunction({'jBcSource', fn, {}});

            fn = @Electrolyte.updateMassConservation;
            model = model.registerPropFunction({'massCons', fn, {'massFlux', 'massSource', 'massAccum'}});

            if model.use_thermal

                fn = @Electrolyte.updateFaceCurrent;
                inputnames = {'j', 'jFaceBc'};
                model = model.registerPropFunction({'jFace', fn, inputnames});

                fn = @Electrolyte.updateFaceBcCurrent;
                inputnames = {};
                model = model.registerPropFunction({'jFaceBc', fn, inputnames});

            end

        end

        function state = updateMassConservation(model, state)

            accum    = state.massAccum;
            source   = state.massSource;
            flux     = state.massFlux;
            bcsource = 0;

            state.massCons = assembleConservationEquation(model, flux, bcsource, source, accum);

        end

        function state = updateChargeConservation(model, state)

            flux     = state.j;
            bcsource = state.jBcSource;
            source   = state.eSource;

            accum    = zeros(model.G.cells.num,1);

            chargeCons = assembleConservationEquation(model, flux, bcsource, source, accum);

            state.chargeCons = chargeCons;

        end

        function state = updateFaceCurrent(model, state)

            G = model.G;
            nf = G.faces.num;
            intfaces = model.operators.internalConn;

            j       = state.j;
            jFaceBc = state.jFaceBc;

            zeroFaceAD = model.AutoDiffBackend.convertToAD(zeros(nf, 1), j);
            jFace = zeroFaceAD + jFaceBc;
            jFace(intfaces) = j;

            state.jFace = jFace;

        end

        function state = updateFaceBcCurrent(model, state)

            state.jFaceBc = 0;

        end

        function model = setTPFVgeometry(model, tPFVgeometry)
        % tPFVgeometry should be instance of TwoPointFiniteVolumeGeometry or MutableTwoPointFiniteVolumeGeometry

            model.G.parentGrid.tPFVgeometry = tPFVgeometry;

        end

        function state = updateAccumTerm(model, state, state0, dt)

            cdotcc  = (state.c - state0.c)/dt;

            effectiveVolumes = model.volumeFraction.*model.G.getVolumes();

            state.massAccum  = effectiveVolumes.*cdotcc;

        end

        function state = updateConductivity(model, state)

            computeConductivity = model.computeConductivityFunc;

            c = state.c;
            T = state.T;

            state.conductivity = computeConductivity(c, T);

        end

        function state = updateDmuDcs(model, state)

            ncomp = model.ncomp; % number of components

            c   = state.c;   % concentration of Li+
            T   = state.T;   % temperature
            phi = state.phi; % potential

            % calculate the concentration derivative of the chemical potential for each species in the electrolyte
            % In the case of a binary electrolyte, we could have simplified those expressions.
            R = model.constants.R;
            dmudcs = cell(2, 1);
            for ind = 1 : ncomp
                dmudcs{ind} = R .* T ./ c;
            end

            state.dmudcs = dmudcs;

        end

        function state = updateDiffusionCoefficient(model, state)

            brcoef = model.bruggemanCoefficient;

            computeD = model.computeDiffusionCoefficientFunc;

            c = state.c;
            T = state.T;

            D = computeD(c, T);

            % set effective coefficient
            state.D = D .* model.volumeFraction .^ brcoef;

        end

        function state = updateCurrentBcSource(model, state)
        % no boundary current fluxes (only volumetric from the reactions)
            state.jBcSource = 0;
            state.jFaceBc = 0;
        end

        function state  = updateCurrent(model, state)

            bg      = model.bruggemanCoefficient;
            con     = model.constants;
            sp      = model.sp;
            volfrac = model.volumeFraction;

            dmudcs       = state.dmudcs;
            phi          = state.phi;
            T            = state.T;
            c            = state.c;
            conductivity = state.conductivity;

            % Compute effective ionic conductivity in porous media
            conductivityeff = conductivity.*volfrac.^bg;

            state.conductivityeff = conductivityeff;
            j = assembleFlux(model, phi, conductivityeff);

            sum_dmudc = dmudcs{1} + dmudcs{2};
            coef = (1/con.F)*(1 - sp.t(1))*conductivityeff.*sum_dmudc;
            jchem = assembleFlux(model, c, coef);

            j = j - jchem;

            state.j = j;

        end

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

            %% 3. Sum the two flux contributions
            flux = diffFlux + fluxE;

            state.massFlux = flux;

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
