classdef Electrolyte < BaseModel

    properties


        %% Input parameters

        % Standard parameters
        species % Structure with following fields
                % - chargeNumber : charge number
                % - transferenceNumber : transference number

        density              % the mass density of the material (symbol: rho)
        ionicConductivity    % a function to determine the ionic conductivity of the electrolyte under given conditions (symbol: kappa)
        diffusionCoefficient % a function to determine the diffusion coefficient of a molecule in the electrolyte under given conditions (symbol: D)
        bruggemanCoefficient % the coefficient for determining effective transport parameters in porous media (symbol: beta)

        %
        % Region support for bruggeman coefficient
        %
        useRegionBruggemanCoefficients % Set to true if the electrolye region for each component should get a specific
                                       % Bruggeman coefficient, as given by regionBruggemanCoefficients. Default value is false
        regionBruggemanCoefficients % Bruggeman coefficients for each region, given as a structure with fields
                                    % - NegativeElectrode 
                                    % - PositiveElectrode 
                                    % - Separator 
        regionTags % Tags that identifies each grid cell for the electrolyte (negative and positive electrodes,
                   % separator). It is used only when useRegionBruggemanCoefficients is true. This property is typically setup by the grid constructor.

        thermalConductivity  % Intrinsic Thermal conductivity of the electrolyte
        specificHeatCapacity % Specific Heat capacity of the electrolyte
        
        nominalEthyleneCarbonateConcentration % only used if a SEI layer is included in the model
        
        % Advanced parameters

        volumeFraction
        effectiveThermalConductivity    % (account for volume fraction)
        effectiveVolumetricHeatCapacity % (account for volume fraction and density)

        %%  helper properties

        constants
        computeConductivityFunc         % Function object (see Utilities/FunctionInterface/Function.m)
        computeConductivity             % function handler (can be called directly)
        computeDiffusionCoefficientFunc % Function object (see Utilities/FunctionInterface/Function.m)
        computeDiffusionCoefficient     % function handler (can be called directly)
        use_thermal
        
        regularisationConstant = 0 % can be useed in update of dmudcs to avoid singular value

    end

    methods

        function model = Electrolyte(inputparams)
        % inputparams is instance of ElectrolyteInputParams or a derived class

            model = model@BaseModel();

            fdnames = {'G'                                    , ...
                       'species'                              , ...
                       'density'                              , ...
                       'ionicConductivity'                    , ...
                       'diffusionCoefficient'                 , ...
                       'bruggemanCoefficient'                 , ...
                       'useRegionBruggemanCoefficients'       , ... 
                       'regionBruggemanCoefficients'          , ... 
                       'regionTags'                           , ... 
                       'thermalConductivity'                  , ...
                       'specificHeatCapacity'                 , ...
                       'volumeFraction'                       , ...
                       'effectiveThermalConductivity'         , ...
                       'effectiveVolumetricHeatCapacity'      , ...
                       'nominalEthyleneCarbonateConcentration', ...
                       'use_thermal'};

            model = dispatchParams(model, inputparams, fdnames);

            if model.useRegionBruggemanCoefficients

                nc    = model.G.getNumberOfCells();
                tags  = model.regionTags;
                bvals = model.regionBruggemanCoefficients;

                b = NaN(nc, 1);
                b(tags == 1) = bvals.NegativeElectrode;
                b(tags == 2) = bvals.PositiveElectrode;
                b(tags == 3) = bvals.Separator;

                model.bruggemanCoefficient = b;
            end
            
            [model.computeConductivity, ...
             model.computeConductivityFunc] = setupFunction(inputparams.ionicConductivity);
            [model.computeDiffusionCoefficient, ...
             model.computeDiffusionCoefficientFunc] = setupFunction(inputparams.diffusionCoefficient);

            model.constants = PhysicalConstants();

            if model.use_thermal

                model.G = model.G.setupCellFluxOperators();

                if isempty(model.effectiveThermalConductivity)

                    vf = model.volumeFraction;
                    bg = model.bruggemanCoefficient;
                    model.effectiveThermalConductivity = vf .^ bg .* model.thermalConductivity;

                end

                if isempty(model.effectiveVolumetricHeatCapacity)

                    vf = model.volumeFraction;
                    rho = model.density;
                    model.effectiveVolumetricHeatCapacity = vf .* rho .* model.specificHeatCapacity;

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
            % Diffusion effective coefficient
            varnames{end + 1} = 'D';
            % Derivative of the chemical potential with respect to the concentrations
            varnames{end + 1} = VarName({}, 'dmudcs', 2);
            % ionic effective conductivity
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

        function model = setTPFVgeometry(model, tPFVgeometry)
        % tPFVgeometry should be instance of TwoPointFiniteVolumeGeometry

            model.G.parentGrid.tPFVgeometry = tPFVgeometry;

        end

        function jsonstruct = exportParams(model)

            jsonstruct = exportParams@BaseModel(model);

            fdnames = {'density'                        , ...
                       'ionicConductivity'              , ...
                       'diffusionCoefficient'           , ...
                       'bruggemanCoefficient'           , ...
                       'thermalConductivity'            , ...
                       'specificHeatCapacity'           , ...
                       'effectiveThermalConductivity'   , ...
                       'effectiveVolumetricHeatCapacity', ...
                       'use_thermal'};
            
            for ifd = 1 : numel(fdnames)
                fdname = fdnames{ifd};
                jsonstruct.(fdname) = model.(fdname);
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

            accum    = zeros(model.G.getNumberOfCells(),1);

            chargeCons = assembleConservationEquation(model, flux, bcsource, source, accum);

            state.chargeCons = chargeCons;

        end

        function state = updateFaceCurrent(model, state)

            G = model.G;
            
            nf       = G.getNumberOfFaces();
            intfaces = G.getIntFaces();

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


        function state = updateAccumTerm(model, state, state0, dt)

            cdotcc  = (state.c - state0.c)/dt;

            effectiveVolumes = model.volumeFraction.*model.G.getVolumes();

            state.massAccum  = effectiveVolumes.*cdotcc;

        end

        function state = updateConductivity(model, state)
            
            brcoef = model.bruggemanCoefficient;
            
            c = state.c;
            T = state.T;

            kappa = model.computeConductivity(c, T);
            state.conductivity = kappa.*model.volumeFraction.^brcoef;

        end

        function state = updateDmuDcs(model, state)

            c   = state.c;   % concentration of Li+
            T   = state.T;   % temperature
            phi = state.phi; % potential

            % calculate the concentration derivative of the chemical potential for each species in the electrolyte
            % In the case of a binary electrolyte, we could have simplified those expressions.
            R = model.constants.R;
            r = model.regularisationConstant;
            
            % We consider only binary electrolytes with two components.
            ncomp = 2;
            dmudcs = cell(ncomp, 1);
            for ind = 1 : ncomp
                dmudcs{ind} = R .* T ./ (c + r);
            end

            state.dmudcs = dmudcs;

        end

        function state = updateDiffusionCoefficient(model, state)

            brcoef = model.bruggemanCoefficient;

            c = state.c;
            T = state.T;

            D = model.computeDiffusionCoefficient(c, T);

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
            sp      = model.species;
            volfrac = model.volumeFraction;

            t = sp.transferenceNumber;

            dmudcs       = state.dmudcs;
            phi          = state.phi;
            T            = state.T;
            c            = state.c;
            conductivity = state.conductivity;

            j = assembleFlux(model, phi, conductivity);

            sum_dmudc = dmudcs{1} + dmudcs{2};
            coef = (1/con.F)*(1 - t)*conductivity.*sum_dmudc;
            jchem = assembleFlux(model, c, coef);

            j = j - jchem;

            state.j = j;

        end

        function state = updateMassFlux(model, state)


            sp = model.species;

            t = sp.transferenceNumber;
            z = sp.chargeNumber;

            % We assume that the current and the diffusion coefficient D has been updated when this function is called
            c = state.c;
            j = state.j;
            D = state.D;

            %% 1. Flux from diffusion
            diffFlux = assembleFlux(model, c, D);
            state.diffFlux = diffFlux;

            %% 2. Flux from electrical forces
            F = model.constants.F;
            fluxE = t ./ (z .* F) .* j;

            %% 3. Sum the two flux contributions
            flux = diffFlux + fluxE;

            state.massFlux = flux;

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
