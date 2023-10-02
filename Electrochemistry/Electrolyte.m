classdef Electrolyte < BaseModel

    properties


        %% Standard parameters
        sp % Structure with following fields
           % - z : charge number
           % - t : transference number


        density              % the mass density of the material (symbol: rho)
        ionicConductivity    % a function to determine the ionic conductivity of the electrolyte under given conditions (symbol: kappa)
        diffusionCoefficient % a function to determine the diffusion coefficient of a molecule in the electrolyte under given conditions (symbol: D)        
        bruggemanCoefficient % the coefficient for determining effective transport parameters in porous media (symbol: beta)


        EffectiveThermalConductivity
        EffectiveVolumetricHeatCapacity

        computeConductivityFunc
        computeDiffusionCoefficientFunc

        %% Advanced parameters
        volumeFraction

        % helper properties
        compnames
        ncomp

        use_thermal
        
    end

    methods

        function model = Electrolyte(paramobj)
        % paramobj is instance of ElectrolyteInputParams or a derived class

            model = model@BaseModel();

            fdnames = {'G'                   , ...
                       'sp'                  , ...
                       'compnames'           , ...
                       'density'             , ...
                       'ionicConductivity'   , ...
                       'diffusionCoefficient', ...
                       'bruggemanCoefficient', ...
                       'volumeFraction'      , ...
                       'use_thermal'};

            model = dispatchParams(model, paramobj, fdnames);

            model.computeConductivityFunc         = str2func(paramobj.ionicConductivity.functionname);
            model.computeDiffusionCoefficientFunc = str2func(paramobj.diffusionCoefficient.functionname);

            model.ncomp = numel(model.compnames);

            % We set the electrolyte volumeFraction based on the porosity of the separator
            % G = model.G;
            % Gp = G.mappings.parentGrid;

            % model.volumeFraction = NaN(G.cells.num, 1);

            % elyte_cells = zeros(Gp.cells.num, 1);
            % elyte_cells(G.mappings.cellmap) = (1 : model.G.cells.num)';
            % elyte_cells_sep = elyte_cells(model.(sep).G.mappings.cellmap);
            % model.volumeFraction(elyte_cells_sep) = model.(sep).porosity;

            if model.use_thermal
                % The effective thermal conductivity in the common region between electrode and electrolyte is setup when the battery is
                % set up. Here we set up the thermal conductivity of the electrolyte in the separator region (we assume for
                % the moment constant values for both porosity and thermal conductivity but this can be changed).

                model.EffectiveThermalConductivity = NaN(G.cells.num, 1);
                model.EffectiveThermalConductivity(elyte_cells_sep) = model.(sep).porosity.*model.thermalConductivity;
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


        function state = updateAccumTerm(model, state, state0, dt)

            cdotcc  = (state.c - state0.c)/dt;

            effectiveVolumes = model.volumeFraction.*model.G.cells.volumes;

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

            brcoef = model.BruggemanCoefficient;

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

            % volume fraction of electrolyte
            volfrac = model.volumeFraction;
            % Compute effective ionic conductivity in porous media
            conductivityeff = conductivity.*volfrac.^brcoef;

            state.conductivityeff = conductivityeff;
            j = assembleFlux(model, phi, conductivityeff);

            sum_dmudc = dmudcs{1} + dmudcs{2};
            coef = (1/F)*(1 - sp.t(1))*conductivityeff.*sum_dmudc;
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
