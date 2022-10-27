classdef Electrolyte < ElectroChemicalComponent

    properties

        Separator

        sp
        compnames
        ncomp

        indchargecarrier
        chargeCarrierName


        volumeFraction

        thermalConductivity % intrinsic thermal conductivity value
        heatCapacity % intrinsic heat capacity value
        density % [kg m^-3] (Note : only of the liquid part, the density of the separator is given there)

        EffectiveThermalConductivity
        EffectiveHeatCapacity

        computeConductivityFunc
        computeDiffusionCoefficientFunc

        BruggemanCoefficient

    end

    methods

        function model = Electrolyte(paramobj)
        % paramobj is instance of ElectrolyteInputParams or a derived class

            model = model@ElectroChemicalComponent(paramobj);

            model.Separator = Separator(paramobj.Separator);

            fdnames = {'sp'                 , ...
                       'compnames'          , ...
                       'chargeCarrierName'  , ...
                       'thermalConductivity', ...
                       'heatCapacity'       , ...
                       'density'            , ...
                       'use_thermal'        , ...
                       'BruggemanCoefficient'};

            model = dispatchParams(model, paramobj, fdnames);

            model.computeConductivityFunc = str2func(paramobj.Conductivity.functionname);
            model.computeDiffusionCoefficientFunc = str2func(paramobj.DiffusionCoefficient.functionname);

            model.ncomp = numel(model.compnames);
            [isok, indchargecarrier] = ismember(model.chargeCarrierName, model.compnames);
            assert(isok, 'charge carrier not found in the list of components');
            model.indchargecarrier = indchargecarrier;

            sep = 'Separator';

            % We set the electrolyte volumeFraction based on the porosity of the separator
            G = model.G;
            Gp = G.mappings.parentGrid;

            model.volumeFraction = NaN(G.cells.num, 1);

            elyte_cells = zeros(Gp.cells.num, 1);
            elyte_cells(G.mappings.cellmap) = (1 : model.G.cells.num)';
            elyte_cells_sep = elyte_cells(model.(sep).G.mappings.cellmap);
            model.volumeFraction(elyte_cells_sep) = model.(sep).porosity;

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

            model = registerVarAndPropfuncNames@ElectroChemicalComponent(model);

            varnames = { VarName({}, 'cs', 2)     , ...
                         'D'                      , ...
                         VarName({}, 'dmudcs', 2) , ...
                         'conductivity'           , ...
                         VarName({}, 'jchems', 2) , ...
                         'diffFlux'};
            model = model.registerVarNames(varnames);

            fn = @Electrolyte.updateConcentrations;
            model = model.registerPropFunction({VarName({}, 'cs', 2), fn, {'c'}});

            fn = @Electrolyte.updateConductivity;
            model = model.registerPropFunction({'conductivity', fn, {'c', 'T'}});

            fn = @Electrolyte.updateChemicalCurrent;
            model = model.registerPropFunction({VarName({}, 'dmudcs', 2), fn, {'c', 'T'}});
            model = model.registerPropFunction({VarName({}, 'jchems', 2), fn, {'c', 'T'}});

            fn = @Electrolyte.updateDiffusionCoefficient;
            model = model.registerPropFunction({'D', fn, {'c', 'T'}});

            fn = @Electrolyte.updateCurrent;
            model = model.registerPropFunction({'j', fn, {'phi', VarName({}, 'jchems', 2), 'conductivity'}});

            fn = @Electrolyte.updateMassFlux;
            model = model.registerPropFunction({'massFlux', fn, {'c', 'j', 'D'}});
            model = model.registerPropFunction({'diffFlux', fn, {'c', 'j', 'D'}});

            fn = @Electrolyte.assembleAccumTerm;
            model = model.registerPropFunction({'massAccum', fn, {'c'}});

            fn = @Electrolyte.updateCurrentBcSource;
            model = model.registerPropFunction({'jBcSource', fn, {}});
            
        end

        function state = updateConcentrations(model, state)

            cs = cell(2, 1);
            cs{1} = state.c;
            cs{2} = state.c;

            state.cs = cs;

        end

        function state = updateAccumTerm(model, state, state0, dt)

            cdotcc  = (state.c - state0.c)/dt;

            effectiveVolumes = model.volumeFraction.*model.G.cells.volumes;

            state.massAccum  = effectiveVolumes.*cdotcc;
            
        end

        function state = updateConductivity(model, state)
            
            computeConductivity = model.computeConductivityFunc;

            c = state.cs{1};
            T = state.T;
            
            state.conductivity = computeConductivity(c, T);
            
        end

        function state = updateChemicalCurrent(model, state)

            ncomp = model.ncomp; % number of components

            cLi          = state.cs{1}; % concentration of Li+
            T            = state.T;     % temperature
            phi          = state.phi;   % potential
            conductivity = state.conductivity;   % potential

            cs  = state.cs;


            % calculate the concentration derivative of the chemical potential for each species in the electrolyte
            R = model.constants.R;
            dmudcs = cell(2, 1);
            for ind = 1 : ncomp
                dmudcs{ind} = R .* T ./ cs{ind};
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

            ind = model.indchargecarrier;

            % We assume that the current and the diffusion coefficient D has been updated when this function is called
            c = state.cs{ind};
            j = state.j;
            D = state.D;

            %% 1. Flux from diffusion
            diffFlux = assembleFlux(model, c, D);
            state.diffFlux = diffFlux;

            %% 2. Flux from electrical forces
            F = model.constants.F;
            fluxE = model.sp.t(ind) ./ (model.sp.z(ind) .* F) .* j;

            %% 3. Sum the two flux contributions
            flux = diffFlux + fluxE;

            state.massFlux = flux;

        end

    end
end




%{
Copyright 2021-2022 SINTEF Industry, Sustainable Energy Technology
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
