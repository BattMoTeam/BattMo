classdef SeaWaterBattery < BaseModel

    properties

        con = PhysicalConstants();

        Cathode
        CathodeActiveMaterial
        Anode
        AnodeActiveMaterial
        Electrolyte

        include_precipitation
        
        couplingTerms    % list of coupling terms (each element is instance of couplingTerm class)

        %% Helper properties (setup at instantiaion)
        couplingCellDict % Dictionary with key = {coupling names} and value = {coupling cells} (derived from couplingTerms)

        %% Fixed variables
        % The following properties are used when the corresponding variables are considered as fixed (depends on the
        % implementation).
        T

        indexAnodeChargeCarrier
        indexCathodeChargeCarrier
    
        stochAnodeElectron
        stochCathodeElectron

        dodebugupdateplot = false;
        dodebugplot = false;
        dodebugtext


        %% Parameter variable for Newton API methods
        
        residualMaxValue = 1e8 % upper bound for the residuals. If residual is higher, we cut time step

        newton_verbose                 = true;
        do_log_capping                 = true
        do_solid_concentration_capping = true;
        do_nucleation_capping          = true;
        do_volume_fraction_capping     = true;
        do_concentration_capping       = false;

        minimumElectrolyteVolumeFraction = 1e-2 % apriori minimum value of the electrolyte volume fraction enforced in Newton update
        
        logVariablesIndices % Indices of the log variables, used in Newton step control.

    end

    methods

        function model = SeaWaterBattery(inputparams)

            model = model@BaseModel();

            %% Setup the model using the input parameters
            fdnames = {'G'                   , ...
                       'couplingTerms'       , ...
                       'include_precipitation', ...
                       'T'};
            model = dispatchParams(model, inputparams, fdnames);

            % Setup Anode
            model = model.setupAnode(inputparams.Anode);
            % Setup Anode Active Material            
            model = model.setupAnodeActiveMaterial(inputparams.AnodeActiveMaterial);
            % Setup Electrolyte
            model = model.setupElectrolyte(inputparams);
            % Setup Cathode Active Material                        
            model.CathodeActiveMaterial = HydrogenActiveMaterial(inputparams.CathodeActiveMaterial);
            % Setup Cathode
            model.Cathode = HydrogenElectrode(inputparams.Cathode);

            elyte = 'Electrolyte';
            ctam  = 'CathodeActiveMaterial';
            anam  = 'AnodeActiveMaterial';            
            
            spdict = model.(elyte).spdict;
            
            model.indexAnodeChargeCarrier   = spdict(model.(anam).chargeCarrierName);
            model.indexCathodeChargeCarrier = spdict(model.(ctam).chargeCarrierName);
            model.stochAnodeElectron        = model.(anam).stochElectron;
            model.stochCathodeElectron      = model.(ctam).stochElectron;

            if ~model.include_precipitation
                % we do not cap the properties that are not primary variables in this case
                model.do_solid_concentration_capping = false;
                model.do_nucleation_capping          = false;
                model.do_volume_fraction_capping     = false;
            end
            
            model = model.setupUtilities();

            model.dodebugtext = '';

        end

        function model = setupAnode(model, inputparams)
            model.Anode = SeaWaterElectrode(inputparams);
        end
            
        function model = registerVarAndPropfuncNames(model)
            
            model = registerVarAndPropfuncNames@BaseModel(model);
            
            ct    = 'Cathode';
            ctam  = 'CathodeActiveMaterial';
            an    = 'Anode';
            anam  = 'AnodeActiveMaterial';
            elyte = 'Electrolyte';

            model = model.registerVarNames({'T'});
            
            % Update electrolyte volume fraction
            if model.include_precipitation
                nqp = model.(elyte).nqp;
                fn = @() SeaWaterBattery.updateElectrolyteVolumeFractionEquation;
                inputnames = {{ct, 'volumeFraction'}, ...
                              {elyte, 'solidVolumeFraction'}, ...
                              {elyte, 'volumeFraction'}, ...
                              {an, 'volumeFraction'}};
                model = model.registerPropFunction({{elyte, 'volumeFractionEquation'}, fn, inputnames});
            end
            
            % Dispatch temperature
            fn = @() SeaWaterBattery.dispatchTemperature;
            inputnames = {'T'};
            model = model.registerPropFunction({{ctam, 'T'}, fn, inputnames});
            model = model.registerPropFunction({{elyte, 'T'}, fn, inputnames});
            model = model.registerPropFunction({{anam, 'T'}, fn, inputnames});
            
            % recover some indices from electrolyte
            electrolyte = model.(elyte);
            spdict = electrolyte.spdict;
            nsp = electrolyte.spdict.length();

            % Dispatch values from cathode to cathode catalyser
            fn = @() SeaWaterBattery.updateCathodeActiveMaterialConcentration;
            indcc     = model.indexCathodeChargeCarrier;
            inputnames = {VarName({elyte}, 'cs', nsp, indcc)};
            model = model.registerPropFunction({{ctam, 'cElectrolyte'}, fn, inputnames});
            
            fn = @() SeaWaterBattery.updateCathodeActiveMaterialPotentials;
            inputnames = {{ct, 'phi'}, ...
                          {elyte, 'phi'}};
            model = model.registerPropFunction({{ctam, 'phiElectrode'}, fn, inputnames});
            model = model.registerPropFunction({{ctam, 'phiElectrolyte'}, fn, inputnames});

            % Dispatch values from anode to anode catalyser
            fn = @() SeaWaterBattery.updateAnodeActiveMaterialConcentration;
            inputnames = {VarName({elyte}, 'cs', nsp, spdict('H+'))};
            model = model.registerPropFunction({{anam, 'cElectrolyte'}, fn, inputnames});
            
            fn = @() SeaWaterBattery.updateAnodeActiveMaterialPotentials;
            inputnames = {{an, 'phi'}            , ...
                          {an, 'volumeFraction'} , ...
                          {elyte, 'phi'}};
            model = model.registerPropFunction({{anam, 'phiElectrode'}, fn, inputnames});
            model = model.registerPropFunction({{anam, 'phiElectrolyte'}, fn, inputnames});
            model = model.registerPropFunction({{anam, 'specificSurfaceArea'}, fn, inputnames});            
            
            % Update anode source term
            fn = @() SeaWaterBattery.updateAnodeSource;
            inputnames = {{anam, 'R'}};
            model = model.registerPropFunction({{an, 'sourceTerm'}, fn, inputnames});
            
            % fn = @() SeaWaterBattery.updateAnodeBcSource; 
            % inputnames = {{an, 'conductivity'}, {an, 'phi'}};
            % model = model.registerPropFunction({{an, 'jBcSource'}, fn, inputnames});
            
            % Update source term for all the quasi particles
            % nqp = model.(elyte).nqp;
            % fn = @() SeaWaterBattery.updateQuasiParticleSource;
            % inputnames = {{ctam, 'R'} , ...
            %               {anam, 'R'}   , ...
            %               {elyte, 'Rprecipitation'}};
            % model = model.registerPropFunction({VarName({elyte}, 'qpSrcTerms', nqp), fn, inputnames});
            
            % Update Electrolyte Current source
            fn = @() SeaWaterBattery.updateElectrolyteCurrentSource;
            inputnames = {{ctam, 'R'}, ...
                          {anam, 'R'}};
            model = model.registerPropFunction({{elyte, 'eSource'}, fn, inputnames});

            fn = @() SeaWaterBattery.updateCathodeSource;
            model = model.registerPropFunction({{ct, 'sourceTerm'}, fn, {}});

            
            % We remove the unused temperature variable in an and ct
            model = model.removeVarNames({{an, 'T'},
                                          {ct, 'T'},
                                         });

            model = model.setAsStaticVarName('T');

            if ~model.include_precipitation
                melyte = model.(elyte);
                model = model.setAsStaticVarNames({{an, 'volumeFraction'},
                                                   {ct, 'volumeFraction'},
                                                   {elyte, 'volumeFraction'},
                                                   VarName({elyte}, 'cs', melyte.indsolidsp, melyte.nsp)});
                varnames = {'accumTerm', 'sourceTerm', 'massCons'};
                eldes = {an, ct};
                for ielde = 1 : numel(eldes)
                    elde = eldes{ielde};
                    for ivarname = 1 : numel(varnames)
                        varname = varnames{ivarname};
                        model = model.removeVarName({elde, varname});
                    end
                end
                    
            end
            
        end
        
        
        function cleanState = addStaticVariables(model, cleanState, state)

            cleanState = addStaticVariables@BaseModel(model, cleanState, state);
            cleanState.time = state.time;
            cleanState.T    = state.T;

            if ~model.include_precipitation

                ct    = 'Cathode';
                an    = 'Anode';
                elyte = 'Electrolyte';

                indsol = model.(elyte).indsolidsp;
                
                cleanState.(ct).volumeFraction    = state.(ct).volumeFraction;
                cleanState.(an).volumeFraction    = state.(an).volumeFraction;
                cleanState.(elyte).volumeFraction = state.(elyte).volumeFraction;
                cleanState.(elyte).cs{indsol}     = state.(elyte).cs{indsol};
                
            end
            
        end
        
        function model = setupUtilities(model)
        % setup utility model properties

            coupdict = containers.Map();
            coupterms = model.couplingTerms;
            for ind = 1 : numel(coupterms)
                coupterm = coupterms{ind};
                coupdict(coupterm.name) = coupterm.couplingcells;
            end
            model.couplingCellDict = coupdict;
        end


        %% Assembly functions
        
        function state = initializeTemperature(model, state)
            state.T = model.T*ones(model.G.getNumberOfCells(), 1);
        end


        function state = updateElectrolyteVolumeFraction(model, state)
        % Update electrolyte volume fraction
        % only used to compute initial value

            coupDict = model.couplingCellDict;

            if model.include_precipitation
                spp = model.Electrolyte.solidPrecipitatePorosity;
            else
                spp = 1;
            end

            catheps  = state.Cathode.volumeFraction;
            aneps    = state.Anode.volumeFraction;
            solideps = state.Electrolyte.solidVolumeFraction;

            % remove volume fraction occupied discharge product
            vf = 1 - (1 - spp).*solideps;

            % remove volume fraction occupied by anode
            coupcells = coupDict('Anode-Electrolyte');
            vf(coupcells(:, 2)) = vf(coupcells(:, 2)) - aneps(coupcells(:, 1));

            % remove volume fraction occupied by cathode
            coupcells = coupDict('Cathode-Electrolyte');
            vf(coupcells(:, 2)) = vf(coupcells(:, 2)) - catheps(coupcells(:, 1));

            state.Electrolyte.volumeFraction = vf;
            
        end

        function state = updateElectrolyteVolumeFractionEquation(model, state)
        % Update electrolyte volume fraction

            coupDict = model.couplingCellDict;
            spp      = model.Electrolyte.solidPrecipitatePorosity;

            catheps  = state.Cathode.volumeFraction;
            aneps    = state.Anode.volumeFraction;
            solideps = state.Electrolyte.solidVolumeFraction;
            vf       = state.Electrolyte.volumeFraction;
            
            % remove volume fraction occupied discharge product
            computed_vf = 1 - (1 - spp).*solideps;
            
            if ~isa(computed_vf, 'ADI')
                computed_vf = model.AutoDiffBackend.convertToAD(computed_vf, aneps);
            end

            % remove volume fraction occupied by anode
            coupcells = coupDict('Anode-Electrolyte');
            computed_vf(coupcells(:, 2)) = computed_vf(coupcells(:, 2)) - aneps(coupcells(:, 1));

            % remove volume fraction occupied by cathode
            coupcells = coupDict('Cathode-Electrolyte');
            computed_vf(coupcells(:, 2)) = computed_vf(coupcells(:, 2)) - catheps(coupcells(:, 1));
            
            state.Electrolyte.volumeFractionEquation = vf - computed_vf;
            
        end

        function state = dispatchTemperature(model, state)
        % Dispatch temperature
            elyte   = 'Electrolyte';
            cath    = 'Cathode';
            cathcat = 'CathodeActiveMaterial';
            an      = 'Anode';
            ancat   = 'AnodeActiveMaterial';

            T = state.T;

            state.(cathcat).T = T(model.(cath).G.mappings.cellmap);
            state.(ancat).T   = T(model.(an).G.mappings.cellmap);
            state.(elyte).T   = T(model.(elyte).G.mappings.cellmap);
        end

        function state = updateCathodeActiveMaterialConcentration(model, state)
        % Dispatch concentration from anode to anode active material
            
            elyte = 'Electrolyte';
            ct    = 'Cathode';
            ctam  = 'CathodeActiveMaterial';

            spdict    = model.(elyte).spdict;
            coupDict  = model.couplingCellDict;
            adbackend = model.AutoDiffBackend;
            indcc     = model.indexCathodeChargeCarrier;
            
            c = state.(elyte).cs{indcc};

            cElectrolyte = zeros(model.(ct).G.getNumberOfCells(), 1);
            cElectrolyte = adbackend.convertToAD(cElectrolyte, c);

            coupcells = coupDict('Cathode-Electrolyte');
            cElectrolyte(coupcells(:, 1)) = c(coupcells(:, 2));

            state.(ctam).cElectrolyte = cElectrolyte;
            
        end

        function state = updateCathodeActiveMaterialPotentials(model, state)
        % Dispatch potentials values from cathode to cathode active material
            elyte = 'Electrolyte';
            ct    = 'Cathode';
            ctam = 'CathodeActiveMaterial';

            spdict = model.(elyte).spdict;
            coupDict = model.couplingCellDict;

            ctphi    = state.(ct).phi;
            elytephi = state.(elyte).phi;

            phiElectrode = ctphi;

            % Initialize values obtained from electrolyte
            phiElectrolyte = 0*ctphi;

            coupcells = coupDict('Cathode-Electrolyte');
            phiElectrolyte(coupcells(:, 1)) = elytephi(coupcells(:, 2));

            state.(ctam).phiElectrode = phiElectrode;
            state.(ctam).phiElectrolyte = phiElectrolyte;
        end

        function state = updateAnodeActiveMaterialConcentration(model, state)
        % Dispatch concentration from anode to anode active material
            
            elyte = 'Electrolyte';
            an    = 'Anode';
            ancat = 'AnodeActiveMaterial';

            spdict    = model.(elyte).spdict;
            coupDict  = model.couplingCellDict;
            adbackend = model.AutoDiffBackend;
            indcc     = model.indexAnodeChargeCarrier;
            
            c = state.(elyte).cs{indcc};

            cElectrolyte = zeros(model.(an).G.getNumberOfCells(), 1);
            cElectrolyte = adbackend.convertToAD(cElectrolyte, c);

            coupcells = coupDict('Anode-Electrolyte');
            cElectrolyte(coupcells(:, 1)) = c(coupcells(:, 2));

            state.(ancat).cElectrolyte = cElectrolyte;
        end

        function state = updateAnodeActiveMaterialPotentials(model, state)
        % Dispatch potentials values from anode to anode active material
            elyte = 'Electrolyte';
            an    = 'Anode';
            ancat = 'AnodeActiveMaterial';

            spdict    = model.(elyte).spdict;
            coupDict  = model.couplingCellDict;
            nc        = model.(an).G.getNumberOfCells();
            adbackend = model.(an).AutoDiffBackend;
            
            anphi    = state.(an).phi;
            elytephi = state.(elyte).phi;
            vf       = state.(an).volumeFraction;
            
            phiElectrode = anphi;

            % Initialize values obtained from electrolyte
            adsample = getSampleAD(anphi, elytephi);
            phiElectrolyte = adbackend.convertToAD(zeros(nc, 1), adsample);

            coupcells = coupDict('Anode-Electrolyte');
            phiElectrolyte(coupcells(:, 1)) = elytephi(coupcells(:, 2));

            specificSurfaceArea = 6./(100e-6).*(0.5.*(1 - cos(2*pi*vf)));
            specificSurfaceArea(value(vf) > 1 | (value(vf) < 0)) = 0;
            
            state.(ancat).phiElectrode   = phiElectrode;
            state.(ancat).phiElectrolyte = phiElectrolyte;
            state.(ancat).specificSurfaceArea = specificSurfaceArea;
        end

        function state = updateAnodeSource(model, state)
        % Update anode source term
            an    = 'Anode';
            elyte = 'Electrolyte';

            V    = model.Anode.V;
            vols = model.Anode.G.getVolumes();
            F    = model.con.F;
            z    = model.stochAnodeElectron;
            
            R = state.AnodeActiveMaterial.R;

            state.Anode.sourceTerm = -R.*vols.*V;
            
        end

        function state = updateCathodeSource(model, state)
        % Update cathode source term

            % V = model.Cathode.V;
            % vols = model.Anode.G.getVolumes();
            % R = state.CathodeActiveMaterial.R;
            % state.Cathode.sourceTerm = -R.*vols.*V;
                
            state.Cathode.sourceTerm = 0;
            
        end

        function state = updateElectrolyteCurrentSource(model, state)
        % Update Electrolyte Current source
            elyte = 'Electrolyte';

            vols     = model.(elyte).G.getVolumes();
            zan      = model.stochAnodeElectron;
            zct      = model.stochCathodeElectron;
            coupDict = model.couplingCellDict;
            F        = model.con.F;
            
            cathR = state.CathodeActiveMaterial.R;
            anR   = state.AnodeActiveMaterial.R;
            phi   = state.(elyte).phi;

            % initialization
            eSource = 0*phi;

            coupcells = coupDict('Cathode-Electrolyte');
            eSource(coupcells(:, 2)) = zct*F*vols(coupcells(:, 2)).*cathR(coupcells(:, 1));

            coupcells = coupDict('Anode-Electrolyte');
            eSource(coupcells(:, 2)) = zan*F*vols(coupcells(:, 2)).*anR(coupcells(:, 1));

            state.Electrolyte.eSource = eSource;
        end


        %% Specialized Newton API methods
        
        function [state, report] = stepFunction(model, state, state0, dt, drivingForces, linsolver, nonlinsolver, iteration, varargin)
        % We stop simulation  if residuals become to high
            
            [state, report] = stepFunction@BaseModel(model, state, state0, dt, drivingForces, linsolver, nonlinsolver, iteration, varargin{:});

            if iteration > 1
                if any(report.Residuals >= model.residualMaxValue)
                    report.Failure = true;
                    report.FailureMsg = 'Too large residual';
                end
            end
            
        end

        function state = addVariables(model, state)
        % Given a state where only the primary variables are defined, this
        % functions add all the additional variables that are computed in the assembly process and have some physical
        % interpretation.
        %
        % To do so, we use getEquations function and sends dummy variable for state0, dt and drivingForces

        % Values that need to be set to get the function getEquations running

            dt = 1;
            state0 = state;
            drivingForces = model.getValidDrivingForces();
            drivingForces.src = @(time) 0;

            % We call getEquations to update state

            [~, state] = model.getEquations(state0, state, dt, drivingForces, 'ResOnly', true);

            % We remove the variables that are not meaningfull in this context, to avoid any confusion
            ct    = 'Cathode';
            ctam  = 'CathodeActiveMaterial';
            an    = 'Anode';
            anam  = 'AnodeActiveMaterial';
            elyte = 'Electrolyte';

            varnames = {{elyte, 'nucleationEquation'}    , ...
                        {elyte, 'volumeFractionEquation'}, ...
                        {elyte, 'dischargeMassCons'}     , ...
                        {elyte, 'chargeCons'}            , ...
                        {elyte, 'qpMassCons', 1}         , ...
                        {elyte, 'qpMassCons', 2}         , ...
                        {elyte, 'qpMassCons', 3}         , ...
                        {elyte, 'atomicMassCons', 1}     , ...
                        {elyte, 'atomicMassCons', 2}     , ...
                        {elyte, 'atomicMassCons', 3}     , ...
                        {an, 'massCons'}                 , ...
                        {an, 'galvanostatic'}            , ...
                        {ct, 'massCons'}                 , ...
                        {ct, 'galvanostatic'}};


            for ivar = 1 : numel(varnames)
                varname = varnames{ivar};
                state = model.setNewProp(state, varname, []);
            end
            
        end

        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)

            [state, report] = updateAfterConvergence@BaseModel(model, state0, state, dt, drivingForces);

            state = model.addVariables(state);
            
        end

        function validforces = getValidDrivingForces(model)
            
            validforces = struct('src', [], 'stopFunction', []);
            
        end

        function model = validateModel(model, varargin)

            if isempty(model.computationalGraph)
                model = model.setupComputationalGraph();
            end

            cg = model.computationalGraph;
            model.primaryVarNames = cg.getPrimaryVariableNames();
            model.funcCallList = cg.getOrderedFunctionCallList();

            primaryvarnames = model.getPrimaryVariableNames();

            logVariablesIndices = [];
            for ip = 1 : numel(primaryvarnames)
                p = primaryvarnames{ip};
                if strcmp(p{1}, 'Electrolyte')
                    if strcmp(p{2}, 'pcs')
                        logVariablesIndices(end + 1) = ip;
                    end
                end
            end

            model.logVariablesIndices = logVariablesIndices;

        end

        function [state, report] = updateState(model, state, problem, dx, drivingForces)

            
            elyte = 'Electrolyte';
            anam  = 'AnodeActiveMaterial';
            ctam  = 'CathodeActiveMaterial';
            ct    = 'Cathode';
            an    = 'Anode';
            
            verbose = model.newton_verbose;

            if model.do_log_capping

                
                %%  modify udpatees for log variables
                elyte = 'Electrolyte';
                loginds = model.logVariablesIndices;
                maxdx = 0;
                for ilog = 1 : numel(loginds)
                    logind = loginds(ilog);
                    maxdx = max(maxdx, max(abs(dx{logind})));
                end
                if (maxdx > 10) & (maxdx < 15)
                    factor = 10/maxdx;
                    for ind = 1 : numel(dx)
                        dx{ind} = factor*dx{ind};
                    end
                    if verbose
                        fprintf('Newton step of max log concentration size %g has been reduced with factor %g\n', maxdx, factor);
                    end
                elseif (maxdx > 15)
                    report.Failure = true;
                    report.FailureMsg = 'Too large updates in the logarithm variables';
                    return
                end


            end
            
            [state, report] = updateState@BaseModel(model, state, problem, dx, drivingForces);

            if model.do_solid_concentration_capping
                
                indsolid = model.(elyte).indsolidsp(1);
                
                if min(state.(elyte).cs{indsolid}) < 0
                    if verbose
                        fprintf('\n** cap solid concentration, min value %g\n\n', min(state.(elyte).cs{indsolid}));
                    end
                    state.(elyte).cs{indsolid} = max(0, state.(elyte).cs{indsolid});
                end
                
            end
            
            if model.do_nucleation_capping
                
                if min(state.(elyte).nucleation) < 0
                    if verbose
                        fprintf('\n** cap nucleation, min value %g\n\n', min(state.(elyte).nucleation));
                    end
                    state.(elyte).nucleation = max(0, state.(elyte).nucleation);
                end
                
            end
            
            if model.do_volume_fraction_capping & model.do_solid_concentration_capping
                
                
                if min(state.(an).volumeFraction) < 0
                    if verbose
                        fprintf('\n** cap anode volume fraction, min value %g\n\n', min(state.(an).volumeFraction));
                    end
                    state.(an).volumeFraction = max(0, state.(an).volumeFraction);
                end
                
                if min(state.(ct).volumeFraction) < 0
                    if verbose
                        fprintf('\n** cap cathode volume fraction, min value %g\n\n', min(state.(ct).volumeFraction));
                    end
                    state.(ct).volumeFraction = max(0, state.(ct).volumeFraction);
                end
                
                min_elyte_vf = model.minimumElectrolyteVolumeFraction;
                
                if min(state.(elyte).volumeFraction) < min_elyte_vf
                    if verbose
                        fprintf('\n** electrolyte volume fraction, min value %g\n\n', min(state.(elyte).volumeFraction));
                    end
                    state.(elyte).volumeFraction = max(min_elyte_vf, state.(elyte).volumeFraction);
                end
                
            end
            
            if model.do_concentration_capping
                
                for ind = 1 : model.(elyte).nlogsp;
                    state.(elyte).pcs{ind} = max(-30, state.(elyte).pcs{ind});
                    state.(elyte).pcs{ind} = min(30, state.(elyte).pcs{ind});
                end
                
            end


            
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
