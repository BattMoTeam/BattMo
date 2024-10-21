classdef Electrolyser < BaseModel

    properties

        con = PhysicalConstants();

        % Components
        IonomerMembrane
        OxygenEvolutionElectrode
        HydrogenEvolutionElectrode

        couplingTerms
        couplingNames


    end

    methods

        function model = Electrolyser(inputparams)

            model = model@BaseModel();

            inputparams = inputparams.validateInputParams();
            
            fdnames = {'G' , ...
                       'couplingTerms'};
            model = dispatchParams(model, inputparams, fdnames);

            model.OxygenEvolutionElectrode   = EvolutionElectrode(inputparams.OxygenEvolutionElectrode);
            model.HydrogenEvolutionElectrode = EvolutionElectrode(inputparams.HydrogenEvolutionElectrode);
            model.IonomerMembrane            = IonomerMembrane(inputparams.IonomerMembrane);

            % setup couplingNames
            model.couplingNames = cellfun(@(x) x.name, model.couplingTerms, 'uniformoutput', false);

            model = model.setupForSimulation();
            
        end

        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@BaseModel(model);

            varnames = {'T', ...
                        VarName({}, 'controlEqs', 2)};

            model = model.registerVarNames(varnames);

            inm = 'IonomerMembrane';
            her = 'HydrogenEvolutionElectrode';
            oer = 'OxygenEvolutionElectrode';

            ctl = 'CatalystLayer';
            exr = 'ExchangeReaction';
            ptl = 'PorousTransportLayer';

            fn = @Electrolyser.dispatchTemperature;
            inputvarnames = {'T'};
            model = model.registerPropFunction({{her, 'T'}, fn, inputvarnames});
            model = model.registerPropFunction({{oer, 'T'}, fn, inputvarnames});
            model = model.registerPropFunction({{inm, 'T'}, fn, inputvarnames});

            fn = @Electrolyser.updateC;

            eldes = {her, oer};
            layers = {ctl, exr};

            fn = @Electrolyser.updateIonomerSources;
            inputvarnames = {};
            for ielde = 1 : numel(eldes)
                elde = eldes{ielde};
                inputvarnames{end + 1} = {elde, exr, 'H2OexchangeRate'};
                inputvarnames{end + 1} = {elde, ctl, 'inmrH2Osource'};
                inputvarnames{end + 1} = {elde, ctl, 'inmrOHsource'};
                inputvarnames{end + 1} = {elde, exr, 'OHexchangeRate'};
            end
            model = model.registerPropFunction({{inm, 'H2OSource'}, fn, inputvarnames});
            model = model.registerPropFunction({{inm, 'OHsource'}, fn, inputvarnames});


            fn = @Electrolyser.dispatchIonomerToReactionLayers;
            inputvarnames = {{inm, 'phi'}, ...
                             {inm, 'cOH'}, ...
                             {inm, 'H2Oa'}};
            for ielde = 1 : numel(eldes)
                for ilayer = 1 : numel(layers)
                    elde = eldes{ielde};
                    layer = layers{ilayer};
                    model = model.registerPropFunction({{elde, layer, 'cOHinmr'}, fn, inputvarnames});
                    model = model.registerPropFunction({{elde, layer, 'phiInmr'}, fn, inputvarnames});
                    model = model.registerPropFunction({{elde, layer, 'H2OaInmr'}, fn, inputvarnames});
                end
            end

            fn = @Electrolyser.setupControl;
            fn = {fn, @PropFunction.drivingForceFuncCallSetupFn};
            inputvarnames = {{oer, ctl, 'I'}, ...
                             {her, ctl, 'E'}
                            };
            model = model.registerPropFunction({VarName({}, 'controlEqs', 2), fn, inputvarnames});

            model = model.setAsStaticVarNames({{inm, 'jBcSource'}     , ...
                                                  {oer, ptl, 'jBcSource'}, ...
                                                  {her, ptl, 'jBcSource'}, ...
                                                  'T'});

            model = model.removeVarNames({{her, ctl, 'I'}, ...
                                          {her, ctl, 'eSource'}});

            
        end


        function [model, state] = setupBcAndInitialState(model)

            oer = 'OxygenEvolutionElectrode';
            her = 'HydrogenEvolutionElectrode';
            ptl = 'PorousTransportLayer';
            ctl = 'CatalystLayer';
            inm = 'IonomerMembrane';

            con = model.con;

            pGas = 101325*Pascal; % pressure of active gas (O2 or H2)
            cOH  = 1000*mol/(meter^3);
            T    = 333.15*Kelvin;
            R    = con.R;
            F    = con.F;

            state.(inm) = model.(inm).setupOHconcentration();

            nc = model.G.getNumberOfCells();
            state.T = T*ones(nc, 1);

            state = model.dispatchTemperature(state);

            eldes = {her, oer};

            for ielde = 1 : numel(eldes)

                elde = eldes{ielde};

                % Compute density and partial molar volume using dedicated functions
                % which are only used at initialisation. The partial molar volumes are kept constant for the remaining of
                % the simulation (at least in current implementation!).
                liqrho = model.(elde).(ptl).density(cOH, T);
                Vs     = model.(elde).(ptl).partialMolarVolume(cOH, liqrho, T);

                nc = model.(elde).(ptl).G.getNumberOfCells();

                fun = @(s) leverett(model.(elde).(ptl).leverettCoefficients, s); % Define Leverett function handle
                sLiquid = fzero(fun, 0.7); % Solve equilibrium liquid saturation

                svf = model.(elde).(ptl).solidVolumeFraction;
                lvf = sLiquid.*(1 - svf);

                state.(elde).(ptl).liqeps    = sLiquid.*(1 - svf);
                state.(elde).(ptl).OHceps    = cOH.*lvf;
                state.(elde).(ptl).liqrhoeps = liqrho*lvf;

                % OBS : the following two values are not the ones that are finally assigned (only sent here so that dispatch functions can be used)
                state.(elde).(ptl).phi = NaN*zeros(nc, 1);
                state.(elde).(ptl).E   = NaN;

                state = model.evalVarName(state, {elde, ptl, 'vaporPressure'});

                H2Ovp = state.(elde).(ptl).vaporPressure;
                H2Ogrho = H2Ovp*model.(elde).(ptl).sp.H2O.molecularWeight./(con.R*T);
                gvf = (1 - lvf - model.(elde).(ptl).solidVolumeFraction); % Gas volume fraction

                state.(elde).(ptl).H2Ogasrhoeps = H2Ogrho.*gvf;

                switch elde

                  case her

                    H2p = pGas;
                    H2rho = H2p.*model.(elde).(ptl).sp.H2.molecularWeight / (con.R * T);
                    state.(elde).(ptl).H2rhoeps = H2rho*gvf;

                  case oer

                    O2p = pGas;
                    O2rho = O2p.*model.(elde).(ptl).sp.O2.molecularWeight / (con.R * T);
                    state.(elde).(ptl).O2rhoeps = O2rho*gvf;

                  otherwise

                    error('electrode not recognized');

                end

                state = model.evalVarName(state, {elde, ctl, 'Eelyte'});

                model.(elde).(ptl).Vs = Vs;

            end

            nc = model.(inm).G.getNumberOfCells();

            % We use water activity in oer to setup activity
            nc_inm = model.(inm).G.getNumberOfCells();

            aw = state.(oer).(ptl).H2Oa(1);
            state.(inm).H2Oa = aw*ones(nc_inm, 1);
            cH2O = IonomerMembrane.groupHydration(model.(inm), aw, T);
            model.(inm).H2O.c0 = cH2O/aw;

            state.(inm).H2Oceps = cH2O.*model.(inm).volumeFraction;

            state.(inm).phi = zeros(nc_inm, 1); % OBS : this is not the value that is finally assigned (only sent here so that dispatchIonomerToReactionLayers can be used).
            state = model.evalVarName(state, {her, ctl, 'Einmr'});
            state = model.evalVarName(state, {oer, ctl, 'Einmr'});

            Eelyte_oer = state.(oer).(ctl).Eelyte(1);
            Einmr_oer  = state.(oer).(ctl).Einmr(1);
            Eelyte_her = state.(her).(ctl).Eelyte(1);
            Einmr_her  = state.(her).(ctl).Einmr(1);

            nc_her = model.(her).(ptl).G.getNumberOfCells();
            nc_oer = model.(oer).(ptl).G.getNumberOfCells();

            state.(her).(ptl).E   = 0;
            state.(her).(ptl).phi = - Eelyte_her*ones(nc_her, 1);
            state.(inm).phi       = - Einmr_her*ones(nc_inm, 1);
            state.(oer).(ptl).E   = Einmr_oer - Einmr_her;
            state.(oer).(ptl).phi = (Einmr_oer - Einmr_her - Eelyte_oer)*ones(nc_oer, 1);

            bd = 'Boundary';

            state = model.evalVarName(state, '.*compGasMasses');
            state = model.evalVarName(state, 'Por.*phasePressures');
            state = model.evalVarName(state, 'Layer.(?!Boundary).*liqrho$');

            for ielde = 1 : numel(eldes)
                elde = eldes{ielde};
                coupterm = model.(elde).(ptl).externalCouplingTerm;
                coupcells = coupterm.couplingcells;
                nc = numel(coupcells);

                gasInd = model.(elde).(ptl).gasInd;
                gvf = state.(elde).(ptl).volumeFractions{model.(elde).(ptl).phaseInd.gas};
                for igas = 1 :  gasInd.ngas
                    compgasdensity = state.(elde).(ptl).compGasMasses{igas}./gvf;
                    compgasdensity = compgasdensity(coupcells);
                    state.(elde).(ptl).(bd).gasDensities{igas} = compgasdensity;
                    controlValues.gasDensities{igas} = compgasdensity;
                end

                mobPhaseInd = model.(elde).(ptl).mobPhaseInd;
                for imobphase = 1 : mobPhaseInd.nmobphase
                    iphase = mobPhaseInd.phaseMap(imobphase);
                    p = state.(elde).(ptl).phasePressures{iphase};
                    p = p(coupcells);
                    state.(elde).(ptl).(bd).phasePressures{iphase} = p;
                    controlValues.mobilePhasePressures{imobphase} = p;
                end

                liqrho = state.(elde).(ptl).liqrho;
                liqrho = liqrho(coupcells);
                state.(elde).(ptl).(bd).liqrho = liqrho;
                controlValues.liqrho = liqrho;

                liquidInd = model.(elde).(ptl).liquidInd;
                cOH = state.(elde).(ptl).concentrations{liquidInd.OH};
                cOH = cOH(coupcells);
                state.(elde).(ptl).(bd).cOH = cOH;
                controlValues.cOH = cOH;

                model.(elde).(ptl).(bd).controlValues = controlValues;

            end

            if model.(oer).(ctl).include_dissolution
                
                dm = 'DissolutionModel';
                nc = model.(oer).(ctl).(dm).G.getNumberOfCells();
                
                state.(oer).(ctl).(dm).volumeFraction = model.(oer).(ctl).(dm).volumeFraction0*ones(nc, 1);
                
            end
            

        end

        function model = setupForSimulation(model)

            shortNames = {'IonomerMembrane'           , 'inm';
                          'HydrogenEvolutionElectrode', 'her';
                          'OxygenEvolutionElectrode'  , 'oer';
                          'CatalystLayer'             , 'ctl';
                          'ExchangeReaction'          , 'exr';
                          'PorousTransportLayer'      , 'ptl';
                          'Boundary'                  , 'bd' ;
                          'DissolutionModel'          , 'dm'};

            model = model.equipModelForComputation('shortNames', shortNames);

            F = model.con.F;

            
            inm = 'IonomerMembrane';
            her = 'HydrogenEvolutionElectrode';
            oer = 'OxygenEvolutionElectrode';

            ctl = 'CatalystLayer';
            exr = 'ExchangeReaction';
            ptl = 'PorousTransportLayer';

            scalings = {{{inm, 'chargeCons'}, F}, ...
                        {{her, ptl, 'chargeCons'}, F}};

            if model.(oer).(ctl).include_dissolution
                scalings = horzcat(scalings, ...
                                   {{{oer, ctl, dm, 'massCons'}, 1e5}});
            end

        end
        
        
        function state = setupControl(model, state, drivingForces)

            oer = 'OxygenEvolutionElectrode';
            her = 'HydrogenEvolutionElectrode';
            ctl = 'CatalystLayer';

            time = state.time;

            controlI = drivingForces.src(time);

            controlEqs{1} = state.(oer).(ctl).I - controlI;
            controlEqs{2} = state.(her).(ctl).E; % we impose zero potential at cathode

            state.controlEqs = controlEqs;

        end

        function state = dispatchTemperature(model, state)

            oer = 'OxygenEvolutionElectrode';
            her = 'HydrogenEvolutionElectrode';
            inm = 'IonomerMembrane';

            state.(oer).T = state.T(model.(oer).G.mappings.cellmap);
            state.(her).T = state.T(model.(her).G.mappings.cellmap);
            state.(inm).T = state.T(model.(inm).G.mappings.cellmap);

        end

        function state = dispatchIonomerToReactionLayers(model, state)

            oer = 'OxygenEvolutionElectrode';
            her = 'HydrogenEvolutionElectrode';
            inm = 'IonomerMembrane';

            ctl = 'CatalystLayer';
            exr = 'ExchangeReaction';

            eldes    = {her, oer};
            layers   = {ctl, exr};
            varnames = {'H2OaInmr', 'cOHinmr', 'phiInmr'};

            for ielde = 1 : numel(eldes)

                elde = eldes{ielde};

                % setup initial value for the variable that also takes care of AD.
                nc = model.(elde).(ctl).G.getNumberOfCells();
                initval = nan(nc, 1);
                [adsample, isAD] = getSampleAD(state.(inm).phi);
                if isAD
                    initval = adsample.convertDouble(initval);
                end

                coupterms = model.couplingTerms;
                coupnames = model.couplingNames;
                coupname  = sprintf('%s-%s', elde, inm);
                coupterm  = getCoupTerm(coupterms, coupname, coupnames);
                coupcells = coupterm.couplingcells;

                for ilayer = 1 : numel(layers)
                    layer = layers{ilayer};
                    for ivar = 1 : numel(varnames)
                        varname = varnames{ivar};
                        state.(elde).(layer).(varname) = initval;
                    end
                    state.(elde).(layer).H2OaInmr(coupcells(:, 1)) = state.(inm).H2Oa(coupcells(:, 2));
                    state.(elde).(layer).cOHinmr(coupcells(:, 1))  = state.(inm).cOH(coupcells(:, 2));
                    state.(elde).(layer).phiInmr(coupcells(:, 1))  = state.(inm).phi(coupcells(:, 2));
                end

            end

        end

        function state = updateIonomerSources(model, state)

            inm = 'IonomerMembrane';
            her = 'HydrogenEvolutionElectrode';
            oer = 'OxygenEvolutionElectrode';

            ctl = 'CatalystLayer';
            exr = 'ExchangeReaction';

            % initialize sources (done in a way that takes care of AD, meaning that OHsource and H2OSource inherits AD
            % structure from state.(inm).H2Oceps)
            OHsource  = 0*state.(inm).H2Oceps;
            H2OSource = 0*state.(inm).H2Oceps;
            
            vols = model.(inm).G.getVolumes();

            eldes = {her, oer};
            
            for ielde = 1 : numel(eldes)

                elde = eldes{ielde};

                coupterms = model.couplingTerms;
                coupnames = model.couplingNames;

                coupname  = sprintf('%s-%s', elde, inm);
                coupterm  = getCoupTerm(coupterms, coupname, coupnames);
                coupcells = coupterm.couplingcells;
                vols      = vols(coupcells(:, 2));

                inmrOHsource  = state.(elde).(ctl).inmrOHsource(coupcells(:, 1));
                inmrH2Osource = state.(elde).(ctl).inmrH2Osource(coupcells(:, 1));
                OHexchR       = state.(elde).(exr).OHexchangeRate(coupcells(:, 1));
                H2OexchR      = state.(elde).(exr).H2OexchangeRate(coupcells(:, 1));

                OHsource(coupcells(:, 2))  = vols.*(inmrOHsource - OHexchR);
                H2OSource(coupcells(:, 2)) = vols.*(inmrH2Osource - H2OexchR);

            end

            state.(inm).OHsource  = OHsource;
            state.(inm).H2OSource = H2OSource;

        end


        function cleanState = addStaticVariables(model, cleanState, state)

            oer = 'OxygenEvolutionElectrode';
            her = 'HydrogenEvolutionElectrode';
            ptl = 'PorousTransportLayer';
            inm = 'IonomerMembrane';

            cleanState = addStaticVariables@BaseModel(model, cleanState, state);
            cleanState.T = state.T;
            cleanState.(inm).jBcSource = 0;
            cleanState.(oer).(ptl).jBcSource = 0;
            cleanState.(her).(ptl).jBcSource = 0;

        end

        function forces = getValidDrivingForces(model)

            forces = getValidDrivingForces@PhysicalModel(model);
            forces.src = [];

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

            [~, state] = getEquations(model, state0, state, dt, drivingForces, 'ResOnly', true);


        end


        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            
            [state, report] = updateState@BaseModel(model, state, problem, dx, drivingForces);

            docapping = true;

            if docapping
                
                inm = 'IonomerMembrane';
                her = 'HydrogenEvolutionElectrode';
                oer = 'OxygenEvolutionElectrode';
                ctl = 'CatalystLayer';
                exr = 'ExchangeReaction';
                ptl = 'PorousTransportLayer';
                dm  = 'DissolutionModel';
                bd  = 'Boundary';

                varnames = {{her, ptl, 'H2rhoeps'}                   , ...
                            {her, ptl, 'liqrhoeps'}                  , ...
                            {her, ptl, 'liqeps'}                     , ...
                            {her, ptl, 'H2Ogasrhoeps'}               , ...
                            {her, ptl, 'OHceps'}                     , ...
                            {her, ptl, 'Boundary', 'cOH'            }, ...
                            {her, ptl, 'Boundary', 'liqrho'         }, ...
                            {her, ptl, 'Boundary', 'gasDensities', 2}, ...
                            {her, ptl, 'Boundary', 'gasDensities', 1}, ...
                            {oer, ptl, 'O2rhoeps'}                   , ...
                            {oer, ptl, 'liqrhoeps'}                  , ...
                            {oer, ptl, 'liqeps'}                     , ...
                            {oer, ptl, 'H2Ogasrhoeps'}               , ...
                            {oer, ptl, 'OHceps'}                     , ...
                            {oer, ptl, 'Boundary', 'cOH'}            , ...
                            {oer, ptl, 'Boundary', 'liqrho'}         , ...
                            {oer, ptl, 'Boundary', 'gasDensities', 2}, ...
                            {oer, ptl, 'Boundary', 'gasDensities', 1}, ...
                            {inm, 'H2Oceps'}};

                if model.(oer).(ctl).include_dissolution
                    varname = {oer, ctl, dm, 'volumeFraction'};
                    val = model.getProp(state, varname);
                    val = max(0, val);
                    state = model.setProp(state, varname, val);
                end
                
                for ivar = 1 : numel(varnames)

                    varname = varnames{ivar};
                    val = model.getProp(state, varname);
                    val = max(0, val);
                    state = model.setProp(state, varname, val);
                    
                end

            end
        end
        
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)

            [state, report] = updateAfterConvergence@BaseModel(model, state0, state, dt, drivingForces);

            state = model.evalVarName(state, 'compGasMasses');

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
