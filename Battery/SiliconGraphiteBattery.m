classdef SiliconGraphiteBattery < Battery

    methods
        
        function model = SiliconGraphiteBattery(paramobj)
            
            model = model@Battery(paramobj);

            assert(model.use_thermal == false, 'thermal model not covered yet');
            assert(model.include_current_collectors == false, 'current collectors  not included yet');
            
        end


        function model = setupSelectedModel(model);
            % defines shorthands for the submodels
            elyte   = 'Electrolyte';
            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            am      = 'ActiveMaterial';
            cc      = 'CurrentCollector';
            itf     = 'Interface';
            sd      = 'SolidDiffusion';
            thermal = 'ThermalModel';
            ctrl    = 'Control';

            gr = 'Graphite';
            si = 'Silicon';

            varEqTypes ={{elyte, 'c'}                 , 'elyte_massCons'         , 'cell'; ...  
                         {elyte, 'phi'}               , 'elyte_chargeCons'       , 'cell'; ...    
                         {ne, am, 'phi'}              , 'ne_am_chargeCons'       , 'cell'; ...    
                         {pe, am, 'phi'}              , 'pe_am_chargeCons'       , 'cell'; ...    
                         {ctrl, 'E'}                  , 'EIeq'                   , 'ctrl'; ...  
                         {ctrl, 'I'}                  , 'controlEq'              , 'ctrl'; ...
                         {ne, am, gr, sd, 'c'}        , 'ne_am_gr_sd_massCons'   , 'cell'; ...
                         {ne, am, gr, sd, 'cSurface'} , 'ne_am_gr_sd_soliddiffeq', 'cell'; ...
                         {ne, am, si, sd, 'c'}        , 'ne_am_si_sd_massCons'   , 'cell'; ...
                         {ne, am, si, sd, 'cSurface'} , 'ne_am_si_sd_soliddiffeq', 'cell'; ...
                         {pe, am, sd, 'c'}            , 'pe_am_sd_massCons'      , 'cell'; ...
                         {pe, am, sd, 'cSurface'}     , 'pe_am_sd_soliddiffeq'   , 'cell' ...
                         };

            addedVariableNames = {};
            addedVariableNames{end + 1} = {ctrl, 'ctrlType'};
            
            if strcmp(model.(ctrl).controlPolicy, 'CCCV')
                addedVariableNames{end + 1} = {ctrl, 'nextCtrlType'};
            end

            primaryVariableNames = varEqTypes(:, 1);
            equationNames        = varEqTypes(:, 2);
            equationTypes        = varEqTypes(:, 3);

            equationIndices = struct();
            for ieq = 1 : numel(equationNames)
                equationIndices.(equationNames{ieq}) = ieq;
            end
            
            model.addedVariableNames   = addedVariableNames;
            model.primaryVariableNames = primaryVariableNames; 
            model.equationNames        = equationNames; 
            model.equationTypes        = equationTypes;
            model.equationIndices      = equationIndices;
            
        end
        
        function model = registerVarAndPropfuncNames(model)
            
        %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)
            
            model = registerVarAndPropfuncNames@Battery(model);
            
            % defines shorthands for the submodels
            elyte   = 'Electrolyte';
            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            am      = 'ActiveMaterial';
            cc      = 'CurrentCollector';
            am      = 'ActiveMaterial';
            itf     = 'Interface';
            sd      = 'SolidDiffusion';
            thermal = 'ThermalModel';
            ctrl    = 'Control';

            gr = 'Graphite';
            si = 'Silicon';

            %%  removing some unused variables
            varnames = {{ne, am, itf, 'phiElectrolyte'},
                        {ne, am, itf, 'cElectrolyte'}};

            model = model.removeVarNames(varnames);
            
            %% Coupling functions
            
            % Dispatch electrolyte concentration and potential in the electrodes
            fn = @Battery.updateElectrodeCoupling;
            inputnames = {{elyte, 'c'}, ...
                          {elyte, 'phi'}};
            model = model.registerPropFunction({{ne, am, gr, itf, 'phiElectrolyte'}, fn, inputnames});
            model = model.registerPropFunction({{ne, am, si, itf, 'phiElectrolyte'}, fn, inputnames});
            model = model.registerPropFunction({{ne, am, gr, itf, 'cElectrolyte'}  , fn, inputnames});
            model = model.registerPropFunction({{ne, am, si, itf, 'cElectrolyte'}  , fn, inputnames});
            model = model.registerPropFunction({{pe, am, itf, 'phiElectrolyte'}, fn, inputnames});
            model = model.registerPropFunction({{pe, am, itf, 'cElectrolyte'}  , fn, inputnames});
            
            % Functions that update the source terms in the electolyte
            fn = @Battery.updateElectrolyteCoupling;
            inputnames = {{ne, am, si, 'Rvol'}, ...
                          {ne, am, gr, 'Rvol'}, ...
                          {pe, am, 'Rvol'}};
            model = model.registerPropFunction({{elyte, 'massSource'}, fn, inputnames});
            model = model.registerPropFunction({{elyte, 'eSource'}, fn, inputnames});
            
        end
        
        function initstate = setupInitialState(model)
        % Setup the values of the primary variables at initial state

            nc = model.G.cells.num;

            T   = model.initT;
            
            bat = model;
            elyte   = 'Electrolyte';
            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            am      = 'ActiveMaterial';
            itf     = 'Interface';
            sd      = 'SolidDiffusion';
            cc      = 'CurrentCollector';
            thermal = 'ThermalModel';
            ctrl    = 'Control';

            gr = 'Graphite';
            si = 'Silicon';

            initstate.(thermal).T = T*ones(nc, 1);

            %% Synchronize temperatures
            initstate = model.updateTemperature(initstate);

            %% setup initial state for negative electrode

            mats = {gr, si};

            for imat = 1 : numel(mats)
                mat = mats{imat};
                % set primary variables
                N = model.(ne).(am).(mat).(sd).N;
                np = model.(pe).(am).(sd).np; % Note : we have by construction np = nc
                cElectrodeInit                          = (model.(ne).(am).(mat).(itf).theta0)*(model.(ne).(am).(mat).(itf).cmax);
                initstate.(ne).(am).(mat).(sd).c        = cElectrodeInit*ones(N*np, 1);
                initstate.(ne).(am).(mat).(sd).cSurface = cElectrodeInit*ones(np, 1);
            end

            initstate.(ne).(am).(gr) = model.(ne).(am).(gr).updateConcentrations(initstate.(ne).(am).(gr));
            initstate.(ne).(am).(gr).(itf) = model.(ne).(am).(gr).(itf).updateOCP(initstate.(ne).(am).(gr).(itf));
            
            ref = initstate.(ne).(am).(gr).(itf).OCP(1);

            nc = bat.(ne).(am).G.cells.num;
            
            initstate.(ne).(am).phi = zeros(nc, 1);
            

            %% setup initial state for positive electrode
            
            pe_itf = bat.(pe).(am).(itf); 

            cElectrodeInit = pe_itf.cmax;
            nc             = pe_itf.G.cells.num;

            initstate.(pe).(am).(sd).cSurface = cElectrodeInit*ones(nc, 1);
            N = model.(pe).(am).(sd).N;
            np = model.(pe).(am).(sd).np; % Note : we have by construction np = nc
            initstate.(pe).(am).(sd).c = c*ones(N*np, 1);

            initstate.(pe).(am) = model.(pe).(am).updateConcentrations(initstate.(pe).(am));
            initstate.(pe).(am).(itf) = pe_itf.updateOCP(initstate.(pe).(am).(itf));

            initstate.(pe).(am).phi = initstate.(pe).(am).(itf).OCP - ref;
            
            %% Setup initial Electrolyte state

            initstate.(elyte).phi = zeros(bat.(elyte).G.cells.num, 1) - ref;
            initstate.(elyte).c   = 1000*ones(bat.(elyte).G.cells.num, 1);

            initstate.(ctrl).E = initstate.(pe).(am).phi(1) - initstate.(ne).(am).phi(1);
            initstate.(ctrl).I = - model.(ctrl).Imax;
            
            switch model.(ctrl).controlPolicy
              case 'CCCV'
                initstate.(ctrl).ctrlType = 'CC_charge1';
                initstate.(ctrl).nextCtrlType = 'CC_charge1';
              case 'IEswitch'
                initstate.(ctrl).ctrlType = 'constantCurrent';
              otherwise
                error('control policy not recognized');
            end
            
        end
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
        % Assembly of the governing equation
            
            opts = struct('ResOnly', false, 'iteration', 0, 'reverseMode', false); 
            opts = merge_options(opts, varargin{:});
            
            time = state0.time + dt;
            if(not(opts.ResOnly) && not(opts.reverseMode))
                state = model.initStateAD(state);
            elseif(opts.reverseMode)
                disp('No AD initatlization in equation old style')
                state0 = model.initStateAD(state0);
            else
                assert(opts.ResOnly);
            end

            %% start replace
            
            state.(elyte)              = model.(elyte).updateCurrentBcSource(state.(elyte));
            state.(ne).(am)            = model.(ne).(am).updateConductivity(state.(ne).(am));
            state.(ne).(am)            = model.(ne).(am).updatejCoupling(state.(ne).(am));
            state.(pe).(am)            = model.(pe).(am).updateConductivity(state.(pe).(am));
            state.(pe).(am)            = model.(pe).(am).updatejCoupling(state.(pe).(am));
            state.(ctrl)               = model.(ctrl).updateControlEquation(state.(ctrl));
            state                      = model.setupEIEquation(state);
            state                      = model.setupExternalCouplingPositiveElectrode(state);
            state.(pe).(am)            = model.(pe).(am).updatejBcSource(state.(pe).(am));
            state.(pe).(am)            = model.(pe).(am).updateCurrent(state.(pe).(am));
            state                      = model.updateTemperature(state);
            state.(pe).(am).(sd)       = model.(pe).(am).(sd).updateMassAccum(state.(pe).(am).(sd), state0.(pe).(am).(sd), dt);
            state.(pe).(am).(sd)       = model.(pe).(am).(sd).updateAverageConcentration(state.(pe).(am).(sd));
            state.(pe).(am)            = model.(pe).(am).updateSOC(state.(pe).(am));
            state.(pe).(am)            = model.(pe).(am).dispatchTemperature(state.(pe).(am));
            state.(pe).(am).(sd)       = model.(pe).(am).(sd).updateDiffusionCoefficient(state.(pe).(am).(sd));
            state.(pe).(am).(sd)       = model.(pe).(am).(sd).updateFlux(state.(pe).(am).(sd));
            state.(pe).(am)            = model.(pe).(am).updateConcentrations(state.(pe).(am));
            state.(pe).(am)            = model.(pe).(am).updatePhi(state.(pe).(am));
            state.(pe).(am).(itf)      = model.(pe).(am).(itf).updateOCP(state.(pe).(am).(itf));
            state                      = model.setupExternalCouplingNegativeElectrode(state);
            state.(ne).(am)            = model.(ne).(am).updatejBcSource(state.(ne).(am));
            state.(ne).(am)            = model.(ne).(am).updateCurrent(state.(ne).(am));
            state.(ne).(am)            = model.(ne).(am).dispatchTemperature(state.(ne).(am));
            state.(ne).(am).(si).(sd)  = model.(ne).(am).(si).(sd).updateMassAccum(state.(ne).(am).(si).(sd), state0.(ne).(am).(si).(sd), dt);
            state.(ne).(am).(si).(sd)  = model.(ne).(am).(si).(sd).updateAverageConcentration(state.(ne).(am).(si).(sd));
            state.(ne).(am).(si)       = model.(ne).(am).(si).updateSOC(state.(ne).(am).(si));
            state.(ne).(am).(si)       = model.(ne).(am).(si).dispatchTemperature(state.(ne).(am).(si));
            state.(ne).(am).(si).(sd)  = model.(ne).(am).(si).(sd).updateDiffusionCoefficient(state.(ne).(am).(si).(sd));
            state.(ne).(am).(si).(sd)  = model.(ne).(am).(si).(sd).updateFlux(state.(ne).(am).(si).(sd));
            state.(ne).(am).(si)       = model.(ne).(am).(si).updateConcentrations(state.(ne).(am).(si));
            state.(ne).(am)            = model.(ne).(am).updatePhi(state.(ne).(am));
            state.(ne).(am).(si).(itf) = model.(ne).(am).(si).(itf).updateOCP(state.(ne).(am).(si).(itf));
            state.(ne).(am).(gr).(sd)  = model.(ne).(am).(gr).(sd).updateMassAccum(state.(ne).(am).(gr).(sd), state0.(ne).(am).(gr).(sd), dt);
            state.(ne).(am).(gr).(sd)  = model.(ne).(am).(gr).(sd).updateAverageConcentration(state.(ne).(am).(gr).(sd));
            state.(ne).(am).(gr)       = model.(ne).(am).(gr).updateSOC(state.(ne).(am).(gr));
            state.(ne).(am).(gr)       = model.(ne).(am).(gr).dispatchTemperature(state.(ne).(am).(gr));
            state.(ne).(am).(gr).(sd)  = model.(ne).(am).(gr).(sd).updateDiffusionCoefficient(state.(ne).(am).(gr).(sd));
            state.(ne).(am).(gr).(sd)  = model.(ne).(am).(gr).(sd).updateFlux(state.(ne).(am).(gr).(sd));
            state.(ne).(am).(gr)       = model.(ne).(am).(gr).updateConcentrations(state.(ne).(am).(gr));
            state.(ne).(am).(gr).(itf) = model.(ne).(am).(gr).(itf).updateOCP(state.(ne).(am).(gr).(itf));
            state.(elyte)              = model.(elyte).updateConcentrations(state.(elyte));
            state.(elyte)              = model.(elyte).assembleAccumTerm(state.(elyte), state0.(elyte), dt);
            state                      = model.updateElectrodeCoupling(state);
            state.(pe).(am).(itf)      = model.(pe).(am).(itf).updateReactionRateCoefficient(state.(pe).(am).(itf));
            state.(pe).(am).(itf)      = model.(pe).(am).(itf).updateEta(state.(pe).(am).(itf));
            state.(pe).(am).(itf)      = model.(pe).(am).(itf).updateReactionRate(state.(pe).(am).(itf));
            state.(pe).(am)            = model.(pe).(am).updateRvol(state.(pe).(am));
            state.(pe).(am)            = model.(pe).(am).updateCurrentSource(state.(pe).(am));
            state.(pe).(am)            = model.(pe).(am).updateChargeConservation(state.(pe).(am));
            state.(pe).(am).(sd)       = model.(pe).(am).(sd).updateMassSource(state.(pe).(am).(sd));
            state.(pe).(am).(sd)       = model.(pe).(am).(sd).assembleSolidDiffusionEquation(state.(pe).(am).(sd));
            state.(pe).(am).(sd)       = model.(pe).(am).(sd).updateMassConservation(state.(pe).(am).(sd));
            state.(ne).(am).(si).(itf) = model.(ne).(am).(si).(itf).updateReactionRateCoefficient(state.(ne).(am).(si).(itf));
            state.(ne).(am).(si).(itf) = model.(ne).(am).(si).(itf).updateEta(state.(ne).(am).(si).(itf));
            state.(ne).(am).(si).(itf) = model.(ne).(am).(si).(itf).updateReactionRate(state.(ne).(am).(si).(itf));
            state.(ne).(am).(si)       = model.(ne).(am).(si).updateRvol(state.(ne).(am).(si));
            state.(ne).(am).(si).(sd)  = model.(ne).(am).(si).(sd).updateMassSource(state.(ne).(am).(si).(sd));
            state.(ne).(am).(si).(sd)  = model.(ne).(am).(si).(sd).assembleSolidDiffusionEquation(state.(ne).(am).(si).(sd));
            state.(ne).(am).(si).(sd)  = model.(ne).(am).(si).(sd).updateMassConservation(state.(ne).(am).(si).(sd));
            state.(ne).(am).(gr).(itf) = model.(ne).(am).(gr).(itf).updateReactionRateCoefficient(state.(ne).(am).(gr).(itf));
            state.(ne).(am).(gr).(itf) = model.(ne).(am).(gr).(itf).updateEta(state.(ne).(am).(gr).(itf));
            state.(ne).(am).(gr).(itf) = model.(ne).(am).(gr).(itf).updateReactionRate(state.(ne).(am).(gr).(itf));
            state.(ne).(am).(gr)       = model.(ne).(am).(gr).updateRvol(state.(ne).(am).(gr));
            state.(ne).(am)            = model.(ne).(am).updateCurrentSource(state.(ne).(am));
            state.(ne).(am)            = model.(ne).(am).updateChargeConservation(state.(ne).(am));
            state                      = model.updateElectrolyteCoupling(state);
            state.(ne).(am).(gr).(sd)  = model.(ne).(am).(gr).(sd).updateMassSource(state.(ne).(am).(gr).(sd));
            state.(ne).(am).(gr).(sd)  = model.(ne).(am).(gr).(sd).assembleSolidDiffusionEquation(state.(ne).(am).(gr).(sd));
            state.(ne).(am).(gr).(sd)  = model.(ne).(am).(gr).(sd).updateMassConservation(state.(ne).(am).(gr).(sd));
            state.(elyte)              = model.(elyte).updateChemicalCurrent(state.(elyte));
            state.(elyte)              = model.(elyte).updateDiffusionCoefficient(state.(elyte));
            state.(elyte)              = model.(elyte).updateConductivity(state.(elyte));
            state.(elyte)              = model.(elyte).updateCurrent(state.(elyte));
            state.(elyte)              = model.(elyte).updateMassFlux(state.(elyte));
            state.(elyte)              = model.(elyte).updateMassConservation(state.(elyte));
            state.(elyte)              = model.(elyte).updateChargeConservation(state.(elyte));

            %% end replace

            eqs = cell(1, numel(model.equationNames));
            ei = model.equationIndices;

            massConsScaling = model.con.F;
            
            eqs{ei.elyte_massCons}          = state.(elyte).massCons*massConsScaling;
            eqs{ei.elyte_chargeCons}        = state.(elyte).chargeCons;
            eqs{ei.ne_am_chargeCons}        = state.(ne).(am).chargeCons;
            eqs{ei.pe_am_chargeCons}        = state.(pe).(am).chargeCons;
            eqs{ei.EIeq}                    = state.(ctrl).EIeq;
            eqs{ei.controlEq}               = state.(ctrl).controlEquation;

            n    = model.(ne).(am).(gr).(itf).n; % number of electron transfer (equal to 1 for Lithium)
            F    = model.con.F;
            vol  = model.(ne).(am).(gr).operators.pv;
            rp   = model.(ne).(am).(gr).(sd).rp;
            vsf  = model.(ne).(am).(gr).(sd).volumetricSurfaceArea;
            surfp = 4*pi*rp^2;
            
            scalingcoef = (vsf*vol(1)*n*F)/surfp;
            
            eqs{ei.ne_am_gr_sd_massCons}    = state.(ne).(am).(gr).(sd).massCons;
            eqs{ei.ne_am_gr_sd_soliddiffeq} = state.(ne).(am).(gr).(sd).solidDiffusionEq;

            n    = model.(ne).(am).(si).(itf).n; % number of electron transfer (equal to 1 for Lithium)
            F    = model.con.F;
            vol  = model.(ne).(am).(si).operators.pv;
            rp   = model.(ne).(am).(si).(sd).rp;
            vsf  = model.(ne).(am).(si).(sd).volumetricSurfaceArea;
            surfp = 4*pi*rp^2;
            
            scalingcoef = (vsf*vol(1)*n*F)/surfp;
                        
            eqs{ei.ne_am_si_sd_massCons}    = state.(ne).(am).(si).(sd).massCons;        
            eqs{ei.ne_am_si_sd_soliddiffeq} = state.(ne).(am).(si).(sd).solidDiffusionEq;

            n    = model.(ne).(am).(itf).n; % number of electron transfer (equal to 1 for Lithium)
            F    = model.con.F;
            vol  = model.(ne).(am).operators.pv;
            rp   = model.(ne).(am).(sd).rp;
            vsf  = model.(ne).(am).(sd).volumetricSurfaceArea;
            surfp = 4*pi*rp^2;
            
            scalingcoef = (vsf*vol(1)*n*F)/surfp;

            eqs{ei.pe_am_sd_massCons}       = state.(ne).(am).(sd).massCons;        
            eqs{ei.pe_am_sd_soliddiffeq}    = state.(ne).(am).(sd).solidDiffusionEq;

            names = model.equationNames;
            types = model.equationTypes;
            
            primaryVars = model.getPrimaryVariables();
            
            %% Setup LinearizedProblem that can be processed by MRST Newton API
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
            
        end
        
        function state = updateTemperature(model, state)
        % Dispatch the temperature in all the submodels

            elyte   = 'Electrolyte';
            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            am      = 'ActiveMaterial';
            cc      = 'CurrentCollector';
            thermal = 'ThermalModel';

            gr = 'Graphite';
            si = 'Silicon';
            
            % (here we assume that the ThermalModel has the "parent" grid)
            state.(elyte).T   = state.(thermal).T(model.(elyte).G.mappings.cellmap);
            state.(ne).(am).T = state.(thermal).T(model.(ne).(am).G.mappings.cellmap);
            state.(pe).(am).T = state.(thermal).T(model.(pe).(am).G.mappings.cellmap);
            
            % Update temperature in the active materials of the electrodes.
            state.(ne).(am) = model.(ne).(am).dispatchTemperature(state.(ne).(am));
            state.(ne).(am).(si) = model.(ne).(am).(si).dispatchTemperature(state.(ne).(am).(si));
            state.(ne).(am).(gr) = model.(ne).(am).(gr).dispatchTemperature(state.(ne).(am).(gr));
            state.(pe).(am) = model.(pe).(am).dispatchTemperature(state.(pe).(am));
            
        end
        
        function state = updateElectrolyteCoupling(model, state)
        % Assemble the electrolyte coupling by adding the ion sources from the electrodes
            
            battery = model;

            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            am    = 'ActiveMaterial';
            itf   = 'Interface';
            
            vols = battery.(elyte).G.cells.volumes;
            F = battery.con.F;
            ne_vsa = battery.(ne).(am).(itf).volumetricSurfaceArea;
            pe_vsa = battery.(pe).(am).(itf).volumetricSurfaceArea;
            
            
            couplingterms = battery.couplingTerms;

            elyte_c_source = zeros(battery.(elyte).G.cells.num, 1);
            elyte_e_source = zeros(battery.(elyte).G.cells.num, 1);
            
            % setup AD 
            phi = state.(elyte).phi;
            if isa(phi, 'ADI')
                adsample = getSampleAD(phi);
                adbackend = model.AutoDiffBackend;
                elyte_c_source = adbackend.convertToAD(elyte_c_source, adsample);
            end
            
            coupnames = model.couplingNames;
            
            ne_si_Rvol = state.(ne).(am).(si).Rvol;
            ne_gr_Rvol = state.(ne).(am).(gr).Rvol;
            if isa(ne_Rvol, 'ADI') & ~isa(elyte_c_source, 'ADI')
                adsample = getSampleAD(ne_Rvol);
                adbackend = model.AutoDiffBackend;
                elyte_c_source = adbackend.convertToAD(elyte_c_source, adsample);
            end
            
            coupterm = getCoupTerm(couplingterms, 'NegativeElectrode-Electrolyte', coupnames);
            elytecells = coupterm.couplingcells(:, 2);
            elyte_c_source(elytecells) = (ne_si_Rvol + ne_gr_Rvol).*vols(elytecells);
            
            pe_Rvol = state.(pe).(am).Rvol;
            if isa(pe_Rvol, 'ADI') & ~isa(elyte_c_source, 'ADI')
                adsample = getSampleAD(pe_Rvol);
                adbackend = model.AutoDiffBackend;
                elyte_c_source = adbackend.convertToAD(elyte_c_source, adsample);
            end
            
            coupterm = getCoupTerm(couplingterms, 'PositiveElectrode-Electrolyte', coupnames);
            elytecells = coupterm.couplingcells(:, 2);
            elyte_c_source(elytecells) = pe_Rvol.*vols(elytecells);
            
            elyte_e_source = elyte_c_source.*battery.(elyte).sp.z(1)*F; 
            
            state.Electrolyte.massSource = elyte_c_source; 
            state.Electrolyte.eSource = elyte_e_source;
            
        end
        

        function state = updateControl(model, state, drivingForces)
            
            ctrl = "Control";
            
            switch model.(ctrl).controlPolicy
              case 'CCCV'
                % nothing to do here
              case 'IEswitch'
                
                E    = state.(ctrl).E;
                I    = state.(ctrl).I;
                time = state.time;
                
                [ctrlVal, ctrltype] = drivingForces.src(time, value(I), value(E));
                
                state.(ctrl).ctrlVal  = ctrlVal;
                state.(ctrl).ctrlType = ctrltype;

              case 'None'
                % nothing done here. This case is only used for addVariables function
              otherwise
                error('control type not recognized');
            end
                
            
        end
        
        function state = updateElectrodeCoupling(model, state)
        % Setup the electrode coupling by updating the potential and concentration of the electrolyte in the active
        % component of the electrodes. There, those quantities are considered as input and used to compute the reaction
        % rate.

            bat = model;
            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            am    = 'ActiveMaterial';
            itf   = 'Interface';

            gr = 'Graphite';
            si = 'Silicon';
            
            phi_elyte = state.(elyte).phi;
            c_elyte = state.(elyte).cs{1};
            
            elyte_cells = zeros(model.G.cells.num, 1);
            elyte_cells(bat.(elyte).G.mappings.cellmap) = (1 : bat.(elyte).G.cells.num)';

            state.(ne).(am).(si).(itf).phiElectrolyte = phi_elyte(elyte_cells(bat.(ne).(am).(si).G.mappings.cellmap));
            state.(ne).(am).(si).(itf).cElectrolyte   = c_elyte(elyte_cells(bat.(ne).(am).(si).G.mappings.cellmap));
            state.(ne).(am).(gr).(itf).phiElectrolyte = phi_elyte(elyte_cells(bat.(ne).(am).(gr).G.mappings.cellmap));
            state.(ne).(am).(gr).(itf).cElectrolyte   = c_elyte(elyte_cells(bat.(ne).(am).(gr).G.mappings.cellmap));s
            
            state.(pe).(am).(itf).phiElectrolyte = phi_elyte(elyte_cells(bat.(pe).(am).G.mappings.cellmap));
            state.(pe).(am).(itf).cElectrolyte = c_elyte(elyte_cells(bat.(pe).(am).G.mappings.cellmap));
            
        end

        function [state, report] = updateState(model, state, problem, dx, drivingForces)

            [state, report] = updateState@BaseModel(model, state, problem, dx, drivingForces);
            
            %% cap concentrations
            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            am    = 'ActiveMaterial';
            sd    = 'SolidDiffusion';
            itf   = 'Interface';

            cmin = model.cmin;
            
            state.(elyte).c = max(cmin, state.(elyte).c);
            
            eldes = {ne, pe};
            for ind = 1 : numel(eldes)
                elde = eldes{ind};
                if strcmp(model.(elde).(am).diffusionModelType, 'simple') | ~model.use_particle_diffusion
                    state.(elde).(am).c = max(cmin, state.(elde).(am).c);
                    cmax = model.(elde).(am).(itf).cmax;
                    state.(elde).(am).c = min(cmax, state.(elde).(am).c);
                else
                    state.(elde).(am).(sd).c = max(cmin, state.(elde).(am).(sd).c);
                    cmax = model.(elde).(am).(itf).cmax;
                    state.(elde).(am).(sd).c = min(cmax, state.(elde).(am).(sd).c);
                end
            end
            
            ctrl = 'Control';            
            state.(ctrl) = model.(ctrl).updateControlState(state.(ctrl));
            
            report = [];
            
        end

        function cleanState = addStaticVariables(model, cleanState, state)
        % Variables that are no AD initiated (but should be "carried over")
            
            cleanState = addStaticVariables@BaseModel(model, cleanState, state);
            
            thermal = 'ThermalModel';
            ctrl = 'Control';
            
            cleanState.(ctrl).ctrlType = state.(ctrl).ctrlType;            
            
            if ~model.use_thermal
                thermal = 'ThermalModel';
                cleanState.(thermal).T = state.(thermal).T;
            end
            
            
        end

        function [model, state] = prepareTimestep(model, state, state0, dt, drivingForces)
            
            [model, state] = prepareTimestep@BaseModel(model, state, state0, dt, drivingForces);
            
            ctrl = 'Control';
            state.(ctrl) = model.(ctrl).prepareStepControl(state.(ctrl), state0.(ctrl), dt, drivingForces);
            
        end
        
        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)
            
            [state, report] = updateAfterConvergence@BaseModel(model, state0, state, dt, drivingForces);
            
            ctrl = 'Control';
            state.(ctrl) = model.(ctrl).updateControlAfterConvergence(state.(ctrl), state0.(ctrl), dt);
        end

        function model = setupCapping(model)
            model.cmin = 0;
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
