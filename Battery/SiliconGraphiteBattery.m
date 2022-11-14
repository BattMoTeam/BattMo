classdef SiliconGraphiteBattery < Battery

    methods
        
        function model = SiliconGraphiteBattery(paramobj)
            
            model = model@Battery(paramobj);

            assert(model.use_thermal == false, 'thermal model not covered yet');
            assert(model.include_current_collectors == false, 'current collectors  not included yet');
            
        end


        function model = setupSelectedModel(model);

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
                         {pe, am, gr, sd, 'c'}        , 'pe_am_gr_sd_massCons'   , 'cell'; ...
                         {pe, am, gr, sd, 'cSurface'} , 'pe_am_gr_sd_soliddiffeq', 'cell'; ...
                         {pe, am, si, sd, 'c'}        , 'pe_am_si_sd_massCons'   , 'cell'; ...
                         {pe, am, si, sd, 'cSurface'} , 'pe_am_si_sd_soliddiffeq', 'cell'; ...
                         };

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
            inputnames = {{ne, am, gr, 'Rvol'}, ...
                          {ne, am, si, 'Rvol'}, ...
                          {pe, am, 'Rvol'}};
            model = model.registerPropFunction({{elyte, 'massSource'}, fn, inputnames});
            model = model.registerPropFunction({{elyte, 'eSource'}, fn, inputnames});
            
        end

        
        function control = setupControl(model, paramobj)

            switch paramobj.controlPolicy
              case "IEswitch"
                control = IEswitchControlModel(paramobj); 
              case "CCCV"
                control = CcCvControlModel(paramobj);
              otherwise
                error('Error controlPolicy not recognized');
            end
            
            C = computeCellCapacity(model);
            CRate = control.CRate;
            
            control.Imax = (C/hour)*CRate;
            
        end
        
        function model = setupMappings(model)
            
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            am    = 'ActiveMaterial';
            cc    = 'CurrentCollector';
            elyte = 'Electrolyte';
            
            G_elyte = model.(elyte).G;
            elytecelltbl.cells = (1 : G_elyte.cells.num)';
            elytecelltbl.globalcells = G_elyte.mappings.cellmap;
            elytecelltbl = IndexArray(elytecelltbl);

            eldes = {ne, pe};

            for ind = 1 : numel(eldes)

                elde = eldes{ind};
                G_elde  = model.(elde).(am).G;
                clear eldecelltbl;
                eldecelltbl.cells = (1 : G_elde.cells.num)';
                eldecelltbl.globalcells = G_elde.mappings.cellmap;
                eldecelltbl = IndexArray(eldecelltbl);
                
                map = TensorMap();
                map.fromTbl = elytecelltbl;
                map.toTbl = eldecelltbl;
                map.replaceFromTblfds = {{'cells', 'elytecells'}};
                map.replaceToTblfds = {{'cells', 'eldecells'}};
                map.mergefds = {'globalcells'};
                
                mappings.(elde) = map.getDispatchInd();
                
            end
            
            model.mappings = mappings;
            
        end
        

        
        function initstate = setupInitialState(model)
        % Setup the values of the primary variables at initial state

            nc = model.G.cells.num;

            SOC = model.SOC;
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
            
            initstate.(thermal).T = T*ones(nc, 1);

            %% Synchronize temperatures
            initstate = model.updateTemperature(initstate);

            
            %% Setup initial state for NegativeElectrode
            
            eldes = {ne, pe};
            
            for ind = 1 : numel(eldes)
                
                elde = eldes{ind};
                
                elde_itf = bat.(elde).(am).(itf); 

                theta = SOC*(elde_itf.theta100 - elde_itf.theta0) + elde_itf.theta0;
                c     = theta*elde_itf.cmax;
                nc    = elde_itf.G.cells.num;

                switch model.(elde).(am).diffusionModelType
                  case 'simple'
                    initstate.(elde).(am).(sd).cSurface = c*ones(nc, 1);
                    initstate.(elde).(am).c = c*ones(nc, 1);
                  case 'full'
                    initstate.(elde).(am).(sd).cSurface = c*ones(nc, 1);
                    N = model.(elde).(am).(sd).N;
                    np = model.(elde).(am).(sd).np; % Note : we have by construction np = nc
                    initstate.(elde).(am).(sd).c = c*ones(N*np, 1);
                end
                
                initstate.(elde).(am) = model.(elde).(am).updateConcentrations(initstate.(elde).(am));
                initstate.(elde).(am).(itf) = elde_itf.updateOCP(initstate.(elde).(am).(itf));

                OCP = initstate.(elde).(am).(itf).OCP;
                if ind == 1
                    % The value in the first cell is used as reference.
                    ref = OCP(1);
                end
                
                initstate.(elde).(am).phi = OCP - ref;
                
            end

            %% Setup initial Electrolyte state

            initstate.(elyte).phi = zeros(bat.(elyte).G.cells.num, 1)-ref;
            initstate.(elyte).c = 1000*ones(bat.(elyte).G.cells.num, 1);

            %% Setup initial Current collectors state

            if model.(ne).include_current_collector
                OCP = initstate.(ne).(am).(itf).OCP;
                OCP = OCP(1) .* ones(bat.(ne).(cc).G.cells.num, 1);
                initstate.(ne).(cc).phi = OCP - ref;
            end
            
            if model.(pe).include_current_collector
                OCP = initstate.(pe).(am).(itf).OCP;
                OCP = OCP(1) .* ones(bat.(pe).(cc).G.cells.num, 1);
                initstate.(pe).(cc).phi = OCP - ref;
            end
            
            initstate.(ctrl).E = OCP(1) - ref;
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

        function state = addVariables(model, state)

            error('not updated');
            
            % Given a state where only the primary variables are defined, this
            % functions add all the additional variables that are computed in the assembly process and have some physical
            % interpretation.
            %
            % To do so, we use getEquations function and sends dummy variable for state0, dt and drivingForces 
            
            % Values that need to be set to get the function getEquations running
            
            dt = 1;
            state0 = state;
            model.Control = ControlModel([]);
            model.Control.controlPolicy = 'None';
            drivingForces = model.getValidDrivingForces();

            % We call getEquations to update state
            
            [~, state] = getEquations(model, state0, state, dt, drivingForces, 'ResOnly', true);

            % we set to empty the fields we know that are not meaningfull

            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            am      = 'ActiveMaterial';
            cc      = 'CurrentCollector';
            elyte   = 'Electrolyte';
            am      = 'ActiveMaterial';
            itf     = 'Interface';
            sd      = "SolidDiffusion";
            thermal = 'ThermalModel';
            ctrl    = 'Control';

            state.(elyte).massAccum = [];
            state.(elyte).massCons = [];
            eldes = {ne, pe};
            for ielde = 1 : numel(eldes)
                elde = eldes{ielde};
                switch model.(elde).(am).diffusionModelType
                  case 'full'
                    state.(elde).(am).(sd).massAccum = [];
                    state.(elde).(am).(sd).massCons = [];
                  case 'simple'
                    % nothing to remove here
                  otherwise
                    error('diffusion model type not recognized');
                end
            end
            
            for ielde = 1 : numel(eldes)
                elde = eldes{ielde};
                state.(elde).(am).(sd) = model.(elde).(am).(sd).updateAverageConcentration(state.(elde).(am).(sd));
                state.(elde).(am) = model.(elde).(am).updateSOC(state.(elde).(am));
                state.(elde).(am) = model.(elde).(am).updateAverageConcentration(state.(elde).(am));
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
            
           
            names = model.equationNames;
            types = model.equationTypes;
            
            %% The equations are reordered in a way that is consitent with the linear iterative solver 
            % (the order of the equation does not matter if we do not use an iterative solver)
            ctrltype = state.Control.ctrlType;
            switch ctrltype
              case {'constantCurrent', 'CC_discharge1', 'CC_discharge2', 'CC_charge1'}
                types{ei.EIeq} = 'cell';   
              case {'constantVoltage', 'CV_charge2'}
                % no changes
              otherwise 
                error('control type not recognized')
            end

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
            
            % (here we assume that the ThermalModel has the "parent" grid)
            state.(elyte).T   = state.(thermal).T(model.(elyte).G.mappings.cellmap);
            state.(ne).(am).T = state.(thermal).T(model.(ne).(am).G.mappings.cellmap);
            state.(pe).(am).T = state.(thermal).T(model.(pe).(am).G.mappings.cellmap);
            if model.include_current_collectors
                state.(ne).(cc).T = state.(thermal).T(model.(ne).(cc).G.mappings.cellmap);
                state.(pe).(cc).T = state.(thermal).T(model.(pe).(cc).G.mappings.cellmap);
            end
            
            % Update temperature in the active materials of the electrodes.
            state.(ne).(am) = model.(ne).(am).dispatchTemperature(state.(ne).(am));
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
            
            ne_Rvol = state.(ne).(am).Rvol;
            if isa(ne_Rvol, 'ADI') & ~isa(elyte_c_source, 'ADI')
                adsample = getSampleAD(ne_Rvol);
                adbackend = model.AutoDiffBackend;
                elyte_c_source = adbackend.convertToAD(elyte_c_source, adsample);
            end
            
            coupterm = getCoupTerm(couplingterms, 'NegativeElectrode-Electrolyte', coupnames);
            elytecells = coupterm.couplingcells(:, 2);
            elyte_c_source(elytecells) = ne_Rvol.*vols(elytecells);
            
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
        
        function state = updateThermalOhmicSourceTerms(model, state)
        % Assemble the ohmic source term :code:`state.jHeatOhmSource`, see :cite:t:`Latz2016`

            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            am      = 'ActiveMaterial';
            cc      = 'CurrentCollector';
            elyte   = 'Electrolyte';
            thermal = 'ThermalModel';
            
            eldes = {ne, pe}; % electrodes
            
            nc = model.G.cells.num;
            
            src = zeros(nc, 1);
                        
            for ind = 1 : numel(eldes)

                elde = eldes{ind};

                if model.include_current_collectors
                    cc_model = model.(elde).(cc);
                    cc_map   = cc_model.G.mappings.cellmap;
                    cc_j     = state.(elde).(cc).jFace;
                    cc_econd = cc_model.EffectiveElectricalConductivity;
                    cc_vols  = cc_model.G.cells.volumes;
                    cc_jsq   = computeCellFluxNorm(cc_model, cc_j); 
                    state.(elde).(cc).jsq = cc_jsq;  %store square of current density
                    src = subsetPlus(src,cc_vols.*cc_jsq./cc_econd, cc_map);
                    %                    src(cc_map) = src(cc_map) + cc_vols.*cc_jsq./cc_econd;
                end
                
                am_model = model.(elde).(am);
                am_map   = am_model.G.mappings.cellmap;
                am_j     = state.(elde).(am).jFace;
                am_econd = am_model.EffectiveElectricalConductivity;
                am_vols  = am_model.G.cells.volumes;
                am_jsq   = computeCellFluxNorm(am_model, am_j);
                state.(elde).(am).jsq = am_jsq;
                
                %src(am_map) = src(am_map) + am_vols.*am_jsq./am_econd;
                src = subsetPlus(src, am_vols.*am_jsq./am_econd, am_map);
            end

            % Electrolyte
            elyte_model    = model.(elyte);
            elyte_map      = elyte_model.G.mappings.cellmap;
            elyte_vf       = elyte_model.volumeFraction;
            elyte_j        = state.(elyte).jFace;
            elyte_bruggman = elyte_model.BruggemanCoefficient;
            elyte_cond     = state.(elyte).conductivity;
            elyte_econd    = elyte_cond.*elyte_vf.^elyte_bruggman;
            elyte_vols     = elyte_model.G.cells.volumes;
            elyte_jsq      = computeCellFluxNorm(elyte_model, elyte_j);
            state.(elyte).jsq = elyte_jsq; %store square of current density
            
            %src(elyte_map) = src(elyte_map) + elyte_vols.*elyte_jsq./elyte_econd;
            src = subsetPlus(src, elyte_vols.*elyte_jsq./elyte_econd, elyte_map);
            state.(thermal).jHeatOhmSource = src;
            
        end
        
        function state = updateThermalChemicalSourceTerms(model, state)
        % Assemble the thermal source term from transport :code:`state.jHeatChemicalSource`, see :cite:t:`Latz2016`
            
            elyte = 'Electrolyte';
            thermal = 'ThermalModel';
            
            % prepare term
            nc = model.G.cells.num;
            src = zeros(nc, 1);
            T = state.(thermal).T;
            phi = state.(elyte).phi;
            nf = model.(elyte).G.faces.num;
            intfaces = model.(elyte).operators.internalConn;
            if isa(T, 'ADI')
                adsample = getSampleAD(T);
                adbackend = model.AutoDiffBackend;
                src = adbackend.convertToAD(src, adsample);
                zeroFace = model.AutoDiffBackend.convertToAD(zeros(nf, 1), phi);
                locstate = state;
            else
                locstate = value(state);
                zeroFace = zeros(nf, 1);
            end

            % Compute chemical heat source in electrolyte
            dmudcs = locstate.(elyte).dmudcs;   % Derivative of chemical potential with respect to concentration
            D      = locstate.(elyte).D;        % Effective diffusion coefficient 
            Dgradc = locstate.(elyte).diffFlux; % Diffusion flux (-D*grad(c))
            DFaceGradc = zeroFace;
            DFaceGradc(intfaces) = Dgradc;
            
            
            % compute norm of square norm of diffusion flux
            elyte_model   = model.(elyte);
            elyte_map     = elyte_model.G.mappings.cellmap;
            elyte_vols    = elyte_model.G.cells.volumes;
            elyte_jchemsq = computeCellFluxNorm(elyte_model, DFaceGradc);
            elyte_src     = elyte_vols.*elyte_jchemsq./D;
            
            % This is a bit hacky for the moment (we should any way consider all the species)
            elyte_src = dmudcs{1}.*elyte_src;
            
            % map to source term at battery level
            src(elyte_map) = src(elyte_map) + elyte_src;
            
            state.(thermal).jHeatChemicalSource = src;
            
        end
        
        
        function state = updateThermalReactionSourceTerms(model, state)
        % Assemble the source term from chemical reaction :code:`state.jHeatReactionSource`, see :cite:t:`Latz2016`            
            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            am      = 'ActiveMaterial';
            itf     = 'Interface';
            thermal = 'ThermalModel';
            
            
            eldes = {ne, pe}; % electrodes
            
            nc = model.G.cells.num;
            
            src = zeros(nc, 1);
           
            T = state.(thermal).T;
            if isa(T, 'ADI')
                adsample = getSampleAD(T);
                adbackend = model.AutoDiffBackend;
                src = adbackend.convertToAD(src, adsample);
                locstate = state;
            else
               locstate = value(state); 
            end
            
            for ind = 1 : numel(eldes)

                elde = eldes{ind};
                
                itf_model = model.(elde).(am).(itf);
                
                F       = itf_model.constants.F;
                n       = itf_model.n;
                itf_map = itf_model.G.mappings.cellmap;
                vsa     = itf_model.volumetricSurfaceArea;
                vols    = model.(elde).(am).G.cells.volumes;

                Rvol = locstate.(elde).(am).Rvol;
                dUdT = locstate.(elde).(am).(itf).dUdT;
                eta  = locstate.(elde).(am).(itf).eta;
                
                itf_src = n*F*vols.*Rvol.*(eta + T(itf_map).*dUdT);
                
                src(itf_map) = src(itf_map) + itf_src;
                
            end

            state.(thermal).jHeatReactionSource = src;

        end
        
        
        function state = updateElectrodeCoupling(model, state)
        % Setup the electrode coupling by updating the potential and concentration of the electrolyte in the active
        % component of the electrodes. There, those quantities are considered as input and used to compute the reaction
        % rate.
        %
        % WARNING : at the moment, we do not pass the concentrations

            
            bat = model;
            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            am    = 'ActiveMaterial';
            itf   = 'Interface';
            cc    = 'CurrentCollector';
            
            eldes = {ne, pe};
            phi_elyte = state.(elyte).phi;
            c_elyte = state.(elyte).cs{1};
            
            elyte_cells = zeros(model.G.cells.num, 1);
            elyte_cells(bat.(elyte).G.mappings.cellmap) = (1 : bat.(elyte).G.cells.num)';

            for ind = 1 : numel(eldes)
                elde = eldes{ind};
                state.(elde).(am).(itf).phiElectrolyte = phi_elyte(elyte_cells(bat.(elde).(am).G.mappings.cellmap));
                state.(elde).(am).(itf).cElectrolyte = c_elyte(elyte_cells(bat.(elde).(am).G.mappings.cellmap));
            end
            
        end

        function state = setupExternalCouplingNegativeElectrode(model, state)
        %
        % Setup external electronic coupling of the negative electrode at the current collector
        %
            ne = 'NegativeElectrode';
            am = 'ActiveMaterial';
            
            if model.(ne).include_current_collector
                
                cc = 'CurrentCollector';

                phi   = state.(ne).(cc).phi;
                sigma = state.(ne).(cc).conductivity;

                [jExternal, jFaceExternal] = setupExternalCoupling(model.(ne).(cc), phi, 0, sigma);
                
                state.(ne).(cc).jExternal = jExternal;
                state.(ne).(cc).jFaceExternal = jFaceExternal;
                state.(ne).(am).jExternal     = 0;
                state.(ne).(am).jFaceExternal = 0;
                
            else
                
                phi   = state.(ne).(am).phi;
                sigma = state.(ne).(am).conductivity;
                
                [jExternal, jFaceExternal] = setupExternalCoupling(model.(ne).(am), phi, 0, sigma);
                
                state.(ne).(am).jExternal = jExternal;
                state.(ne).(am).jFaceExternal = jFaceExternal;
                
            end
            
            
        end
        
        function state = setupExternalCouplingPositiveElectrode(model, state)
        %
        % Setup external electronic coupling of the positive electrode at the current collector
        %            
            pe   = 'PositiveElectrode';
            ctrl = 'Control';
            am   = 'ActiveMaterial';
            
            E   = state.(ctrl).E;

            if model.(pe).include_current_collector
                
                cc   = 'CurrentCollector';
                
                phi   = state.(pe).(cc).phi;
                sigma = state.(pe).(cc).conductivity;
                
                [jExternal, jFaceExternal] = setupExternalCoupling(model.(pe).(cc), phi, E, sigma);
                
                state.(pe).(cc).jExternal = jExternal;
                state.(pe).(cc).jFaceExternal = jFaceExternal;
                state.(pe).(am).jExternal     = 0;
                state.(pe).(am).jFaceExternal = 0;
            else
                
                phi   = state.(pe).(am).phi;
                sigma = state.(pe).(am).conductivity;
                
                [jExternal, jFaceExternal] = setupExternalCoupling(model.(pe).(am), phi, E, sigma);
                
                state.(pe).(am).jExternal = jExternal;
                state.(pe).(am).jFaceExternal = jFaceExternal;
                
            end
            
        end
        
        function state = setupEIEquation(model, state)
            
            pe   = 'PositiveElectrode';
            ctrl = 'Control';
            
            I = state.(ctrl).I;
            E = state.(ctrl).E;
            
            if model.include_current_collectors
                cc   = 'CurrentCollector';
                
                phi = state.(pe).(cc).phi;
                
                coupterm = model.(pe).(cc).externalCouplingTerm;
                faces    = coupterm.couplingfaces;
                cond_pcc = model.(pe).(cc).EffectiveElectricalConductivity;
                [trans_pcc, cells] = model.(pe).(cc).operators.harmFaceBC(cond_pcc, faces);
                
                state.Control.EIequation = sum(trans_pcc.*(state.(pe).(cc).phi(cells) - E)) - I;

            else
                
                am = 'ActiveMaterial';
                
                phi = state.(pe).(am).phi;
                
                coupterm = model.(pe).(am).externalCouplingTerm;
                faces    = coupterm.couplingfaces;
                cond_pcc = model.(pe).(am).EffectiveElectricalConductivity;
                [trans_pcc, cells] = model.(pe).(am).operators.harmFaceBC(cond_pcc, faces);
                
                state.Control.EIequation = sum(trans_pcc.*(state.(pe).(am).phi(cells) - E)) - I;

            end
            
        end
        
        function state = initStateAD(model, state)
        % initialize a new cleaned-up state with AD variables

            % initStateAD in BaseModel erase all fields
            newstate = initStateAD@BaseModel(model, state);

            % add the variable that we want to add on state
            % TODO : pass those using addStaticVariables method
            addedvarnames = model.addedVariableNames;
            for i = 1 : numel(addedvarnames)
                var = model.getProp(state, addedvarnames{i});
                assert(isnumeric(var) | ischar(var));
                newstate = model.setNewProp(newstate, addedvarnames{i}, var);
            end

            newstate.time = state.time;
            
            state = newstate;
            
        end 
        

        function primaryvarnames = getPrimaryVariables(model)

            primaryvarnames = model.primaryVariableNames;
            
        end
        
        function forces = getValidDrivingForces(model)
            
            forces = getValidDrivingForces@PhysicalModel(model);
            
            ctrl = 'Control';
            switch model.(ctrl).controlPolicy
              case 'CCCV'
                forces.CCCV = true;
              case 'IEswitch'
                forces.IEswitch = true;
                forces.src = [];
              case 'None'
                % used only in addVariables
              otherwise
                error('Error controlPolicy not recognized');
            end
            % TODO this is a hack to get thing go
            forces.Imax = [];
            
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
        
        
        function outputvars = extractGlobalVariables(model, states)
            
            ns = numel(states);

            for i = 1 : ns
                E    = states{i}.Control.E;
                I    = states{i}.Control.I;
                time = states{i}.time;
                
                outputvars{i} = struct('E'   , E   , ...
                                       'I'   , I   , ...
                                       'time', time);
                if model.use_thermal
                    T    = states{i}.ThermalModel.T; 
                    outputvars{i}.Tmax = max(T);
                end
            
            end
        end

    end
    
    methods(Static)
        
        function [found, varind] = getVarIndex(varname, pvarnames)

            varname = strjoin(varname, '_');
            pvarnames = cellfun(@(name) strjoin(name, '_'), pvarnames, 'uniformoutput', false);

            [found, varind] = ismember(varname, pvarnames);

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
