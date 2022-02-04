classdef IsothermalBattery < BaseModel
% 
% The battery model consists of 
%
% * an Electrolyte model given in :attr:`Electrolyte` property
% * a Negative Electrode Model given in :attr:`NegativeElectrode` property
% * a Positive Electrode Model given in :attr:`PositiveElectrode` property
% * a Thermal model given in :attr:`ThermalModel` property
%
    properties
        
        con = PhysicalConstants();

        Electrolyte       % Electrolyte model, instance of :class:`Electrolyte <Electrochemistry.Electrodes.Electrolyte>`
        NegativeElectrode % Negative Electrode Model, instance of :class:`Electrode <Electrochemistry.Electrodes.Electrode>`
        PositiveElectrode % Positive Electrode Model, instance of :class:`Electrode <Electrochemistry.Electrodes.Electrode>`
        
        SOC % State Of Charge

        initT % Initial temperature
        
        Ucut % Voltage cut
        
        couplingTerms % Coupling terms
        cmin % mininum concentration used in capping

        couplingNames 
        
        mappings
        
        use_solid_diffusion
    end
    
    methods
        
        function model = IsothermalBattery(paramobj,varargin)
            opt = struct('use_solid_diffusion',true);
            opt = merge_options(opt,varargin{:});
            model = model@BaseModel();
            
            % All the submodels should have same backend (this is not assigned automaticallly for the moment)
            model.AutoDiffBackend = SparseAutoDiffBackend('useBlocks', false);
            
            %% Setup the model using the input parameters
            fdnames = {'G'             , ...
                       'couplingTerms' , ...
                       'initT'         , ...
                       'SOC'           , ...
                       'I'             , ...
                       'Ucut'};
            
            model = dispatchParams(model, paramobj, fdnames);
            
            % Assign the components : Electrolyte, NegativeElectrode, PositiveElectrode
            model.NegativeElectrode = model.setupElectrode(paramobj.NegativeElectrode);
            model.PositiveElectrode = model.setupElectrode(paramobj.PositiveElectrode);
            model.Electrolyte = model.setupElectrolyte(paramobj.Electrolyte);
            
            % setup Electrolyte model (setup electrolyte volume fractions in the different regions)
            model = model.setupElectrolyteModel();            
            
            % setup couplingNames
            model.couplingNames = cellfun(@(x) x.name, model.couplingTerms, 'uniformoutput', false);
            
            % setup some mappings (mappings from electrodes to electrolyte)
            model = model.setupMappings();
            
            % setup capping
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            eac   = 'ElectrodeActiveComponent';
            am    = 'ActiveMaterial';
            cmax_ne = model.(ne).(eac).(am).Li.cmax;
            cmax_pe = model.(pe).(eac).(am).Li.cmax;
            model.cmin = 1e-5*max(cmax_ne, cmax_pe);
            model.use_solid_diffusion = opt.use_solid_diffusion; 
        end

        
        function electrode = setupElectrode(model, paramobj)
        % Setup the electrode models (both :attr:`NegativeElectrode` and :attr:`PositiveElectrode`). Here, :code:`paramobj`
        % is instance of :class:`ElectrodeInputParams <Electrochemistry.Electrodes.ElectrodeInputParams>`
            electrode = Electrode(paramobj);
        end
        
        function electrolyte = setupElectrolyte(model, paramobj)
        % Setup the electrolyte model :attr:`Electrolyte`. Here, :code:`paramobj` is instance of
        % :class:`ElectrolyteInputParams <Electrochemistry.ElectrolyteInputParams>`
            switch paramobj.electrolyteType
              case 'binary'
                electrolyte = Electrolyte(paramobj);
              otherwise
                % binary is default
                electrolyte = Electrolyte(paramobj)
            end
        end
        
        function model = setupMappings(model)
            
            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            eac     = 'ElectrodeActiveComponent';
            cc      = 'CurrentCollector';
            elyte   = 'Electrolyte';
            
            
            G_elyte = model.(elyte).G;
            elytecelltbl.cells = (1 : G_elyte.cells.num)';
            elytecelltbl.globalcells = G_elyte.mappings.cellmap;
            elytecelltbl = IndexArray(elytecelltbl);

            eldes = {ne, pe};

            for ind = 1 : numel(eldes)

                elde = eldes{ind};
                G_elde  = model.(elde).(eac).G;
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
        
        function model = setupElectrolyteModel(model)
        % Assign the electrolyte volume fractions in the different regions

            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            eac   = 'ElectrodeActiveComponent';
            sep   = 'Separator';

            elyte_cells = zeros(model.G.cells.num, 1);
            elyte_cells(model.(elyte).G.mappings.cellmap) = (1 : model.(elyte).G.cells.num)';

            model.(elyte).volumeFraction = NaN(model.(elyte).G.cells.num, 1);
            model.(elyte).volumeFraction(elyte_cells(model.(ne).(eac).G.mappings.cellmap))  = model.(ne).(eac).porosity;
            model.(elyte).volumeFraction(elyte_cells(model.(pe).(eac).G.mappings.cellmap))  = model.(pe).(eac).porosity;
            model.(elyte).volumeFraction(elyte_cells(model.(elyte).(sep).G.mappings.cellmap)) = model.(elyte).(sep).porosity;

        end
        
        function initstate = setupInitialState(model)
        % Setup the initial state

            nc = model.G.cells.num;

            SOC = model.SOC;
            T   = model.initT;
            
            initstate.SOC = SOC*ones(nc, 1);
            initstate.T = T*ones(nc, 1);
            
            bat = model;
            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            am    = 'ActiveMaterial';
            eac   = 'ElectrodeActiveComponent';
            cc    = 'CurrentCollector';
                        
            %% synchronize temperatures
            initstate = model.updateTemperature(initstate);
            
            %% setup initial NegativeElectrode state
            
            % shortcut
            % negAm : ActiveMaterial of the negative electrode
            
            negAm = bat.(ne).(eac).(am); 
            
            m = (1 ./ (negAm.theta100 - negAm.theta0));
            b = -m .* negAm.theta0;
            theta = (SOC - b) ./ m;
            c = theta .* negAm.Li.cmax;
            c = c*ones(negAm.G.cells.num, 1);

            initstate.(ne).(eac).c = c;
            % We bypass the solid diffusion equation to set directly the particle surface concentration (this is a bit hacky)
            initstate.(ne).(eac).(am).cElectrode = c;
            initstate.(ne).(eac).(am) = negAm.updateOCP(initstate.(ne).(eac).(am));

            OCP = initstate.(ne).(eac).(am).OCP;
            ref = OCP(1);
            
            initstate.(ne).(eac).phi = OCP - ref;

            %% setup initial PositiveElectrode state

            % shortcut
            % posAm : ActiveMaterial of the positive electrode
            
            posAm = bat.(pe).(eac).(am);
            
            m = (1 ./ (posAm.theta100 - posAm.theta0));
            b = -m .* posAm.theta0;
            theta = (SOC - b) ./ m;
            c = theta .* posAm.Li.cmax;
            c = c*ones(posAm.G.cells.num, 1);

            initstate.(pe).(eac).c = c;
            % We bypass the solid diffusion equation to set directly the particle surface concentration (this is a bit hacky)
            initstate.(pe).(eac).(am).cElectrode = c;
            initstate.(pe).(eac).(am) = posAm.updateOCP(initstate.(pe).(eac).(am));
            
            OCP = initstate.(pe).(eac).(am).OCP;
            initstate.(pe).(eac).phi = OCP - ref;

            %% setup initial Electrolyte state

            initstate.(elyte).phi = zeros(bat.(elyte).G.cells.num, 1)-ref;
            cs = cell(2,1);
            initstate.(elyte).cs = cs;
            initstate.(elyte).cs{1} = 1000*ones(bat.(elyte).G.cells.num, 1);

            %% setup initial Current collectors state

            OCP = initstate.(ne).(eac).(am).OCP;
            OCP = OCP(1) .* ones(bat.(ne).(cc).G.cells.num, 1);
            initstate.(ne).(cc).phi = OCP - ref;

            OCP = initstate.(pe).(eac).(am).OCP;
            OCP = OCP(1) .* ones(bat.(pe).(cc).G.cells.num, 1);
            initstate.(pe).(cc).phi = OCP - ref;
            
            initstate.(pe).(cc).E = OCP(1) - ref;
            initstate.(pe).(cc).I = 0;
            
        end
        
        function [problem, state] = getEquations(model, state0, state,dt, drivingForces, varargin)
        % Assembly of the governing equation
            opts = struct('ResOnly', false, 'iteration', 0); 
            opts = merge_options(opts, varargin{:});
            
            time = state0.time + dt;
            if(not(opts.ResOnly))
                state = model.initStateAD(state);
            end
            
            %% for now temperature and SOC are kept constant
            nc = model.G.cells.num;
            state.SOC = model.SOC*ones(nc, 1);
            
            % Shortcuts used in this function
            battery = model;
            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            eac     = 'ElectrodeActiveComponent';
            cc      = 'CurrentCollector';
            elyte   = 'Electrolyte';
            am      = 'ActiveMaterial';
            
            electrodes = {ne, pe};
            electrodecomponents = {eac, cc};

            %% Synchronization across components

            % temperature
            state = battery.updateTemperature(state);

            state.(elyte) = battery.(elyte).updateConcentrations(state.(elyte));
            
            for ind = 1 : numel(electrodes)
                elde = electrodes{ind};
                % potential and concentration between active material and electode active component
                state.(elde).(eac) = battery.(elde).(eac).updatePhi(state.(elde).(eac));
                if(model.use_solid_diffusion)
                    state.(elde).(eac) = battery.(elde).(eac).updateChargeCarrier(state.(elde).(eac));
                else
                    state.(elde).(eac).(am).cElectrode = state.(elde).(eac).c;
                end
                
            end
            
            %% Accumulation terms

            state = battery.updateAccumTerms(state, state0, dt);

            %% Update Electrolyte -> Electrodes coupling 
            
            state = battery.updateElectrodeCoupling(state); 

            %% Update reaction rates in both electrodes

            for ind = 1 : numel(electrodes)
                elde = electrodes{ind};
                state.(elde).(eac).(am) = battery.(elde).(eac).(am).updateReactionRateCoefficient(state.(elde).(eac).(am));
                %if(model.use_solid_diffusion)
                    state.(elde).(eac).(am) = battery.(elde).(eac).(am).updateDiffusionCoefficient(state.(elde).(eac).(am));
                %end
                state.(elde).(eac).(am) = battery.(elde).(eac).(am).updateOCP(state.(elde).(eac).(am));
                state.(elde).(eac).(am) = battery.(elde).(eac).(am).updateReactionRate(state.(elde).(eac).(am));
            end

            %% Update Electrodes -> Electrolyte  coupling

            state = battery.updateElectrolyteCoupling(state);
            
            %% Update coupling within electrodes and external coupling
            
            state.(ne) = battery.(ne).updateCoupling(state.(ne));
            state.(pe) = battery.(pe).updateCoupling(state.(pe));

            state.(ne).(eac) = battery.(ne).(eac).updatejBcSource(state.(ne).(eac));
            state.(pe).(eac) = battery.(pe).(eac).updatejBcSource(state.(pe).(eac));
            
            state = model.setupExternalCouplingNegativeElectrode(state);
            state = model.setupExternalCouplingPositiveElectrode(state);
            
            state.(ne).(cc) = battery.(ne).(cc).updatejBcSource(state.(ne).(cc));
            state.(pe).(cc) = battery.(pe).(cc).updatejBcSource(state.(pe).(cc));
            
            %% elyte charge conservation

            state.(elyte) = battery.(elyte).updateCurrentBcSource(state.(elyte));
            state.(elyte) = battery.(elyte).updateConductivity(state.(elyte));
            state.(elyte) = battery.(elyte).updateChemicalCurrent(state.(elyte));
            state.(elyte) = battery.(elyte).updateCurrent(state.(elyte));
            state.(elyte) = battery.(elyte).updateChargeConservation(state.(elyte));

            %% Electrodes charge conservation - Active material part

            for ind = 1 : numel(electrodes)
                elde = electrodes{ind};
                state.(elde).(eac) = battery.(elde).(eac).updateIonAndCurrentSource(state.(elde).(eac));
                state.(elde).(eac) = battery.(elde).(eac).updateCurrent(state.(elde).(eac));
                state.(elde).(eac) = battery.(elde).(eac).updateChargeConservation(state.(elde).(eac));
            end
            
            %% elyte mass conservation

            state.(elyte) = battery.(elyte).updateDiffusionCoefficient(state.(elyte));
            state.(elyte) = battery.(elyte).updateChargeCarrierFlux(state.(elyte));
            state.(elyte) = battery.(elyte).updateMassConservation(state.(elyte));

            
            for ind = 1 : numel(electrodes)
                elde = electrodes{ind};
                
                %% Electrodes mass conservation
                state.(elde).(eac) = battery.(elde).(eac).updateChargeCarrierFlux(state.(elde).(eac));
                state.(elde).(eac) = battery.(elde).(eac).updateMassConservation(state.(elde).(eac));
                
                %% Electrodes charge conservation - current collector part
                state.(elde).(cc) = battery.(elde).(cc).updateCurrent(state.(elde).(cc));
                state.(elde).(cc) = battery.(elde).(cc).updateChargeConservation(state.(elde).(cc));

            end

            %% update solid diffustion equations
            if(model.use_solid_diffusion)
                for ind = 1 : numel(electrodes)
                    elde = electrodes{ind};
                    state.(elde).(eac).(am) = battery.(elde).(eac).(am).assembleSolidDiffusionEquation(state.(elde).(eac).(am));
                end
            end
            
            %% setup relation between E and I at positive current collectror
            
            state = model.setupEIEquation(state);
            
            %% Set up the governing equations
            
            eqs = {};
            
            %% We collect mass and charge conservation equations for the electrolyte and the electrodes


            massConsScaling = model.con.F;    
                        
            eqs{end + 1} = state.(elyte).massCons*massConsScaling;
            eqs{end + 1} = state.(elyte).chargeCons;
            
            resistance = 1/model.(ne).(eac).EffectiveElectricalConductivity(1);
            eqs{end + 1} = state.(ne).(eac).massCons*massConsScaling;
            eqs{end + 1} = state.(ne).(eac).chargeCons;
            if(model.use_solid_diffusion)
                eqs{end + 1} = state.(ne).(eac).(am).solidDiffusionEq;
            end
            
            resistance = 1/model.(pe).(eac).EffectiveElectricalConductivity(1);
            eqs{end + 1} = state.(pe).(eac).massCons*massConsScaling;
            eqs{end + 1} = state.(pe).(eac).chargeCons;
            if(model.use_solid_diffusion)
                eqs{end + 1} = massConsScaling*state.(pe).(eac).(am).solidDiffusionEq;
            end
            
            resistance = 1/model.(ne).(cc).EffectiveElectricalConductivity(1);
            eqs{end + 1} = state.(ne).(cc).chargeCons;
            eqs{end + 1} = state.(pe).(cc).chargeCons;
            
            eqs{end + 1} = state.EIeq;
            
            % we add the control equation
            I = state.(pe).(cc).I;
            E = state.(pe).(cc).E;
            [val, ctrltype] = drivingForces.src(time, value(I), value(E));
            switch ctrltype
              case 'I'
                eqs{end + 1} = I - val;
              case 'E'
                eqs{end + 1} = (E - val)*max(model.NegativeElectrode.CurrentCollector.operators.T_all);
            end

            %% Give type and names to equations and names of the primary variables (for book-keeping)
            
            
            if(model.use_solid_diffusion)                
                names = {'elyte_massCons'   , ...
                     'elyte_chargeCons' , ...
                     'ne_eac_massCons'  , ...
                     'ne_eac_chargeCons', ...
                     'ne_eac_am_soliddiffeq', ...
                     'pe_eac_massCons'  , ...
                     'pe_eac_chargeCons', ...
                     'pe_eac_am_soliddiffeq', ...
                     'ne_cc_chargeCons' , ...
                     'pe_cc_chargeCons' , ...
                     'EIeq', ...
                     'controlEq'};
                  types = {'cell','cell','cell','cell', 'sdiff','cell','cell','sdiff','cell','cell', 'cntrl', 'cntrl'};
                  switch ctrltype
                     case 'I'
                        types{end-1} = 'cell';   
                    case 'E'    
                        order = [1:9,11,10];
                        eqs = {eqs{order}};
                        names = {names{order}};
                      otherwise 
                          error()
                  end
            else
                names = {'elyte_massCons'   , ...
                     'elyte_chargeCons' , ...
                     'ne_eac_massCons'  , ...
                     'ne_eac_chargeCons', ...
                     'pe_eac_massCons'  , ...
                     'pe_eac_chargeCons', ...
                     'ne_cc_chargeCons' , ...
                     'pe_cc_chargeCons' , ...
                     'EIeq', ...
                     'controlEq'};
                switch ctrltype
                    case 'I'
                        types = {'cell','cell','cell','cell', 'cell','cell','cell','cell', 'cell', 'cntrl'};
                    case 'E'
                        types = {'cell','cell','cell','cell', 'cell','cell','cell','cell', 'cntrl', 'cntrl'};
                        order = [1:8,10,9];
                        eqs = {eqs{order}};
                        names = {names{order}};
                    otherwise
                        error()
                end
            end
            primaryVars = model.getPrimaryVariables();

            %% setup LinearizedProblem that can be processed by MRST Newton API
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
            
        end

        function state = updateTemperature(model, state)
        % Dispatch the temperature in all the submodels

            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            eac   = 'ElectrodeActiveComponent';
            cc    = 'CurrentCollector';
            
            % (here we assume that the ThermalModel has the "parent" grid)
            state.(elyte).T    = state.T(model.(elyte).G.mappings.cellmap);
            state.(ne).(eac).T = state.T(model.(ne).(eac).G.mappings.cellmap);
            state.(ne).(cc).T  = state.T(model.(ne).(cc).G.mappings.cellmap);
            state.(pe).(eac).T = state.T(model.(pe).(eac).G.mappings.cellmap);
            state.(pe).(cc).T  = state.T(model.(pe).(cc).G.mappings.cellmap);
            
            % Update temperature in the active materials of the electrodes.
            state.(ne).(eac) = model.(ne).(eac).updateTemperature(state.(ne).(eac));
            state.(pe).(eac) = model.(pe).(eac).updateTemperature(state.(pe).(eac));
            
        end
        
        
        function state = updateElectrolyteCoupling(model, state)
        % Assemble the electrolyte coupling by adding the ion sources from the electrodes
            
            battery = model;
            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            am    = 'ActiveMaterial';
            eac   = 'ElectrodeActiveComponent';
            
            vols = battery.(elyte).G.cells.volumes;
            F = battery.con.F;
            
            ccSourceName = battery.(elyte).chargeCarrierSourceName;
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
            
            ne_R = state.(ne).(eac).(am).R;
            coupterm = getCoupTerm(couplingterms, 'NegativeElectrode-Electrolyte', coupnames);
            elytecells = coupterm.couplingcells(:, 2);
            elyte_c_source(elytecells) = ne_R.*vols(elytecells); % we divide with F later
            
            pe_R = state.(pe).(eac).(am).R;
            coupterm = getCoupTerm(couplingterms, 'PositiveElectrode-Electrolyte', coupnames);
            elytecells = coupterm.couplingcells(:, 2);
            elyte_c_source(elytecells) = pe_R.*vols(elytecells);
            
            elyte_e_source = elyte_c_source.*battery.(elyte).sp.z(1)*F; 
            
            state.Electrolyte.(ccSourceName) = elyte_c_source; 
            state.Electrolyte.eSource = elyte_e_source;
            
        end
        
        function state = updateAccumTerms(model, state, state0, dt)
        % Assemble the accumulation terms for transport equations (in electrolyte and electrodes)
            
            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            am    = 'ActiveMaterial';
            eac   = 'ElectrodeActiveComponent';
            
            ccAccumName = model.(elyte).chargeCarrierAccumName;
            
            cdotcc  = (state.(elyte).cs{1} - state0.(elyte).cs{1})/dt;
            effectiveVolumes = model.(elyte).volumeFraction.*model.(elyte).G.cells.volumes;
            ccAccum  = effectiveVolumes.*cdotcc;
            state.(elyte).(ccAccumName) = ccAccum;
            
            names = {ne, pe};
            for i = 1 : numel(names)
                elde = names{i}; % electrode name
                cdotcc   = (state.(elde).(eac).c - state0.(elde).(eac).c)/dt;
                effectiveVolumes = model.(elde).(eac).volumeFraction.*model.(elde).(eac).G.cells.volumes;
                ccAccum  = effectiveVolumes.*cdotcc;
                state.(elde).(eac).(ccAccumName) = ccAccum;
            end
            
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
            eac   = 'ElectrodeActiveComponent';
            cc    = 'CurrentCollector';
            
            eldes = {ne, pe};
            phi_elyte = state.(elyte).phi;
            c_elyte = state.(elyte).cs{1};
            
            elyte_cells = zeros(model.G.cells.num, 1);
            elyte_cells(bat.(elyte).G.mappings.cellmap) = (1 : bat.(elyte).G.cells.num)';

            for ind = 1 : numel(eldes)
                elde = eldes{ind};
                state.(elde).(eac).(am).phiElectrolyte = phi_elyte(elyte_cells(bat.(elde).(eac).G.mappings.cellmap));
                state.(elde).(eac).(am).cElectrolyte = c_elyte(elyte_cells(bat.(elde).(eac).G.mappings.cellmap));
            end
            
        end

        function state = setupExternalCouplingNegativeElectrode(model, state)
        %
        % Setup external electronic coupling of the negative electrode at the current collector
        %
            ne = 'NegativeElectrode';
            cc = 'CurrentCollector';
            
            phi = state.(ne).(cc).phi;

            jExternal = model.(ne).(cc).setupExternalCoupling(phi, 0);
            
            state.(ne).(cc).jExternal = jExternal;
            
        end
        
        function state = setupExternalCouplingPositiveElectrode(model, state)
        %
        % Setup external electronic coupling of the positive electrode at the current collector
        %            
            pe = 'PositiveElectrode';
            cc = 'CurrentCollector';
            
            phi = state.(pe).(cc).phi;
            E = state.(pe).(cc).E;
            
            jExternal = model.(pe).(cc).setupExternalCoupling(phi, E);
            
            state.(pe).(cc).jExternal = jExternal;
            
        end

        function state = setupEIEquation(model, state)
            
            pe = 'PositiveElectrode';
            cc = 'CurrentCollector';
            
            I = state.(pe).(cc).I;
            E = state.(pe).(cc).E;
            phi = state.(pe).(cc).phi;
            
            coupterm = model.(pe).(cc).couplingTerm;
            faces = coupterm.couplingfaces;
            cond_pcc = model.(pe).(cc).EffectiveElectricalConductivity;
            [trans_pcc, cells] = model.(pe).(cc).operators.harmFaceBC(cond_pcc, faces);
            state.EIeq = sum(trans_pcc.*(state.(pe).(cc).phi(cells) - E)) - I;

        end
        
        
        function state = initStateAD(model,state)
            
            bat = model;
            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            am    = 'ActiveMaterial';
            eac   = 'ElectrodeActiveComponent';
            cc    = 'CurrentCollector';
            
            adbackend = model.AutoDiffBackend();
            if(model.use_solid_diffusion)
                [state.(elyte).cs{1}  , ...
                    state.(elyte).phi    , ...
                    state.(ne).(eac).c   , ...
                    state.(ne).(eac).phi , ...
                    state.(ne).(eac).(am).cElectrode , ...
                    state.(pe).(eac).c   , ...
                    state.(pe).(eac).phi , ...
                    state.(pe).(eac).(am).cElectrode , ...
                    state.(ne).(cc).phi  , ...
                    state.(pe).(cc).phi  , ...
                    state.(pe).(cc).E, ...
                    state.(pe).(cc).I] = ...
                    adbackend.initVariablesAD(...
                    state.(elyte).cs{1}  , ...
                    state.(elyte).phi    , ...
                    state.(ne).(eac).c   , ...
                    state.(ne).(eac).phi , ...
                    state.(ne).(eac).(am).cElectrode , ...
                    state.(pe).(eac).c   , ...
                    state.(pe).(eac).phi , ...
                    state.(pe).(eac).(am).cElectrode , ...
                    state.(ne).(cc).phi  , ...
                    state.(pe).(cc).phi  , ...
                    state.(pe).(cc).E, ...
                    state.(pe).(cc).I);
            else
                [state.(elyte).cs{1}  , ...
                    state.(elyte).phi    , ...
                    state.(ne).(eac).c   , ...
                    state.(ne).(eac).phi , ...
                    state.(pe).(eac).c   , ...
                    state.(pe).(eac).phi , ...
                    state.(ne).(cc).phi  , ...
                    state.(pe).(cc).phi  , ...
                    state.(pe).(cc).E, ...
                    state.(pe).(cc).I] = ...
                    adbackend.initVariablesAD(...
                    state.(elyte).cs{1}  , ...
                    state.(elyte).phi    , ...
                    state.(ne).(eac).c   , ...
                    state.(ne).(eac).phi , ...
                    state.(pe).(eac).c   , ...
                    state.(pe).(eac).phi , ...
                    state.(ne).(cc).phi  , ...
                    state.(pe).(cc).phi  , ...
                    state.(pe).(cc).E, ...
                    state.(pe).(cc).I);
            end
            % PRIMARY variables
        end
        
        function p = getPrimaryVariables(model)
            
            bat = model;
            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            am    = 'ActiveMaterial';
            eac   = 'ElectrodeActiveComponent';
            cc    = 'CurrentCollector';
            if(model.use_solid_diffusion)
            p = {{elyte, 'cs', 1} , ...
                 {elyte, 'phi'}   , ...   
                 {ne, eac, 'c'}   , ...    
                 {ne, eac, 'phi'} , ...   
                 {ne, eac, am, 'cElectrode'} , ...
                 {pe, eac, 'c'}   , ...    
                 {pe, eac, 'phi'} , ...   
                 {pe, eac, am, 'cElectrode'} , ...
                 {ne, cc, 'phi'}  , ...    
                 {pe, cc, 'phi'}  , ...
                 {pe, cc, 'E'}, ...
                 {pe, cc, 'I'}};
            else
                 p = {{elyte, 'cs', 1} , ...
                 {elyte, 'phi'}   , ...   
                 {ne, eac, 'c'}   , ...    
                 {ne, eac, 'phi'} , ...   
                 {pe, eac, 'c'}   , ...    
                 {pe, eac, 'phi'} , ...   
                 {ne, cc, 'phi'}  , ...    
                 {pe, cc, 'phi'}  , ...
                 {pe, cc, 'E'}, ...
                 {pe, cc, 'I'}};
            end
            
        end

        function scales = getScales(model)
            prims = model.getPrimaryVariables();
            scales = struct();
            for i=1:numel(prims)
                p = prims{i};
                scale = struct('min',-inf,'max',inf', ...
                               'relchangemax',1,'abschangemax',inf);
                if( strcmp(p{end},'c'))
                    % assumes consentration
                    cmin = model.cmin;
                    mname=p;
                    mname{end} = 'ActiveMaterial';
                    mname{end+1}= 'Li';
                    mname{end+1} = 'cmax';
                    cmax = model.getProp(model,mname);
                    scale = struct('min',0.001*cmax,...
                                   'max',cmax, ...
                                   'relchangemax',inf, ...
                                   'abschangemax',0.1*cmax);
                elseif(strcmp(p{end},'cElectrode'))
                    % assumes consentration
                    cmin = model.cmin;
                    mname=p;
                    mname{end}= 'Li';
                    mname{end+1} = 'cmax';
                    cmax = model.getProp(model,mname);
                    scale = struct('min',0.001*cmax,...
                                   'max',cmax, ...
                                   'relchangemax',inf,'abschangemax',0.1*cmax);
    
                elseif(isnumeric(p{end}) && strcmp(p{1},'Electrolyte'))
                    cmin = model.cmin;
                    mname=p;
                    %cmax = model.getProp(model,mname);
                    scale = struct('min',100,'max',2000, ...
                               'relchangemax',inf,'abschangemax',100);
                elseif(strcmp(p{end},'phi') || strcmp(p{end},'E'))
                    scale = struct('min',-10,'max',10, ...
                               'relchangemax',inf,'abschangemax',0.2);
                elseif(strcmp(p{end},'I'))
                    scale = struct('min',-inf,'max',inf, ...
                               'relchangemax',inf,'abschangemax',inf);
                else
                    error();
                end
                   
                scales = model.setNewProp(scales,p,scale);
                
            end
            %            scales =[];
        end
        
        
        function validforces = getValidDrivingForces(model)
        
            validforces=struct('src', [], 'stopFunction', []); 
            
        end
        
        function model = validateModel(model, varargin)
            mnames = {{'Electrolyte'}, ...
                      {'PositiveElectrode','ElectrodeActiveComponent'}, ...
                      {'NegativeElectrode','ElectrodeActiveComponent'}, ...
                      {'NegativeElectrode','CurrentCollector'}, ...
                      {'PositiveElectrode','CurrentCollector'}};
            model.Electrolyte.AutoDiffBackend=model.AutoDiffBackend;
            model.Electrolyte=model.Electrolyte.validateModel(varargin{:});
            model.PositiveElectrode.ElectrodeActiveComponent.AutoDiffBackend= model.AutoDiffBackend;
            model.PositiveElectrode.ElectrodeActiveComponent = model.PositiveElectrode.ElectrodeActiveComponent.validateModel(varargin{:});
            model.NegativeElectrode.ElectrodeActiveComponent.AutoDiffBackend= model.AutoDiffBackend;
            model.NegativeElectrode.ElectrodeActiveComponent = model.NegativeElectrode.ElectrodeActiveComponent.validateModel(varargin{:});
            model.NegativeElectrode.CurrentCollector.AutoDiffBackend= model.AutoDiffBackend;
            model.NegativeElectrode.CurrentCollector=model.NegativeElectrode.CurrentCollector.validateModel(varargin{:});
            model.PositiveElectrode.CurrentCollector.AutoDiffBackend=model.AutoDiffBackend;
            model.PositiveElectrode.CurrentCollector= model.PositiveElectrode.CurrentCollector.validateModel(varargin{:});
        
        end
        

        function [state, report] = updateState(model, state, problem, dx, drivingForces)

            [state, report] = updateState@BaseModel(model, state, problem, dx, drivingForces);
            
            %% cap concentrations
            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            eac   = 'ElectrodeActiveComponent';
            am    = 'ActiveMaterial';

            cmin = model.cmin;
            
            state.(elyte).cs{1} = max(cmin, state.(elyte).cs{1});
            
            eldes = {ne, pe};
            for ind = 1 : numel(eldes)
                elde = eldes{ind};
                state.(elde).(eac).c = max(cmin, state.(elde).(eac).c);
                cmax = model.(elde).(eac).(am).Li.cmax;
                state.(elde).(eac).c = min(cmax, state.(elde).(eac).c);
            end
            
            report = [];
            
        end
        
    end
    
end
