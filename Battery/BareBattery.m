classdef BareBattery < BaseModel
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
        
        
    end
    
    methods
        
        function model = BareBattery(paramobj)

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
            am    = 'ActiveMaterial';
            cmax_ne = model.(ne).(am).Li.cmax;
            cmax_pe = model.(pe).(am).Li.cmax;
            model.cmin = 1e-5*max(cmax_ne, cmax_pe);
            
        end

        
        function electrode = setupElectrode(model, paramobj)
        % Setup the electrode models (both :attr:`NegativeElectrode` and :attr:`PositiveElectrode`). Here, :code:`paramobj`
        % is instance of :class:`ElectrodeInputParams <Electrochemistry.Electrodes.ElectrodeInputParams>`
            electrode = ElectrodeActiveComponent(paramobj);
            
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
            elyte   = 'Electrolyte';
            
            G_elyte = model.(elyte).G;
            elytecelltbl.cells = (1 : G_elyte.cells.num)';
            elytecelltbl.globalcells = G_elyte.mappings.cellmap;
            elytecelltbl = IndexArray(elytecelltbl);

            eldes = {ne, pe};

            for ind = 1 : numel(eldes)

                elde = eldes{ind};
                G_elde  = model.(elde).G;
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
            sep   = 'Separator';

            elyte_cells = zeros(model.G.cells.num, 1);
            elyte_cells(model.(elyte).G.mappings.cellmap) = (1 : model.(elyte).G.cells.num)';

            model.(elyte).volumeFraction = NaN(model.(elyte).G.cells.num, 1);
            model.(elyte).volumeFraction(elyte_cells(model.(ne).G.mappings.cellmap))  = model.(ne).porosity;
            model.(elyte).volumeFraction(elyte_cells(model.(pe).G.mappings.cellmap))  = model.(pe).porosity;
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
                        
            %% synchronize temperatures
            initstate = model.updateTemperature(initstate);
            
            %% setup initial NegativeElectrode state
            
            % shortcut
            % negAm : ActiveMaterial of the negative electrode
            
            negAm = bat.(ne).(am); 
            
            m = (1 ./ (negAm.theta100 - negAm.theta0));
            b = -m .* negAm.theta0;
            theta = (SOC - b) ./ m;
            c = theta .* negAm.Li.cmax;
            c = c*ones(negAm.G.cells.num, 1);

            initstate.(ne).c = c;
            % We bypass the solid diffusion equation to set directly the particle surface concentration (this is a bit hacky)
            initstate.(ne).(am).cElectrode = c;
            initstate.(ne).(am) = negAm.updateOCP(initstate.(ne).(am));

            OCP = initstate.(ne).(am).OCP;
            ref = OCP(1);
            
            initstate.(ne).phi = OCP - ref;

            %% setup initial PositiveElectrode state

            % shortcut
            % posAm : ActiveMaterial of the positive electrode
            
            posAm = bat.(pe).(am);
            
            m = (1 ./ (posAm.theta100 - posAm.theta0));
            b = -m .* posAm.theta0;
            theta = (SOC - b) ./ m;
            c = theta .* posAm.Li.cmax;
            c = c*ones(posAm.G.cells.num, 1);

            initstate.(pe).c = c;
            % We bypass the solid diffusion equation to set directly the particle surface concentration (this is a bit hacky)
            initstate.(pe).(am).cElectrode = c;
            initstate.(pe).(am) = posAm.updateOCP(initstate.(pe).(am));
            
            OCP = initstate.(pe).(am).OCP;
            initstate.(pe).phi = OCP - ref;

            %% setup initial Electrolyte state

            initstate.(elyte).phi = zeros(bat.(elyte).G.cells.num, 1)-ref;
            cs = cell(2,1);
            initstate.(elyte).cs = cs;
            initstate.(elyte).cs{1} = 1000*ones(bat.(elyte).G.cells.num, 1);

            %% setup initial Current collectors state

            initstate.(pe).E = OCP(1) - ref;
            initstate.(pe).I = 0;
            
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
            elyte   = 'Electrolyte';
            am      = 'ActiveMaterial';
            
            electrodes = {ne, pe};

            %% Synchronization across components

            % temperature
            state = battery.updateTemperature(state);

            state.(elyte) = battery.(elyte).updateConcentrations(state.(elyte));
            
            for ind = 1 : numel(electrodes)
                elde = electrodes{ind};
                % potential and concentration between active material and electode active component
                state.(elde) = battery.(elde).updatePhi(state.(elde));
                state.(elde) = battery.(elde).updateChargeCarrier(state.(elde));
                
            end
            
            %% Accumulation terms

            state = battery.updateAccumTerms(state, state0, dt);

            %% Update Electrolyte -> Electrodes coupling 
            
            state = battery.updateElectrodeCoupling(state); 

            %% Update reaction rates in both electrodes

            for ind = 1 : numel(electrodes)
                elde = electrodes{ind};
                state.(elde).(am) = battery.(elde).(am).updateReactionRateCoefficient(state.(elde).(am));
                state.(elde).(am) = battery.(elde).(am).updateDiffusionCoefficient(state.(elde).(am));
                state.(elde).(am) = battery.(elde).(am).updateOCP(state.(elde).(am));
                state.(elde).(am) = battery.(elde).(am).updateReactionRate(state.(elde).(am));
            end

            %% Update Electrodes -> Electrolyte  coupling

            state = battery.updateElectrolyteCoupling(state);
            
            %% Update  external coupling
            
            state = model.setupExternalCouplingNegativeElectrode(state);
            state = model.setupExternalCouplingPositiveElectrode(state);
            
            %% elyte charge conservation

            state.(elyte) = battery.(elyte).updateCurrentBcSource(state.(elyte));
            state.(elyte) = battery.(elyte).updateConductivity(state.(elyte));
            state.(elyte) = battery.(elyte).updateChemicalCurrent(state.(elyte));
            state.(elyte) = battery.(elyte).updateCurrent(state.(elyte));
            state.(elyte) = battery.(elyte).updateChargeConservation(state.(elyte));

            %% Electrodes charge conservation - Active material part

            for ind = 1 : numel(electrodes)
                elde = electrodes{ind};
                state.(elde) = battery.(elde).updateIonAndCurrentSource(state.(elde));
                state.(elde) = battery.(elde).updateCurrent(state.(elde));
                state.(elde) = battery.(elde).updateChargeConservation(state.(elde));
            end
            
            %% elyte mass conservation

            state.(elyte) = battery.(elyte).updateDiffusionCoefficient(state.(elyte));
            state.(elyte) = battery.(elyte).updateChargeCarrierFlux(state.(elyte));
            state.(elyte) = battery.(elyte).updateMassConservation(state.(elyte));
            
            for ind = 1 : numel(electrodes)
                elde = electrodes{ind};
                
                %% Electrodes mass conservation
                state.(elde) = battery.(elde).updateChargeCarrierFlux(state.(elde));
                state.(elde) = battery.(elde).updateMassConservation(state.(elde));
                
            end

            %% update solid diffustion equations
            for ind = 1 : numel(electrodes)
                elde = electrodes{ind};
                state.(elde).(am) = battery.(elde).(am).assembleSolidDiffusionEquation(state.(elde).(am));
            end
            
            %% setup relation between E and I at positive current collectror
            
            state = model.setupEIEquation(state);
            
            %% Set up the governing equations
            
            eqs = {};
            
            %% We collect mass and charge conservation equations for the electrolyte and the electrodes

            eqs{end + 1} = state.(elyte).massCons;
            eqs{end + 1} = state.(elyte).chargeCons;
            
            eqs{end + 1} = state.(ne).massCons;
            eqs{end + 1} = state.(ne).chargeCons;
            eqs{end + 1} = state.(ne).(am).solidDiffusionEq;
            
            eqs{end + 1} = state.(pe).massCons;
            eqs{end + 1} = state.(pe).chargeCons;
            eqs{end + 1} = state.(pe).(am).solidDiffusionEq;
            
            eqs{end + 1} = state.EIeq;
            
            % we add the control equation
            val = drivingForces.src(time);
            eqs{end + 1} = state.(pe).I - val;

            %% Give type and names to equations and names of the primary variables (for book-keeping)
            
            types = {'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell', 'cell'};
            
            names = {'elyte_massCons'   , ...
                     'elyte_chargeCons' , ...
                     'ne_massCons'      , ...
                     'ne_chargeCons'    , ...
                     'ne_am_soliddiffeq', ...
                     'pe_massCons'      , ...
                     'pe_chargeCons'    , ...
                     'pe_am_soliddiffeq', ...
                     'EIeq', ...
                     'controlEq'};
            
            primaryVars = model.getPrimaryVariables();

            %% setup LinearizedProblem that can be processed by MRST Newton API
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
            
        end

        function state = updateTemperature(model, state)
        % Dispatch the temperature in all the submodels

            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            cc    = 'CurrentCollector';
            
            % (here we assume that the ThermalModel has the "parent" grid)
            state.(elyte).T    = state.T(model.(elyte).G.mappings.cellmap);
            state.(ne).T = state.T(model.(ne).G.mappings.cellmap);
            state.(pe).T = state.T(model.(pe).G.mappings.cellmap);
            
            % Update temperature in the active materials of the electrodes.
            state.(ne) = model.(ne).updateTemperature(state.(ne));
            state.(pe) = model.(pe).updateTemperature(state.(pe));
            
        end
        
        
        function state = updateElectrolyteCoupling(model, state)
        % Assemble the electrolyte coupling by adding the ion sources from the electrodes
            
            battery = model;
            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            am    = 'ActiveMaterial';
            
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
            
            ne_R = state.(ne).(am).R;
            coupterm = getCoupTerm(couplingterms, 'NegativeElectrode-Electrolyte', coupnames);
            elytecells = coupterm.couplingcells(:, 2);
            elyte_c_source(elytecells) = ne_R.*vols(elytecells); % we divide with F later
            
            pe_R = state.(pe).(am).R;
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
            
            ccAccumName = model.(elyte).chargeCarrierAccumName;
            
            cdotcc  = (state.(elyte).cs{1} - state0.(elyte).cs{1})/dt;
            effectiveVolumes = model.(elyte).volumeFraction.*model.(elyte).G.cells.volumes;
            ccAccum  = effectiveVolumes.*cdotcc;
            state.(elyte).(ccAccumName) = ccAccum;
            
            names = {ne, pe};
            for i = 1 : numel(names)
                elde = names{i}; % electrode name
                cdotcc   = (state.(elde).c - state0.(elde).c)/dt;
                effectiveVolumes = model.(elde).volumeFraction.*model.(elde).G.cells.volumes;
                ccAccum  = effectiveVolumes.*cdotcc;
                state.(elde).(ccAccumName) = ccAccum;
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
            
            eldes = {ne, pe};
            phi_elyte = state.(elyte).phi;
            c_elyte = state.(elyte).cs{1};
            
            elyte_cells = zeros(model.G.cells.num, 1);
            elyte_cells(bat.(elyte).G.mappings.cellmap) = (1 : bat.(elyte).G.cells.num)';

            for ind = 1 : numel(eldes)
                elde = eldes{ind};
                state.(elde).(am).phiElectrolyte = phi_elyte(elyte_cells(bat.(elde).G.mappings.cellmap));
                state.(elde).(am).cElectrolyte = c_elyte(elyte_cells(bat.(elde).G.mappings.cellmap));
            end
            
        end

        function state = setupExternalCouplingNegativeElectrode(model, state)
        %
        % Setup external electronic coupling of the negative electrode
        %
         
            battery = model;
            ne = 'NegativeElectrode';

            phi = state.(ne).phi;
            
            couplingterms = battery.couplingTerms;
            coupnames = battery.couplingNames;
            coupterm = getCoupTerm(couplingterms, 'Exterior-NegativeElectrode', coupnames);
            
            jBcSource = phi*0.0; %NB hack to initialize zero ad
            sigmaeff = model.(ne).EffectiveElectricalConductivity;
            faces = coupterm.couplingfaces;
            % We impose potential equal to zero at negative electrode
            bcval = 0;
            [t, cells] = model.(ne).operators.harmFaceBC(sigmaeff, faces);
            jBcSource(cells) = jBcSource(cells) + t.*(bcval - phi(cells));
            
            state.(ne).jBcSource = jBcSource;
            
        end
        
        function state = setupExternalCouplingPositiveElectrode(model, state)
        %
        % Setup external electronic coupling of the positive electrode at the current collector
        %            
         
            battery = model;
            pe = 'PositiveElectrode';
            
            phi = state.(pe).phi;
            E = state.(pe).E;
            
            couplingterms = battery.couplingTerms;
            coupnames = battery.couplingNames;
            coupterm = getCoupTerm(couplingterms, 'Exterior-PositiveElectrode', coupnames);
            
            jBcSource = phi*0.0; %NB hack to initialize zero ad
            sigmaeff = model.(pe).EffectiveElectricalConductivity;
            faces = coupterm.couplingfaces;
            % We impose potential equal to value given by E at the positive electrode
            bcval = E;
            [t, cells] = model.(pe).operators.harmFaceBC(sigmaeff, faces);
            jBcSource(cells) = jBcSource(cells) + t.*(bcval - phi(cells));
            
            state.(pe).jBcSource = jBcSource;
        
        end

        function state = setupEIEquation(model, state)
            
            battery = model;
            pe = 'PositiveElectrode';
            
            couplingterms = battery.couplingTerms;
            coupnames = battery.couplingNames;
            coupterm = getCoupTerm(couplingterms, 'Exterior-PositiveElectrode', coupnames);
            
            I = state.(pe).I;
            E = state.(pe).E;
            phi = state.(pe).phi;
            
            faces = coupterm.couplingfaces;
            cond_pcc = model.(pe).EffectiveElectricalConductivity;
            [trans_pcc, cells] = model.(pe).operators.harmFaceBC(cond_pcc, faces);
            state.EIeq = sum(trans_pcc.*(state.(pe).phi(cells) - E)) - I;

        end
        
        
        function state = initStateAD(model,state)
            
            bat = model;
            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            am    = 'ActiveMaterial';
            
            adbackend = model.AutoDiffBackend();
            useMex=false;
            if(isprop(adbackend, 'useMex'))
                useMex = adbackend.useMex; 
            end
            opts=struct('types',[1,1,2,2,3,3,4,5,6,7,8],'useMex',useMex);
            [state.(elyte).cs{1}  , ...
             state.(elyte).phi    , ...   
             state.(ne).c   , ...   
             state.(ne).phi , ...   
             state.(ne).(am).cElectrode , ...
             state.(pe).c   , ...    
             state.(pe).phi , ...   
             state.(pe).(am).cElectrode , ...
             state.(pe).E, ...
             state.(pe).I] = ...
                adbackend.initVariablesAD(...
                    state.(elyte).cs{1}  , ...
                    state.(elyte).phi    , ...   
                    state.(ne).c   , ...    
                    state.(ne).phi , ...   
                    state.(ne).(am).cElectrode , ...
                    state.(pe).c   , ...    
                    state.(pe).phi , ...   
                    state.(pe).(am).cElectrode , ...
                    state.(pe).E, ...
                    state.(pe).I, ...
                    opts); 
            % PRIMARY variables
        end
        
        function p = getPrimaryVariables(model)
            
            bat = model;
            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            am    = 'ActiveMaterial';
            p = {{elyte, 'cs', 1} , ...
                 {elyte, 'phi'}   , ...   
                 {ne, 'c'}   , ...    
                 {ne, 'phi'} , ...   
                 {ne, am, 'cElectrode'} , ...
                 {pe, 'c'}   , ...    
                 {pe, 'phi'} , ...   
                 {pe, am, 'cElectrode'} , ...
                 {pe, 'E'}, ...
                 {pe, 'I'}};
            
        end
        
        
        function validforces = getValidDrivingForces(model)
        
            validforces=struct('src', [], 'stopFunction', []); 
            
        end
        
        function model = validateModel(model, varargin)

            model.Electrolyte.AutoDiffBackend = model.AutoDiffBackend;
            model.Electrolyte = model.Electrolyte.validateModel(varargin{:});
            
            model.PositiveElectrode.AutoDiffBackend = model.AutoDiffBackend;
            model.PositiveElectrode                 = model.PositiveElectrode.validateModel(varargin{:});
            model.NegativeElectrode.AutoDiffBackend = model.AutoDiffBackend;
            model.NegativeElectrode                 = model.NegativeElectrode.validateModel(varargin{:});
        
        end
        

        function [state, report] = updateState(model, state, problem, dx, drivingForces)

            [state, report] = updateState@BaseModel(model, state, problem, dx, drivingForces);
            
            %% cap concentrations
            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            am    = 'ActiveMaterial';

            cmin = model.cmin;
            
            state.(elyte).cs{1} = max(cmin, state.(elyte).cs{1});
            
            eldes = {ne, pe};
            for ind = 1 : numel(eldes)
                elde = eldes{ind};
                state.(elde).c = max(cmin, state.(elde).c);
                cmax = model.(elde).(am).Li.cmax;
                state.(elde).c = min(cmax, state.(elde).c);
            end
            
            report = [];
            
        end
        
    end
    
end
