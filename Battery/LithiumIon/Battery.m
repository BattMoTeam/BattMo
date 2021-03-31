classdef Battery < PhysicalModel

    properties
        
        con = PhysicalConstants();

        % Temperature and SOC
        % for the moment here, for convenience. Will be moved
        T
        SOC

        % Input current
        J
        % Voltage cut
        Ucut

        % Components
        Electrolyte
        NegativeElectrode
        PositiveElectrode
        
        couplingTerms
        couplingNames
        
    end
    
    methods
        
        function model = Battery(paramobj)
        % Shorcuts used here
        % elyte : Electrolyte
        % ne : NegativeElectrode
        % pe : PositiveElectrode
            
            model = model@PhysicalModel([]);
            
            % OBS : All the submodels should have same backend (this is not assigned automaticallly for the moment)
            model.AutoDiffBackend = SparseAutoDiffBackend('useBlocks', true);
            
            %% Setup the model using the input parameters
            fdnames = {'G', ...
                       'couplingTerms', ...
                       'T'  , ...
                       'SOC', ...
                       'J'  , ...
                       'Ucut'};
            
            model = dispatchParams(model, paramobj, fdnames);
            
            % Assign the components : Electrolyte, NegativeElectrode, PositiveElectrode
            model.Electrolyte       = Electrolyte(paramobj.elyte);
            model.NegativeElectrode = Electrode(paramobj.ne);
            model.PositiveElectrode = Electrode(paramobj.pe);

            % setup couplingNames
            model.couplingNames = cellfun(@(x) x.name, model.couplingTerms, 'uniformoutput', false);
            
            % setup Electrolyte porosity
            model = model.setElectrolytePorosity();
            
        end
        
        function [problem, state] = getEquations(model, state0, state,dt, drivingForces, varargin)
            
            time = state0.time + dt;
            state = model.initStateAD(state);
            
            %% for now temperature and SOC are kept constant
            nc = model.G.cells.num;
            state.T   = model.T*ones(nc, 1);
            state.SOC = model.SOC*ones(nc, 1);
            
            % Shortcuts
            battery = model;
            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            eac     = 'ElectrodeActiveComponent';
            cc      = 'CurrentCollector';
            elyte   = 'Electrolyte';
            am      = 'ActiveMaterial';
            
            electrodes = {'NegativeElectrode', 'PositiveElectrode'};

            %% Synchronization across components

            % temperature
            state = battery.updateT(state);
            
            
            for ind = 1 : numel(electrodes)
                elde = electrodes{ind};
                
                % temperature
                state.(elde) = battery.(elde).updateT(state.(elde));
                state.(elde).(eac) = battery.(elde).(eac).updateT(state.(elde).(eac));

                %% potential and concentration between active material and electode active component
                state.(elde).(eac) = battery.(elde).(eac).updatePhi(state.(elde).(eac));
                state.(elde).(eac) = battery.(elde).(eac).updateChargeCarrier(state.(elde).(eac));
                
            end
            
            %% Accumulation terms

            state = battery.updateAccumTerms(state, state0, dt);

            %% Update Electrolyte -> Electrodes coupling 
            
            state = battery.setupElectrodeCoupling(state); 

            %% Update reaction rates in both electrodes

            for ind = 1 : numel(electrodes)
                elde = electrodes{ind};
                state.(elde).(eac).(am) = battery.(elde).(eac).(am).updateMaterialProperties(state.(elde).(eac).(am));
                state.(elde).(eac).(am) = battery.(elde).(eac).(am).updateReactionRate(state.(elde).(eac).(am));
            end

            %% Update Electrodes -> Electrolyte  coupling

            state = battery.setupElectrolyteCoupling(state);
            
            %% Update Current collectors <-> Electrode active components couplings
            
            state.(ne) = battery.(ne).setupCoupling(state.(ne));
            state.(pe) = battery.(pe).setupCoupling(state.(pe));

            %% elyte charge conservation

            state.(elyte) = battery.(elyte).updateCurrentBcSource(state.(elyte));
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
                state.(elde).(eac) = battery.(elde).(eac).updateDiffusionCoefficient(state.(elde).(eac));
                state.(elde).(eac) = battery.(elde).(eac).updateChargeCarrierFlux(state.(elde).(eac));
                state.(elde).(eac) = battery.(elde).(eac).updateMassConservation(state.(elde).(eac));
                
                %% Electrodes charge conservation - current collector part
                state.(elde).(cc) = battery.(elde).(cc).updateCurrent(state.(elde).(cc));
                state.(elde).(cc) = battery.(elde).(cc).updateChargeConservation(state.(elde).(cc));
            
            end
            
            %% Set up the governing equations
            
            eqs={};
            
            %% We collect mass and charge conservation equations for the electrolyte and the electrodes

            eqs{end + 1} = state.(elyte).massCons;
            eqs{end + 1} = state.(elyte).chargeCons;
            
            eqs{end + 1} = state.(ne).(eac).massCons;
            eqs{end + 1} = state.(ne).(eac).chargeCons;
            eqs{end + 1} = state.(ne).(cc).chargeCons;
            
            eqs{end + 1} = state.(pe).(eac).massCons;
            eqs{end + 1} = state.(pe).(eac).chargeCons;
            eqs{end + 1} = state.(pe).(cc).chargeCons;
            
            %% We setup and add the control equation (fixed total current at PositiveCurrentCollector)
            
            src = drivingForces.src(time);
            coupterms = battery.(ne).couplingTerms;
            coupnames = battery.(ne).couplingNames;
            coupterm = getCoupTerm(coupterms, 'bc-CurrentCollector', coupnames);
            faces = coupterm.couplingfaces;
            bcval = state.(pe).(cc).E;
            cond_pcc = battery.(pe).(cc).EffectiveElectronicConductivity;
            [trans_pcc, cells] = battery.(pe).(cc).operators.harmFaceBC(cond_pcc, faces);
            control = src - sum(trans_pcc.*(bcval - state.(pe).(cc).phi(cells)));
            
            eqs{end + 1} = -control;

            %% Give type and names to equations and names of the primary variables (for book-keeping)
            
            types = {'cell','cell','cell','cell',...
                     'cell','cell','cell','cell','cell'};
            
            names = {'elyte_massCons'   , ...
                     'elyte_chargeCons' , ...
                     'ne_eac_massCons'  , ...
                     'ne_eac_chargeCons', ...
                     'ne_cc_chargeCons' , ...
                     'ne_eac_massCons'  , ...
                     'ne_eac_chargeCons', ...
                     'ne_cc_chargeCons' , ...
                     'control'};
            
            primaryVars = model.getPrimaryVariables();

            %% setup LinearizedProblem that can be processed by MRST Newton API
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
        
        end

        
        function state = updateT(model, state)
            names = {'NegativeElectrode', 'PositiveElectrode', 'Electrolyte'};
            for ind = 1 : numel(names);
                name = names{ind};
                nc = model.(name).G.cells.num;
                state.(name).T = state.T(1)*ones(nc, 1);
            end
        end
        
        function model = setElectrolytePorosity(model)
        % Abbreviations used in this function 
        % elyte : Electrolyte
        % ne    : NegativeElectrode
        % pe    : PositiveElectrode
            
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
        % Setup initial state
        %
        % Abbreviations used in this function 
        % elyte : Electrolyte
        % ne    : NegativeElectrode
        % pe    : PositiveElectrode
        % eac   : ElectrodeActiveComponent
        % cc    : CurrentCollector
            
            nc = model.G.cells.num;

            SOC = model.SOC;
            T   = model.T;
            
            initstate.T   =  T*ones(nc, 1);
            initstate.SOC =  SOC*ones(nc, 1);
            
            bat = model;
            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            am    = 'ActiveMaterial';
            eac   = 'ElectrodeActiveComponent';
            cc    = 'CurrentCollector';
            
            %% synchronize temperatures
            initstate = model.updateT(initstate);
            initstate.(ne) = bat.(ne).updateT(initstate.(ne));
            initstate.(ne).(eac) = bat.(ne).(eac).updateT(initstate.(ne).(eac));
            initstate.(pe) = bat.(pe).updateT(initstate.(pe));
            initstate.(pe).(eac) = bat.(pe).(eac).updateT(initstate.(pe).(eac));
            
            %% setup initial NegativeElectrode state
            
            % shortcut
            % negAm : ActiveMaterial of the negative electrode
            
            negAm = bat.(ne).(eac).(am); 
            
            m = (1 ./ (negAm.theta100 - negAm.theta0));
            b = -m .* negAm.theta0;
            theta = (SOC - b) ./ m;
            c = theta .* negAm.Li.cmax;
            c = c*ones(negAm.G.cells.num, 1);

            initstate.(ne).(eac).(am).c = c;
            initstate.(ne).(eac).(am) = negAm.updateMaterialProperties(initstate.(ne).(eac).(am));

            OCP = initstate.(ne).(eac).(am).OCP;
            initstate.(ne).(eac).(am).phi = OCP;

            %% setup initial PositiveElectrode state

            % shortcut
            % posAm : ActiveMaterial of the positive electrode
            
            posAm = bat.(pe).(eac).(am);
            
            m = (1 ./ (posAm.theta100 - posAm.theta0));
            b = -m .* posAm.theta0;
            theta = (SOC - b) ./ m;
            c = theta .* posAm.Li.cmax;
            c = c*ones(posAm.G.cells.num, 1);

            initstate.(pe).(eac).(am).c = c;
            initstate.(pe).(eac).(am) = posAm.updateMaterialProperties(initstate.(pe).(eac).(am));

            OCP = initstate.(pe).(eac).(am).OCP;
            initstate.(pe).(eac).(am).phi = OCP;

            %% setup initial Electrolyte state

            initstate.(elyte).phi = zeros(bat.(elyte).G.cells.num, 1);
            cs = cell(2,1);
            initstate.(elyte).cs = cs;
            initstate.(elyte).cs{1} = 1000*ones(bat.(elyte).G.cells.num, 1);

            %% setup initial Current collectors state

            OCP = initstate.(ne).(eac).(am).OCP;
            OCP = OCP(1) .* ones(bat.(ne).(cc).G.cells.num, 1);
            initstate.(ne).(cc).phi = OCP;

            OCP = initstate.(pe).(eac).(am).OCP;
            OCP = OCP(1) .* ones(bat.(pe).(cc).G.cells.num, 1);
            initstate.(pe).(cc).phi = OCP;
            
            initstate.(pe).(cc).E = OCP(1);
            
        end
        
        function state = setupElectrolyteCoupling(model, state)
        % Setup the electrolyte coupling by adding ion sources from the electrodes
        % shortcuts:
        % c_source : Source term for charge carrier.
                        
            bat = model;
            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            am    = 'ActiveMaterial';
            eac   = 'ElectrodeActiveComponent';
            
            ccSourceName = bat.(elyte).chargeCarrierSourceName;
            couplingterms = bat.couplingTerms;

            phi = state.(elyte).phi;

            elyte_c_source = zeros(bat.(elyte).G.cells.num, 1);
            
            if isa(phi, 'ADI')
                adsample = getSampleAD(phi);
                adbackend = model.AutoDiffBackend;
                elyte_c_source = adbackend.convertToAD(elyte_c_source, adsample);
            end
            
            coupnames = model.couplingNames;
            
            ne_R = state.(ne).(eac).(am).R;
            coupterm = getCoupTerm(couplingterms, 'NegativeElectrode-Electrolyte', coupnames);
            elytecells = coupterm.couplingcells(:, 2);
            elyte_c_source(elytecells) = ne_R;            
            
            pe_R = state.(pe).(eac).(am).R;
            coupterm = getCoupTerm(couplingterms, 'PositiveElectrode-Electrolyte', coupnames);
            elytecells = coupterm.couplingcells(:, 2);
            elyte_c_source(elytecells) = pe_R;
            
            state.Electrolyte.(ccSourceName) = elyte_c_source; 
        
        end
        
        function state = updateAccumTerms(model, state, state0, dt)
                    
            bat = model;
            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            am    = 'ActiveMaterial';
            eac   = 'ElectrodeActiveComponent';
            
            ccAccumName = bat.(elyte).chargeCarrierAccumName;
            
            cdotcc  = (state.(elyte).cs{1} - state0.(elyte).cs{1})/dt;
            ccAccum = bat.(elyte).volumeFraction.*cdotcc;
            state.(elyte).(ccAccumName) = ccAccum;
            
            names = {ne, pe};
            for i = 1 : numel(names)
                elde = names{i}; % electrode name
                cdotcc   = (state.(elde).(eac).(am).c - state0.(elde).(eac).(am).c)/dt;
                ccAccum  = bat.(elde).(eac).volumeFraction.*cdotcc;
                state.(elde).(eac).(ccAccumName) = ccAccum;
            end
            
        end

        
        function state = setupElectrodeCoupling(model, state)
        % Setup electrode coupling by updating the potential and concentration of the electrolyte in the active component of the
        % electrode. There, those quantities are considered as input and used to compute the reaction rate.
        %
        %
        % WARNING : at the moment, we do not pass the concentrations
        %
        % shortcuts:
        % elyte : Electrolyte
        % neac  : NegativeElectrode.ElectrodeActiveComponent 
        % peac  : PositiveElectrode.ElectrodeActiveComponent
            
            bat = model;
            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            am    = 'ActiveMaterial';
            eac   = 'ElectrodeActiveComponent';
            cc    = 'CurrentCollector';
            
            
            phi_elyte = state.(elyte).phi;
            
            elyte_cells = zeros(model.G.cells.num, 1);
            elyte_cells(bat.(elyte).G.mappings.cellmap) = (1 : bat.(elyte).G.cells.num)';

            phi_elyte_neac = phi_elyte(elyte_cells(bat.(ne).(eac).G.mappings.cellmap));
            phi_elyte_peac = phi_elyte(elyte_cells(bat.(ne).(eac).G.mappings.cellmap));

            state.(ne).(eac).(am).phiElectrolyte = phi_elyte_neac;
            state.(pe).(eac).(am).phiElectrolyte = phi_elyte_peac;
            
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
            [state.(elyte).cs{1}      , ...
             state.(elyte).phi        , ...   
             state.(ne).(eac).(am).c  , ...   
             state.(ne).(eac).(am).phi, ...   
             state.(pe).(eac).(am).c  , ...    
             state.(pe).(eac).(am).phi, ...   
             state.(ne).(cc).phi      , ...    
             state.(pe).(cc).phi      , ...    
             state.(pe).(cc).E] = ...
                adbackend.initVariablesAD(...
                    state.(elyte).cs{1}      , ...
                    state.(elyte).phi        , ...   
                    state.(ne).(eac).(am).c  , ...    
                    state.(ne).(eac).(am).phi, ...   
                    state.(pe).(eac).(am).c  , ...    
                    state.(pe).(eac).(am).phi, ...   
                    state.(ne).(cc).phi      , ...    
                    state.(pe).(cc).phi      , ...    
                    state.(pe).(cc).E);       
        end
        
        function p = getPrimaryVariables(model)
            bat = model;
            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            am    = 'ActiveMaterial';
            eac   = 'ElectrodeActiveComponent';
            cc    = 'CurrentCollector';
            
            p ={{elyte, 'cs', 1}       , ...
                {elyte, 'phi'}         , ...   
                {ne, eac, am, 'c'} , ...    
                {ne, eac, am, 'phi'}   , ...   
                {pe, eac, am, 'c'} , ...    
                {pe, eac, am, 'phi'}   , ...   
                {ne, cc, 'phi'}        , ...    
                {pe, cc, 'phi'}        , ...    
                {pe, cc, 'E'}
               };
            
        end
        
        function state = setProp(model, state, names, val)
            if iscell(names) & (numel(names) > 1)
                name = names{1};
                names = names(2 : end);
                state.(name) = model.setProp(state.(name), names, val);
            elseif iscell(names) & (numel(names) == 1)
                name = names{1};
                if isnumeric(name)
                    state{name} = val;
                else
                    state.(name) = val;
                end
            else
                error('format not recognized');                
            end                
        end
        
        function var = getProp(model, state, names)
            if iscell(names) && (numel(names) > 1)
                name = names{1};
                names = names(2 : end);
                var = model.getProp(state.(name), names);
            elseif iscell(names) & (numel(names) == 1)
                name = names{1};
                if isnumeric(name)
                    var = state{name};
                else
                    var = state.(name);
                end
            else
                error('format not recognized');
            end                
        end
        
        function submod = getSubmodel(model, names)
            submod = model.(names{1});
            for i=2:numel(names)
                submod = submod.(names{i});
            end
        end

        function validforces = getValidDrivingForces(model)
            validforces=struct('src', [], 'stopFunction', []); 
        end
         
        
        function [state, report] = updateState(model, state, problem, dx, drivingForces)
            p = model.getPrimaryVariables();
            for i=2:numel(dx)
                val = model.getProps(state, p{i});
                val = val + dx{i};
                state = model.setProp(state, p{i}, val);
            end
            %% not sure how to handle cells
            state.Electrolyte.cs{1} =  state.Electrolyte.cs{1} + dx{1};
            report = [];
        end
        
        
        function state = reduceState(model, state, removeContainers)
        % Reduce state to double (do not use property containers)
            state = value(state);
        end

        
    end
    
end
