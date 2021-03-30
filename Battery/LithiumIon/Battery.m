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
            
            model.AutoDiffBackend = SparseAutoDiffBackend('useBlocks',true);
            
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
            
            automaticAssembly;
            
            %% Set up the governing equations
            
            eqs={};
            
            %% We collect mass and charge conservation equations for the electrolyte and the electrodes

            names = {'Electrolyte', 'NegativeElectrode', 'PositiveElectrode'};
            
            for i = 1 : numel(names)
                eqs{end + 1} = model.getProp(state,{names{i}, 'massCons'});
                eqs{end + 1} = model.getProp(state,{names{i}, 'chargeCons'});
            end
            
            %% We collect charge conservation equations for the current collectors
            
            names = {'NegativeCurrentCollector', 'PositiveCurrentCollector'};
            for i = 1 : numel(names)
                eqs{end + 1} = model.getProp(state, {names{i}, 'chargeCons'});
            end
            
            %% We setup and add the control equation (fixed total current at PositiveCurrentCollector)
            
            src = drivingForces.src(time);
            coupterm = model.getCoupTerm('bc-PositiveCurrentCollector');
            faces = coupterm.couplingfaces;
            bcval = state.PositiveCurrentCollector.E;
            cond_pcc = model.PositiveCurrentCollector.EffectiveElectronicConductivity;
            [trans_pcc, cells] = model.PositiveCurrentCollector.operators.harmFaceBC(cond_pcc, faces);
            control = src - sum(trans_pcc.*(bcval - state.PositiveCurrentCollector.phi(cells)));
            
            eqs{end+1} = -control;


            
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

            initstate.(ne).(eac).(am).Li = c;
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

            initstate.(pe).(eac).(am).Li = c;
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
            ccName = bat.(elyte).chargeCarrierName;
            
            cdotcc  = (state.(elyte).cs{1} - state0.(elyte).cs{1})/dt;
            ccAccum = bat.(elyte).volumeFraction.*cdotcc;
            state.(elyte).(ccAccumName) = ccAccum;
            
            names = {ne, pe};
            for i = 1 : numel(names)
                elde = names{i}; % electrode name
                cdotcc   = (state.(elde).(eac).(am).(ccName) - state0.(elde).(eac).(am).(ccName))/dt;
                ccAccum  = bat.(elde).(eac).volumeFraction.*cdotcc;
                state.(elde).(ccAccumName) = ccAccum;
            end
            
        end

        
        function state = setupElectrodeCoupling(model, state)
        % Setup electrod coupling by updating the potential and concentration of the electrolyte in the active component of the
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
            [state.(elyte).cs{1}      ,...
             state.(elyte).phi        ,...   
             state.(ne).(eac).(am).Li ,...    
             state.(ne).(eac).(am).phi,...   
             state.(pe).(eac).(am).Li ,...    
             state.(pe).(eac).(am).phi,...   
             state.(ne).(cc).phi      ,...    
             state.(pe).(cc).phi      ,...    
             state.(pe).(cc).E] = ...
                adbackend.initVariablesAD(...
                    state.(elyte).cs{1}      ,...
                    state.(elyte).phi        ,...   
                    state.(ne).(eac).(am).Li ,...    
                    state.(ne).(eac).(am).phi,...   
                    state.(pe).(eac).(am).Li ,...    
                    state.(pe).(eac).(am).phi,...   
                    state.(ne).(cc).phi      ,...    
                    state.(pe).(cc).phi      ,...    
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
            
            p ={{elyte,'cs'}         , ...
                {elyte,'phi'}        , ...   
                {ne, eac, am, 'Li'}  , ...    
                {ne, eac, am, 'phi'} , ...   
                {pe, eac, am, 'Li'}  , ...    
                {pe, eac, am, 'phi'} , ...   
                {ne, cc,'phi'}       , ...    
                {pe, cc,'phi'}       , ...    
                {pe, cc,'E'}
               };
            
        end
        
        function state = setProp(model,state,names,val)
            nname=numel(names);
            if(nname==1)
                state.(names{1})=val;
            elseif(nname==2)
                state.(names{1}).(names{2})=val;
            elseif(nname==3)
                state.(names{1}).(names{2}).(names{3}) = val;
            else
                error('not implmented')
            end
        end
        
        function var = getProp(model, state, names)
            var = state.(names{1});
            for i = 2 : numel(names)
                var = var.(names{i});
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
         
        
        function [state, report] = updateState(model,state, problem, dx, drivingForces)
            p = model.getPrimaryVariables();
            for i=2:numel(dx)
                val = model.getProps(state,p{i});
                val = val + dx{i};
                state = model.setProp(state,p{i},val);
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
