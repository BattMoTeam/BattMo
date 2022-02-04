classdef Battery < BaseModel
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
        ThermalModel      % Thermal model, instance of :class:`ThermalComponent <Electrochemistry.ThermalComponent>`
        
        SOC % State Of Charge

        initT % Initial temperature
        
        Ucut % Voltage cut
        
        couplingTerms % Coupling terms
        cmin % mininum concentration used in capping

        couplingNames 
        
        mappings
        use_solid_diffusion
        use_thermal
    end
    
    methods
        
        function model = Battery(paramobj,varargin)
            opt = struct('use_solid_diffusion',true,'use_thermal',true);
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
            model.ThermalModel = ThermalComponent(paramobj.ThermalModel);
            
            % setup Electrolyte model (setup electrolyte volume fractions in the different regions)
            model = model.setupElectrolyteModel();            
            
            % setup Thermal Model by assigning the effective heat capacity and conductivity, which is computed from the sub-models.
            model = model.setupThermalModel();
            
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
            model.use_thermal = opt.use_thermal;
        end

        function model = setupThermalModel(model, paramobj)
        % Setup the thermal model :attr:`ThermalModel`. Here, :code:`paramobj` is instance of
        % :class:`ThermalComponentInputParams <Electrochemistry.ThermalComponentInputParams>`
            
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            eac   = 'ElectrodeActiveComponent';
            cc    = 'CurrentCollector';
            elyte = 'Electrolyte';
            sep   = 'Separator';
            
            eldes = {ne, pe}; % electrodes
            
            G = model.G;
            nc = G.cells.num;
            
            hcap = zeros(nc, 1); % effective heat capacity
            hcond = zeros(nc, 1); % effective heat conductivity
            
            for ind = 1 : numel(eldes)
                elde = eldes{ind};
                
                % The effecive and intrinsic thermal parameters for the current collector are the same.
                cc_map = model.(elde).(cc).G.mappings.cellmap;
                cc_hcond = model.(elde).(cc).thermalConductivity;
                cc_hcap = model.(elde).(cc).heatCapacity;

                hcap(cc_map) = hcap(cc_map) + cc_hcap;
                hcond(cc_map) = hcond(cc_map) + cc_hcond;
                
                % Effective parameters from the Electrode Active Component region.
                eac_map = model.(elde).(eac).G.mappings.cellmap;
                eac_hcond = model.(elde).(eac).thermalConductivity;
                eac_hcap = model.(elde).(eac).heatCapacity;
                eac_volfrac = model.(elde).(eac).volumeFraction;
                
                eac_hcap = eac_hcap.*eac_volfrac;
                eac_hcond = eac_hcond.*eac_volfrac.^1.5;
                
                hcap(eac_map) = hcap(eac_map) + eac_hcap;
                hcond(eac_map) = hcond(eac_map) + eac_hcond;
                
            end

            % Electrolyte
            
            elyte_map = model.(elyte).G.mappings.cellmap;
            elyte_hcond = model.(elyte).thermalConductivity;
            elyte_hcap = model.(elyte).heatCapacity;
            elyte_volfrac = model.(elyte).volumeFraction;
            
            elyte_hcap = elyte_hcap.*elyte_volfrac;
            elyte_hcond = elyte_hcond.*elyte_volfrac.^1.5;
            
            hcap(elyte_map) = hcap(elyte_map) + elyte_hcap;
            hcond(elyte_map) = hcond(elyte_map) + elyte_hcond;            
            
            % Separator
            
            sep_map = model.(elyte).(sep).G.mappings.cellmap;
            
            sep_hcond = model.(elyte).(sep).thermalConductivity;
            sep_hcap = model.(elyte).(sep).heatCapacity;
            sep_volfrac = model.(elyte).(sep).volumeFraction;
            
            sep_hcap = sep_hcap.*sep_volfrac;
            sep_hcond = sep_hcond.*sep_volfrac.^1.5;
            
            hcap(sep_map) = hcap(sep_map) + sep_hcap;
            hcond(sep_map) = hcond(sep_map) + sep_hcond;            

            model.ThermalModel.EffectiveHeatCapacity = hcap;
            model.ThermalModel.EffectiveThermalConductivity = hcond;
            
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
        
        function [SOCN,SOCP] = calculateSOC(model,state)
            bat = model;
            %elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            am    = 'ActiveMaterial';
            eac   = 'ElectrodeActiveComponent';
            %cc    = 'CurrentCollector';
            %thermal = 'ThermalModel';
            
            negAm = bat.(ne).(eac).(am);
            c =state.NegativeElectrode.ElectrodeActiveComponent.c;%ActiveMaterial.cElectrode
            theta = c/negAm.Li.cmax;
            negAm = bat.(ne).(eac).(am); 
            m = (1 ./ (negAm.theta100 - negAm.theta0));
            b = -m .* negAm.theta0;
            SOCN = theta*m+b;
            vol = model.NegativeElectrode.ElectrodeActiveComponent.volumeFraction;
            SOCN = sum(SOCN.*vol)/sum(vol);
            
            posAm = bat.(pe).(eac).(am);
            c =state.PositiveElectrode.ElectrodeActiveComponent.c;
            theta = c/posAm.Li.cmax;
            m = (1 ./ (posAm.theta100 - posAm.theta0));
            b = -m .* posAm.theta0;
            SOCP = theta*m+b;
            vol = model.PositiveElectrode.ElectrodeActiveComponent.volumeFraction;
            SOCP = sum(SOCP.*vol)/sum(vol);
        end
        
        function initstate = setupInitialState(model)
        % Setup the initial state

            nc = model.G.cells.num;

            SOC = model.SOC;
            T   = model.initT;
            
            %initstate.SOC = SOC*ones(nc, 1);
            
            bat = model;
            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            am    = 'ActiveMaterial';
            eac   = 'ElectrodeActiveComponent';
            cc    = 'CurrentCollector';
            thermal = 'ThermalModel';
            
            initstate.(thermal).T = T*ones(nc, 1);
            
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
            initstate.(elyte).c = 1000*ones(bat.(elyte).G.cells.num, 1);

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
            %state.SOC = model.SOC*ones(nc, 1);
            
            % Shortcuts used in this function
            battery = model;
            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            eac     = 'ElectrodeActiveComponent';
            cc      = 'CurrentCollector';
            elyte   = 'Electrolyte';
            am      = 'ActiveMaterial';
            thermal = 'ThermalModel';
            
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
                    state.(elde).(eac).c = state.(elde).(eac).(am).cElectrode;
                    state.(elde).(eac) = battery.(elde).(eac).updateChargeCarrier(state.(elde).(eac));                    
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
                state.(elde).(eac).(am) = battery.(elde).(eac).(am).updateDiffusionCoefficient(state.(elde).(eac).(am));
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
            state.(elyte) = battery.(elyte).updateMassFlux(state.(elyte));
            state.(elyte) = battery.(elyte).updateMassConservation(state.(elyte));

            
            for ind = 1 : numel(electrodes)
                elde = electrodes{ind};
                
                %% Electrodes mass conservation
                state.(elde).(eac) = battery.(elde).(eac).updateMassFlux(state.(elde).(eac));
                state.(elde).(eac) = battery.(elde).(eac).updateMassConservation(state.(elde).(eac));
                
                %% Electrodes charge conservation - current collector part
                state.(elde).(cc) = battery.(elde).(cc).updateCurrent(state.(elde).(cc));
                state.(elde).(cc) = battery.(elde).(cc).updateChargeConservation(state.(elde).(cc));

            end

            %% update solid diffustion equations
            for ind = 1 : numel(electrodes)
                elde = electrodes{ind};
                state.(elde).(eac).(am) = battery.(elde).(eac).(am).assembleSolidDiffusionEquation(state.(elde).(eac).(am));
            end
            
            %% update Thermal source term from electrical resistance

            state = battery.updateThermalOhmicSourceTerms(state);
            state = battery.updateThermalChemicalSourceTerms(state);
            state = battery.updateThermalReactionSourceTerms(state);
            
            state.(thermal) = battery.(thermal).updateHeatSourceTerm(state.(thermal));
            state.(thermal) = battery.(thermal).updateThermalBoundarySourceTerms(state.(thermal));
            
            %% update Accumulation terms for the energy equation
            
            state = battery.updateThermalAccumTerms(state, state0, dt);
            
            %% Update energy conservation residual term
            
            state.(thermal) = model.(thermal).updateEnergyConservation(state.(thermal));
            
            %% setup relation between E and I at positive current collectror
            
            state = model.setupEIEquation(state);
            
            %% Set up the governing equations
            
            eqs = {};
            
            %% We collect mass and charge conservation equations for the electrolyte and the electrodes
            massConsScaling = model.con.F;
            
            eqs{end + 1} = state.(elyte).massCons*massConsScaling;
            eqs{end + 1} = state.(elyte).chargeCons;
            
            eqs{end + 1} = state.(ne).(eac).massCons*massConsScaling;
            eqs{end + 1} = state.(ne).(eac).chargeCons;
            eqs{end + 1} = state.(ne).(eac).(am).solidDiffusionEq.*massConsScaling.*battery.(ne).(eac).(am).G.cells.volumes/dt;
            
            eqs{end + 1} = state.(pe).(eac).massCons*massConsScaling;
            eqs{end + 1} = state.(pe).(eac).chargeCons;
            eqs{end + 1} = state.(pe).(eac).(am).solidDiffusionEq.*massConsScaling.*battery.(pe).(eac).(am).G.cells.volumes/dt;
            
            eqs{end + 1} = state.(ne).(cc).chargeCons;
            eqs{end + 1} = state.(pe).(cc).chargeCons;
            
            eqs{end + 1} = state.(thermal).energyCons;
            
            eqs{end + 1} = -state.EIeq;
            
            % we add the control equation
            I = state.(pe).(cc).I;
            E = state.(pe).(cc).E;
            [val, ctrltype] = drivingForces.src(time, value(I), value(E));
            switch ctrltype
              case 'I'
                eqs{end + 1} = I - val;
              case 'E'
                eqs{end + 1} = (E - val)*1e5;
            end

            eqs{1} = eqs{1} - model.Electrolyte.sp.t(1)*eqs{2};
            
            %% Give type and names to equations and names of the primary variables (for book-keeping)
            
 
            
            types = {'cell','cell','cell','cell', 'sdiff','cell','cell','cdiff','cell','cell', 'cell', 'cntrl', 'cntrl'};
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
                     'energyCons'       , ...
                     'EIeq', ...
                     'controlEq'};
            neq=numel(eqs);     
            keep = model.getEquationsToUses(neq);     
 
            if(not(all(keep)))
                ind = find(keep);
                eqs={eqs{ind}};
                types = {types{ind}};
                names = {names{ind}};
            end
            switch ctrltype
                case 'I'
                   types{end-1} = 'cell';   
                case 'E'
                   neqs = numel(types);
                   order = [1:neqs-2,neqs,neqs-1];
                   types = { types{order} };
                   eqs = {eqs{order}};
                   names = {names{order}};
                otherwise 
                          error()
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
            thermal = 'ThermalModel';
            
            % (here we assume that the ThermalModel has the "parent" grid)
            state.(elyte).T    = state.(thermal).T(model.(elyte).G.mappings.cellmap);
            state.(ne).(eac).T = state.(thermal).T(model.(ne).(eac).G.mappings.cellmap);
            state.(ne).(cc).T  = state.(thermal).T(model.(ne).(cc).G.mappings.cellmap);
            state.(pe).(eac).T = state.(thermal).T(model.(pe).(eac).G.mappings.cellmap);
            state.(pe).(cc).T  = state.(thermal).T(model.(pe).(cc).G.mappings.cellmap);
            
            % Update temperature in the active materials of the electrodes.
            state.(ne).(eac) = model.(ne).(eac).updateTemperature(state.(ne).(eac));
            state.(pe).(eac) = model.(pe).(eac).updateTemperature(state.(pe).(eac));
            
        end
        
        function keep = getEquationsToUses(model,neq)
             keep = true(neq,1);     
            if(not(model.use_solid_diffusion))
                keep(5)=false;
                keep(8)=false;
            end
            if(not(model.use_thermal))
                keep(11)=false;
            end
        end
        
        function keep = getVariablesToUses(model,neq)
             %% reduction will not work if not this is equalt to equations: if need one need to change LinearSolverAD.m
             keep = model.getEquationsToUses(neq);
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
            
            state.Electrolyte.massSource = elyte_c_source; 
            state.Electrolyte.eSource = elyte_e_source;
            
        end
        
        function state = updateAccumTerms(model, state, state0, dt)
        % Assemble the accumulation terms for transport equations (in electrolyte and electrodes)
            
            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            am    = 'ActiveMaterial';
            eac   = 'ElectrodeActiveComponent';
            
            cdotcc  = (state.(elyte).c - state0.(elyte).c)/dt;
            effectiveVolumes = model.(elyte).volumeFraction.*model.(elyte).G.cells.volumes;
            massAccum  = effectiveVolumes.*cdotcc;
            state.(elyte).massAccum = massAccum;
            
            names = {ne, pe};
            for i = 1 : numel(names)
                elde = names{i}; % electrode name
                cdotcc   = (state.(elde).(eac).c - state0.(elde).(eac).c)/dt;
                effectiveVolumes = model.(elde).(eac).volumeFraction.*model.(elde).(eac).G.cells.volumes;
                massAccum  = effectiveVolumes.*cdotcc;
                state.(elde).(eac).massAccum = massAccum;
            end
            
        end

        
        function state = updateThermalAccumTerms(model, state, state0, dt)
        % Assemble the accumulation term for the energy equation
            thermal = 'ThermalModel';
            
            hcap = model.(thermal).EffectiveHeatCapacity;
            
            T = state.(thermal).T;
            T0 = state0.(thermal).T;

            % (here we assume that the ThermalModel has the "parent" grid)
            vols = model.G.cells.volumes;
            
            state.(thermal).accumHeat = hcap.*vols.*(T - T0)/dt;
            
        end


        function state = updateThermalOhmicSourceTerms(model, state)
        % Assemble the ohmic source term :code:`state.jHeatOhmSource`, see :cite:t:`Latz2016`

            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            eac     = 'ElectrodeActiveComponent';
            cc      = 'CurrentCollector';
            elyte   = 'Electrolyte';
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
                
                cc_model = model.(elde).(cc);
                cc_map   = cc_model.G.mappings.cellmap;
                cc_j     = locstate.(elde).(cc).j;
                cc_econd = cc_model.EffectiveElectricalConductivity;
                cc_vols  = cc_model.G.cells.volumes;
                cc_jsq   = computeCellFluxNorm(cc_model, cc_j); 
                state.(elde).(cc).jsq = cc_jsq;  %store square of current density
                
                src(cc_map) = src(cc_map) + cc_vols.*cc_jsq./cc_econd;

                eac_model = model.(elde).(eac);
                eac_map   = eac_model.G.mappings.cellmap;
                eac_j     = locstate.(elde).(eac).j;
                eac_econd = eac_model.EffectiveElectricalConductivity;
                eac_vols   = eac_model.G.cells.volumes;
                eac_jsq   = computeCellFluxNorm(eac_model, eac_j);
                state.(elde).(eac).jsq = eac_jsq;
                
                src(eac_map) = src(eac_map) + eac_vols.*eac_jsq./eac_econd;
                
            end

            % Electrolyte
            elyte_model = model.(elyte);
            elyte_map   = elyte_model.G.mappings.cellmap;
            elyte_vf    = elyte_model.volumeFraction;
            elyte_j     = locstate.(elyte).j;
            elyte_cond  = locstate.(elyte).conductivity;
            elyte_econd = elyte_cond.*elyte_vf.^1.5;
            elyte_vols  = elyte_model.G.cells.volumes;
            elyte_jsq   = computeCellFluxNorm(elyte_model, elyte_j);
            state.(elyte).jsq = elyte_jsq; %store square of current density
            
            src(elyte_map) = src(elyte_map) + elyte_vols.*elyte_jsq./elyte_econd;
            
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
            if isa(T, 'ADI')
                adsample = getSampleAD(T);
                adbackend = model.AutoDiffBackend;
                src = adbackend.convertToAD(src, adsample);
                locstate = state;
            else
                locstate = value(state);
            end

            % Compute chemical heat source in electrolyte
            dmudcs = locstate.(elyte).dmudcs;   % Derivative of chemical potential with respect to concentration
            D      = locstate.(elyte).D;        % Effective diffusion coefficient 
            Ddc    = locstate.(elyte).diffFlux; % Diffusion flux (-D*grad(c))
            
            % compute norm of square norm of diffusion flux
            elyte_model = model.(elyte);
            elyte_map   = elyte_model.G.mappings.cellmap;
            elyte_vols  = elyte_model.G.cells.volumes;
            elyte_jchemsq = computeCellFluxNorm(elyte_model, Ddc);
            elyte_src = elyte_vols.*elyte_jchemsq./D;
            
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
            eac     = 'ElectrodeActiveComponent';
            am      = 'ActiveMaterial';
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
                am_model = model.(elde).(eac).(am);
                am_map   = am_model.G.mappings.cellmap;
                R = locstate.(elde).(eac).(am).R;
                eta = locstate.(elde).(eac).(am).eta;
                vols = model.(elde).(eac).G.cells.volumes;
                am_src = vols.*R.*eta;
                
                src(am_map) = src(am_map) + am_src;
                
            end

            % We multiply by volumes
            src = model.G.cells.volumes.*src;
            
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
       function state = initStateAD(model, state)
            [pnames,extras]  = model.getPrimaryVariables();
            vars = cell(numel(pnames),1);
            for i=1:numel(pnames)
                vars{i} = model.getProp(state,pnames{i});
            end
            % Get the AD state for this model           
            [vars{:}] = model.AutoDiffBackend.initVariablesAD(vars{:});
            newstate =struct();
            for i=1:numel(pnames)
               newstate = model.setNewProp(newstate,pnames{i},vars{i});
            end
            
            for i=1:numel(extras)
                var = model.getProp(state,extras{i});
                assert(isnumeric(var));
                newstate = model.setNewProp(newstate,extras{i},var);
            end
            time = state.time;
            state = newstate;
            state.time = time;
        end 
%         
%         function state = initStateAD(model,state)
%             
%             bat = model;
%             elyte = 'Electrolyte';
%             ne    = 'NegativeElectrode';
%             pe    = 'PositiveElectrode';
%             am    = 'ActiveMaterial';
%             eac   = 'ElectrodeActiveComponent';
%             cc    = 'CurrentCollector';
%             thermal = 'ThermalModel';
%             
%             adbackend = model.AutoDiffBackend();
%             useMex=false;
%             if(isprop(adbackend,'useMex'))
%                 useMex = adbackend.useMex; 
%             end
%             if(model.use_solid_diffusion)
%             opts=struct('types',[1,1,2,2,2,3,3,3,4,5,6,7,8],'useMex',useMex);
%             [state.(elyte).c  , ...
%              state.(elyte).phi    , ...   
%              state.(ne).(eac).c   , ...   
%              state.(ne).(eac).phi , ...   
%              state.(ne).(eac).(am).cElectrode , ...
%              state.(pe).(eac).c   , ...    
%              state.(pe).(eac).phi , ...   
%              state.(pe).(eac).(am).cElectrode , ...
%              state.(ne).(cc).phi  , ...    
%              state.(pe).(cc).phi  , ...    
%              state.(thermal).T    , ...
%              state.(pe).(cc).E, ...
%              state.(pe).(cc).I] = ...
%                 adbackend.initVariablesAD(...
%                     state.(elyte).c  , ...
%                     state.(elyte).phi    , ...   
%                     state.(ne).(eac).c   , ...    
%                     state.(ne).(eac).phi , ...   
%                     state.(ne).(eac).(am).cElectrode , ...
%                     state.(pe).(eac).c   , ...    
%                     state.(pe).(eac).phi , ...   
%                     state.(pe).(eac).(am).cElectrode , ...
%                     state.(ne).(cc).phi  , ...    
%                     state.(pe).(cc).phi  , ...    
%                     state.(thermal).T    , ...
%                     state.(pe).(cc).E, ...
%                     state.(pe).(cc).I, ...
%                     opts);
%             else
%             opts=struct('types',[1,1,2,2,3,3,4,5,6,7,8],'useMex',useMex);
%             [state.(elyte).c  , ...
%              state.(elyte).phi    , ...   
%              state.(ne).(eac).c   , ...   
%              state.(ne).(eac).phi , ...   
%              state.(pe).(eac).c   , ...    
%              state.(pe).(eac).phi , ...   
%              state.(ne).(cc).phi  , ...    
%              state.(pe).(cc).phi  , ...    
%              state.(thermal).T    , ...
%              state.(pe).(cc).E, ...
%              state.(pe).(cc).I] = ...
%                 adbackend.initVariablesAD(...
%                     state.(elyte).c  , ...
%                     state.(elyte).phi    , ...   
%                     state.(ne).(eac).c   , ...    
%                     state.(ne).(eac).phi , ...   
%                     state.(pe).(eac).c   , ...    
%                     state.(pe).(eac).phi , ...   
%                     state.(ne).(cc).phi  , ...    
%                     state.(pe).(cc).phi  , ...    
%                     state.(thermal).T    , ...
%                     state.(pe).(cc).E, ...
%                     state.(pe).(cc).I, ...
%                     opts);
%                 
%             end
%            % PRIMARY variables
%        end
        
        function [p,extra] = getPrimaryVariables(model)
            
            bat = model;
            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            am    = 'ActiveMaterial';
            eac   = 'ElectrodeActiveComponent';
            cc    = 'CurrentCollector';
            thermal = 'ThermalModel';
            %% chaged order of c cElectrode to get easier reduction
              p = {{elyte, 'c'}                , ...
                 {elyte, 'phi'}              , ...   
                 {ne, eac, am, 'cElectrode'}              , ...    
                 {ne, eac, 'phi'}            , ...   
                 {ne, eac, 'c'} , ...
                 {pe, eac, am, 'cElectrode'}              , ...    
                 {pe, eac, 'phi'}            , ...   
                 {pe, eac, 'c'} , ...
                 {ne, cc, 'phi'}             , ...    
                 {pe, cc, 'phi'}             , ...
                 {thermal, 'T'}              , ...
                 {pe, cc, 'E'}               , ...
                 {pe, cc, 'I'}};          
           
%             p = {{elyte, 'c'}                , ...
%                  {elyte, 'phi'}              , ...   
%                  {ne, eac, 'c'}              , ...    
%                  {ne, eac, 'phi'}            , ...   
%                  {ne, eac, am, 'cElectrode'} , ...
%                  {pe, eac, 'c'}              , ...    
%                  {pe, eac, 'phi'}            , ...   
%                  {pe, eac, am, 'cElectrode'} , ...
%                  {ne, cc, 'phi'}             , ...    
%                  {pe, cc, 'phi'}             , ...
%                  {thermal, 'T'}              , ...
%                  {pe, cc, 'E'}               , ...
%                  {pe, cc, 'I'}};
            neq = numel(p);
            keep = model.getVariablesToUses(neq);
            extra = {p{find(not(keep))}};
           if(not(all(keep)))
               p = {p{find(keep)}};
           end
            
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
            %for i=1:numel(mnames)
            %    mname=mnames{i}
            %    submodel=model.getSubmodel(mname);
            %    submodel.AutoDiffBackend = model.AutoDiffBackend;
            %    submodel=submodel.validateModel(varargin{:});
            %    model  = model.setProp(model,mname,submodel);
            %end
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
            
            state.(elyte).c = max(cmin, state.(elyte).c);
            
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
