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
        Control           % Control Model
        
        SOC % State Of Charge

        initT % Initial temperature
        
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
            fdnames = {'G'            , ...
                       'couplingTerms', ...
                       'initT'        , ...
                       'SOC'          , ...
                       'I'};
            
            model = dispatchParams(model, paramobj, fdnames);

            % Assign the components : Electrolyte, NegativeElectrode, PositiveElectrode
            model.NegativeElectrode = model.setupElectrode(paramobj.NegativeElectrode);
            model.PositiveElectrode = model.setupElectrode(paramobj.PositiveElectrode);
            model.Electrolyte       = model.setupElectrolyte(paramobj.Electrolyte);
            model.ThermalModel      = ThermalComponent(paramobj.ThermalModel);
            model.Control           = model.setupControl(paramobj.Control);
            
            % define shorthands
            elyte   = 'Electrolyte';
            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            am      = 'ActiveMaterial';
            cc      = 'CurrentCollector';
            am      = 'ActiveMaterial';
            itf     = 'Interface';
            thermal = 'ThermalModel';
           
            % setup Electrolyte model (setup electrolyte volume fractions in the different regions)
            model = model.setupElectrolyteModel();            
            
            % setup Thermal Model by assigning the effective heat capacity and conductivity, which is computed from the sub-models.
            model = model.setupThermalModel();
            
            % setup couplingNames
            model.couplingNames = cellfun(@(x) x.name, model.couplingTerms, 'uniformoutput', false);
            
            % setup some mappings (mappings from electrodes to electrolyte)
            model = model.setupMappings();
            
            % setup capping
            cmax_ne = model.(ne).(am).(itf).cmax;
            cmax_pe = model.(pe).(am).(itf).cmax;
            model.cmin = 1e-5*max(cmax_ne, cmax_pe);
            model.use_solid_diffusion = opt.use_solid_diffusion;
            model.use_thermal = opt.use_thermal;
            
        end

        
        function model = registerVarAndPropfuncNames(model)
            
            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)
            
            model = registerVarAndPropfuncNames@BaseModel(model);
            
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
            
            varnames = {'SOC', ...
                        'controlEq'};
            
            model = model.registerVarNames(varnames);
            
            %% Temperature dispatch functions
            fn = @Battery.updateTemperature;
            
            inputnames = {{thermal, 'T'}};
            model = model.registerPropFunction({{ne, am, 'T'} , fn, inputnames});
            model = model.registerPropFunction({{ne, cc , 'T'}, fn, inputnames});
            model = model.registerPropFunction({{pe, am, 'T'} , fn, inputnames});
            model = model.registerPropFunction({{pe, cc , 'T'}, fn, inputnames});  
            model = model.registerPropFunction({{elyte, 'T'}  , fn, inputnames});
                  
            %% Coupling functions
            
            % dispatch electrolyte concentration and potential in the electrodes
            fn = @Battery.updateElectrodeCoupling;
            inputnames = {{elyte, 'c'}, ...
                          {elyte, 'phi'}};
            model = model.registerPropFunction({{ne, am, itf, 'phiElectrolyte'}, fn, inputnames});
            model = model.registerPropFunction({{ne, am, itf, 'cElectrolyte'}  , fn, inputnames});
            model = model.registerPropFunction({{pe, am, itf, 'phiElectrolyte'}, fn, inputnames});
            model = model.registerPropFunction({{pe, am, itf, 'cElectrolyte'}  , fn, inputnames});
            
            % Functions that update the source terms in the electolyte
            
            fn = @Battery.updateElectrolyteCoupling;
            
            inputnames = {{ne, am, itf, 'R'}, ...
                          {pe, am, itf, 'R'}};
            model = model.registerPropFunction({{elyte, 'cSource'}, fn, inputnames});
            model = model.registerPropFunction({{elyte, 'eSource'}, fn, inputnames});
            
            % Function that assemble the control equation
            
            fn = @Batter.setupEIEquation;
            inputnames = {{'pe', 'cc', 'E'}, ...
                          {'pe', 'cc', 'I'}, ...
                          {'pe', 'cc', 'phi'}, ...
                         };
            model = model.registerPropFunction({{'controlEq'}, fn, inputnames});
            
            % Function that update Thermal accumulation terms
            
            fn = @Battery.updateThermalAccumTerms;
            inputnames = {'T'};
            model = model.registerPropFunction({{thermal, 'accumHeat'}, fn, inputnames});
            
            % Function that update the Thermal Ohmic Terms
            
            fn = @Battery.updateThermalOhmicSourceTerms;
            inputnames = {{elyte, 'j'}   , ...
                          {ne, cc, 'j'}  , ...
                          {ne, am, 'j'} , ...
                          {pe, cc, 'j'}  , ...
                          {pe, am, 'j'}};
            model = model.registerPropFunction({{thermal, 'jHeatOhmSource'}, fn, inputnames});
            model = model.registerPropFunction({{thermal, 'jHeatBcSource'} , fn, inputnames});
            
            %% Function that updates the Thermal Chemical Terms
            
            fn = @Battery.updateThermalChemicalSourceTerms;
            inputnames = {{elyte, 'diffFlux'}, ...
                          {elyte, 'D'}       , ...
                          {elyte, 'dmudcs'}};
            model = model.registerPropFunction({{thermal, 'jHeatChemicalSource'}, fn, inputnames});
                          
            %% Functio that updates Thermal Reaction Terms
            
            fn = @Battery.updateThermalReactionSourceTerms;
            inputnames = {{ne, am, itf, 'R'}  , ...
                          {ne, am, itf, 'eta'}, ...
                          {pe, am, itf, 'R'}  , ...
                          {pe, am, itf, 'eta'}};
            model = model.registerPropFunction({{thermal, 'jHeatReactionSource'}, fn, inputnames});
                                                    
            %% Function that setup external coupling at positive and negative electrodes
            
            fn = @Battery.setupExternalCouplingNegativeElectrode;
            model = model.registerPropFunction({{ne, cc, 'jExternal'}, fn, {'phi'}});
            
            fn = @Battery.setupExternalCouplingPositiveElectrode;
            model = model.registerPropFunction({{pe, cc, 'jExternal'}, fn, {'phi', 'E'}});
        end

        function control = setupControl(model, paramobj)


            switch paramobj.controlPolicy
              case "IEswitch"
                control = ControlModel(paramobj); 
              case "CCCV"
                control = CcCvControlModel(paramobj);
              otherwise
                error('Error controlPolicy not recognized');
            end
            
            C = computeCellCapacity(model);
            CRate = control.CRate;
            
            control.Imax = (C/hour)*CRate;
            
        end
        
        function model = setupThermalModel(model, paramobj)
        % Setup the thermal model :attr:`ThermalModel`. Here, :code:`paramobj` is instance of
        % :class:`ThermalComponentInputParams <Electrochemistry.ThermalComponentInputParams>`
            
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            am    = 'ActiveMaterial';
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
                am_map = model.(elde).(am).G.mappings.cellmap;
                am_hcond = model.(elde).(am).thermalConductivity;
                am_hcap = model.(elde).(am).heatCapacity;
                am_volfrac = model.(elde).(am).volumeFraction;
                
                am_hcap = am_hcap.*am_volfrac;
                am_hcond = am_hcond.*am_volfrac.^1.5;
                
                hcap(am_map) = hcap(am_map) + am_hcap;
                hcond(am_map) = hcond(am_map) + am_hcond;
                
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
            electrolyte = Electrolyte(paramobj);
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
        
        function model = setupElectrolyteModel(model)
        % Assign the electrolyte volume fractions in the different regions

            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            am    = 'ActiveMaterial';
            sep   = 'Separator';

            elyte_cells = zeros(model.G.cells.num, 1);
            elyte_cells(model.(elyte).G.mappings.cellmap) = (1 : model.(elyte).G.cells.num)';

            model.(elyte).volumeFraction = NaN(model.(elyte).G.cells.num, 1);
            model.(elyte).volumeFraction(elyte_cells(model.(ne).(am).G.mappings.cellmap)) = model.(ne).(am).porosity;
            model.(elyte).volumeFraction(elyte_cells(model.(pe).(am).G.mappings.cellmap)) = model.(pe).(am).porosity;
            model.(elyte).volumeFraction(elyte_cells(model.(elyte).(sep).G.mappings.cellmap)) = model.(elyte).(sep).porosity;

        end
        
        function [SOCN, SOCP] = calculateSOC(model, state)
            
            bat = model;

            ne  = 'NegativeElectrode';
            pe  = 'PositiveElectrode';
            am  = 'ActiveMaterial';
            itf = 'Interface';
            
            negAm = bat.(ne).(am).(itf);
            c     = state.(ne).(am).c;
            theta = c/negAm.cmax;
            m     = (1 ./ (negAm.theta100 - negAm.theta0));
            b     = -m .* negAm.theta0;
            SOCN  = theta*m + b;
            vol   = model.(ne).(am).volumeFraction;
            
            SOCN = sum(SOCN.*vol)/sum(vol);
            
            posAm = bat.(pe).(am).(itf);
            c     = state.(pe).(am).c;
            theta = c/posAm.cmax;
            m     = (1 ./ (posAm.theta100 - posAm.theta0));
            b     = -m .* posAm.theta0;
            SOCP  = theta*m + b;
            vol   = model.(pe).(am).volumeFraction;
            
            SOCP = sum(SOCP.*vol)/sum(vol);
            
        end
        
        function initstate = setupInitialState(model)
        % Setup the initial state

            nc = model.G.cells.num;

            SOC = model.SOC;
            T   = model.initT;
            
            %initstate.SOC = SOC*ones(nc, 1);
            
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
            
            %% synchronize temperatures
            
            initstate = model.updateTemperature(initstate);
            
            %% setup initial NegativeElectrode state
            
            % shortcut
            % negItf : ActiveMaterial of the negative electrode
            
            negItf = bat.(ne).(am).(itf); 
            
            m     = (1 ./ (negItf.theta100 - negItf.theta0));
            b     = -m .* negItf.theta0;
            theta = (SOC - b) ./ m;
            c     = theta .* negItf.cmax;
            c     = c*ones(negItf.G.cells.num, 1);

            initstate.(ne).(am).c = c;
            initstate.(ne).(am).(sd).cSurface = c;
            % We bypass the solid diffusion equation to set directly the particle surface concentration (this is a bit hacky)
            initstate.(ne).(am) = model.(ne).(am).updateConcentrations(initstate.(ne).(am));
            initstate.(ne).(am).(itf) = negItf.updateOCP(initstate.(ne).(am).(itf));

            OCP = initstate.(ne).(am).(itf).OCP;
            ref = OCP(1);
            
            initstate.(ne).(am).phi = OCP - ref;

            %% setup initial PositiveElectrode state

            % shortcut
            % posItf : ActiveMaterial of the positive electrode
            
            posItf = bat.(pe).(am).(itf);
            
            m     = (1 ./ (posItf.theta100 - posItf.theta0));
            b     = -m .* posItf.theta0;
            theta = (SOC - b) ./ m;
            c     = theta .* posItf.cmax;
            c     = c*ones(posItf.G.cells.num, 1);

            initstate.(pe).(am).c = c;
            initstate.(pe).(am).(sd).cSurface = c;
            % We bypass the solid diffusion equation to set directly the particle surface concentration (this is a bit hacky)
            initstate.(pe).(am) = model.(pe).(am).updateConcentrations(initstate.(pe).(am));
            initstate.(pe).(am).(itf) = posItf.updateOCP(initstate.(pe).(am).(itf));
            
            OCP = initstate.(pe).(am).(itf).OCP;
            initstate.(pe).(am).phi = OCP - ref;

            %% setup initial Electrolyte state

            initstate.(elyte).phi = zeros(bat.(elyte).G.cells.num, 1)-ref;
            initstate.(elyte).c = 1000*ones(bat.(elyte).G.cells.num, 1);

            %% setup initial Current collectors state

            OCP = initstate.(ne).(am).(itf).OCP;
            OCP = OCP(1) .* ones(bat.(ne).(cc).G.cells.num, 1);
            initstate.(ne).(cc).phi = OCP - ref;

            OCP = initstate.(pe).(am).(itf).OCP;
            OCP = OCP(1) .* ones(bat.(pe).(cc).G.cells.num, 1);
            initstate.(pe).(cc).phi = OCP - ref;
            
            initstate.(ctrl).E = OCP(1) - ref;
            initstate.(ctrl).I = - model.(ctrl).Imax;
            
            switch model.(ctrl).controlPolicy
              case 'CCCV'
                initstate.(ctrl).ctrlType = 'CC_charge';
                initstate.(ctrl).nextCtrlType = 'CC_charge';
              case 'IEswitch'
                initstate.(ctrl).ctrlType = 'constantCurrent';
              otherwise
                error('control policy not recognized');
            end
            
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
            
            % Shorthands used in this function
            battery = model;
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
            
            electrodes = {ne, pe};
            electrodecomponents = {am, cc};

            %% Synchronization across components

            % temperature
            state = battery.updateTemperature(state);

            state.(elyte) = battery.(elyte).updateConcentrations(state.(elyte));
            
            for ind = 1 : numel(electrodes)
                elde = electrodes{ind};
                % potential and concentration between interface and active material
                state.(elde).(am) = battery.(elde).(am).updatePhi(state.(elde).(am));
                if (model.use_solid_diffusion)
                    state.(elde).(am) = battery.(elde).(am).updateConcentrations(state.(elde).(am));
                else
                    state.(elde).(am).c = state.(elde).(am).(itf).cElectrode;
                    state.(elde).(am) = battery.(elde).(am).updateConcentrations(state.(elde).(am));
                end              
            end
            
            %% Accumulation terms

            state = battery.updateAccumTerms(state, state0, dt);

            %% Update Electrolyte -> Electrodes coupling 
            
            state = battery.updateElectrodeCoupling(state); 

            %% Update reaction rates in both electrodes

            for ind = 1 : numel(electrodes)
                elde = electrodes{ind};
                state.(elde).(am).(itf) = battery.(elde).(am).(itf).updateReactionRateCoefficient(state.(elde).(am).(itf));
                state.(elde).(am).(itf) = battery.(elde).(am).(itf).updateOCP(state.(elde).(am).(itf));
                state.(elde).(am).(itf) = battery.(elde).(am).(itf).updateReactionRate(state.(elde).(am).(itf));
            end

            %% Update Electrodes -> Electrolyte  coupling

            state = battery.updateElectrolyteCoupling(state);
            
            %% Update coupling within electrodes and external coupling
            
            state.(ne) = battery.(ne).updateCoupling(state.(ne));
            state.(pe) = battery.(pe).updateCoupling(state.(pe));

            state.(ne).(am) = battery.(ne).(am).updatejBcSource(state.(ne).(am));
            state.(pe).(am) = battery.(pe).(am).updatejBcSource(state.(pe).(am));
            
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
                state.(elde).(am) = battery.(elde).(am).updateCurrentSource(state.(elde).(am));
                state.(elde).(am) = battery.(elde).(am).updateCurrent(state.(elde).(am));
                state.(elde).(am) = battery.(elde).(am).updateChargeConservation(state.(elde).(am));
            end
            
            %% elyte mass conservation

            state.(elyte) = battery.(elyte).updateDiffusionCoefficient(state.(elyte));
            state.(elyte) = battery.(elyte).updateMassFlux(state.(elyte));
            state.(elyte) = battery.(elyte).updateMassConservation(state.(elyte));

            
            for ind = 1 : numel(electrodes)
                elde = electrodes{ind};
                
                %% Electrodes mass conservation
                state.(elde).(am) = battery.(elde).(am).updateMassFlux(state.(elde).(am));
                state.(elde).(am) = battery.(elde).(am).updateMassSource(state.(elde).(am));
                state.(elde).(am) = battery.(elde).(am).updateMassConservation(state.(elde).(am));
                
                %% Electrodes charge conservation - current collector part
                state.(elde).(cc) = battery.(elde).(cc).updateCurrent(state.(elde).(cc));
                state.(elde).(cc) = battery.(elde).(cc).updateChargeConservation(state.(elde).(cc));

            end

            %% update solid diffustion equations
            for ind = 1 : numel(electrodes)
                elde = electrodes{ind};
                state.(elde).(am).(sd) = battery.(elde).(am).(sd).updateDiffusionCoefficient(state.(elde).(am).(sd));
                state.(elde).(am) = battery.(elde).(am).dispatchRate(state.(elde).(am));
                state.(elde).(am).(sd) = battery.(elde).(am).(sd).assembleSolidDiffusionEquation(state.(elde).(am).(sd));
            end
            
            %% update Face fluxes
            for ind = 1 : numel(electrodes)
                elde = electrodes{ind};
                state.(elde).(am) = battery.(elde).(am).updateFaceCurrent(state.(elde).(am));
                state.(elde).(cc) = battery.(elde).(cc).updateFaceCurrent(state.(elde).(cc));
            end
            state.(elyte) = battery.(elyte).updateFaceCurrent(state.(elyte));
            
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
            state = model.updateControl(state, drivingForces);
            state.(ctrl) = model.(ctrl).updateControlEquation(state.(ctrl));
            
            %% Set up the governing equations
            
            eqs = {};
            
            %% We collect the governing equations
            % The governing equations are the mass and charge conservation equations for the electrolyte and the
            % electrodes and the solid diffusion model equations and the control equations. The equations are scaled to
            % a common value.

            massConsScaling = model.con.F;
            
            eqs{end + 1} = state.(elyte).massCons*massConsScaling;
            eqs{end + 1} = state.(elyte).chargeCons;
            
            eqs{end + 1} = state.(ne).(am).massCons*massConsScaling;
            eqs{end + 1} = state.(ne).(am).chargeCons;
            eqs{end + 1} = state.(ne).(am).(sd).solidDiffusionEq.*massConsScaling.*battery.(ne).(am).(itf).G.cells.volumes/dt;
            
            eqs{end + 1} = state.(pe).(am).massCons*massConsScaling;
            eqs{end + 1} = state.(pe).(am).chargeCons;
            eqs{end + 1} = state.(pe).(am).(sd).solidDiffusionEq.*massConsScaling.*battery.(pe).(am).(itf).G.cells.volumes/dt;
            
            eqs{end + 1} = state.(ne).(cc).chargeCons;
            eqs{end + 1} = state.(pe).(cc).chargeCons;
            
            eqs{end + 1} = state.(thermal).energyCons;
            
            eqs{end + 1} = -state.(ctrl).EIequation;
            eqs{end + 1} = state.(ctrl).controlEquation;

            eqs{1} = eqs{1} - model.Electrolyte.sp.t(1)*eqs{2};
            
            %% Give type and names to equations and names of the primary variables (for book-keeping)
            
            types = {'cell','cell','cell','cell', 'sdiff','cell','cell','cdiff','cell','cell', 'cell', 'cntrl', 'cntrl'};
            names = {'elyte_massCons'   , ...
                     'elyte_chargeCons' , ...
                     'ne_am_massCons'  , ...
                     'ne_am_chargeCons', ...
                     'ne_am_am_soliddiffeq', ...
                     'pe_am_massCons'  , ...
                     'pe_am_chargeCons', ...
                     'pe_am_am_soliddiffeq', ...
                     'ne_cc_chargeCons' , ...
                     'pe_cc_chargeCons' , ...
                     'energyCons'       , ...
                     'EIeq', ...
                     'controlEq'};
            neq=numel(eqs);     
            keep = model.getEquationsToUses(neq);     
 
            if(not(all(keep)))
                ind   = find(keep);
                eqs   = {eqs{ind}};
                types = {types{ind}};
                names = {names{ind}};
            end
            
            isfixed = false;
            if isfixed
                switch ctrltype
                  case 'I'
                    types{end-1} = 'cell';   
                  case 'E'
                    neqs  = numel(types);
                    order = [1:neqs-2,neqs,neqs-1];
                    types = { types{order} };
                    eqs   = {eqs{order}};
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

            elyte   = 'Electrolyte';
            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            am      = 'ActiveMaterial';
            cc      = 'CurrentCollector';
            thermal = 'ThermalModel';
            
            % (here we assume that the ThermalModel has the "parent" grid)
            state.(elyte).T   = state.(thermal).T(model.(elyte).G.mappings.cellmap);
            state.(ne).(am).T = state.(thermal).T(model.(ne).(am).G.mappings.cellmap);
            state.(ne).(cc).T = state.(thermal).T(model.(ne).(cc).G.mappings.cellmap);
            state.(pe).(am).T = state.(thermal).T(model.(pe).(am).G.mappings.cellmap);
            state.(pe).(cc).T = state.(thermal).T(model.(pe).(cc).G.mappings.cellmap);
            
            % Update temperature in the active materials of the electrodes.
            state.(ne).(am) = model.(ne).(am).dispatchTemperature(state.(ne).(am));
            state.(pe).(am) = model.(pe).(am).dispatchTemperature(state.(pe).(am));
            
        end
        
        function keep = getEquationsToUses(model,neq)
            keep = true(neq, 1);     
            if(not(model.use_solid_diffusion))
                keep(5) = false;
                keep(8) = false;
            end
            if(not(model.use_thermal))
                keep(11) = false;
            end
        end
        
        function keep = getVariablesToUses(model, neq)
             %% reduction will not work if not this is equal to equations: if need one need to change LinearSolverAD.m
             keep = model.getEquationsToUses(neq);
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
            
            ne_R = state.(ne).(am).(itf).R;
            coupterm = getCoupTerm(couplingterms, 'NegativeElectrode-Electrolyte', coupnames);
            elytecells = coupterm.couplingcells(:, 2);
            elyte_c_source(elytecells) = ne_R.*vols(elytecells); % we divide with F later
            
            pe_R = state.(pe).(am).(itf).R;
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
            
            cdotcc  = (state.(elyte).c - state0.(elyte).c)/dt;
            effectiveVolumes = model.(elyte).volumeFraction.*model.(elyte).G.cells.volumes;
            massAccum  = effectiveVolumes.*cdotcc;
            state.(elyte).massAccum = massAccum;
            
            names = {ne, pe};
            for i = 1 : numel(names)
                elde = names{i}; % electrode name
                cdotcc   = (state.(elde).(am).c - state0.(elde).(am).c)/dt;
                effectiveVolumes = model.(elde).(am).volumeFraction.*model.(elde).(am).G.cells.volumes;
                massAccum  = effectiveVolumes.*cdotcc;
                state.(elde).(am).massAccum = massAccum;
            end
            
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
            am      = 'ActiveMaterial';
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
                cc_j     = locstate.(elde).(cc).jFace;
                cc_econd = cc_model.EffectiveElectricalConductivity;
                cc_vols  = cc_model.G.cells.volumes;
                cc_jsq   = computeCellFluxNorm(cc_model, cc_j); 
                state.(elde).(cc).jsq = cc_jsq;  %store square of current density
                
                src(cc_map) = src(cc_map) + cc_vols.*cc_jsq./cc_econd;

                am_model = model.(elde).(am);
                am_map   = am_model.G.mappings.cellmap;
                am_j     = locstate.(elde).(am).jFace;
                am_econd = am_model.EffectiveElectricalConductivity;
                am_vols  = am_model.G.cells.volumes;
                am_jsq   = computeCellFluxNorm(am_model, am_j);
                state.(elde).(am).jsq = am_jsq;
                
                src(am_map) = src(am_map) + am_vols.*am_jsq./am_econd;
                
            end

            % Electrolyte
            elyte_model = model.(elyte);
            elyte_map   = elyte_model.G.mappings.cellmap;
            elyte_vf    = elyte_model.volumeFraction;
            elyte_j     = locstate.(elyte).jFace;
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
                vols    = model.(elde).(am).G.cells.volumes;

                R      = locstate.(elde).(am).(itf).R;
                dUdT   = locstate.(elde).(am).(itf).dUdT;
                eta    = locstate.(elde).(am).(itf).eta;
                
                itf_src = n*F*vols.*R.*(eta + T(itf_map).*dUdT);
                
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
            cc = 'CurrentCollector';
            
            phi = state.(ne).(cc).phi;

            [jExternal, jFaceExternal] = model.(ne).(cc).setupExternalCoupling(phi, 0);
            
            state.(ne).(cc).jExternal = jExternal;
            state.(ne).(cc).jFaceExternal = jFaceExternal;
            
        end
        
        function state = setupExternalCouplingPositiveElectrode(model, state)
        %
        % Setup external electronic coupling of the positive electrode at the current collector
        %            
            pe   = 'PositiveElectrode';
            cc   = 'CurrentCollector';
            ctrl = 'Control';
            
            phi = state.(pe).(cc).phi;
            E   = state.(ctrl).E;
            
            [jExternal, jFaceExternal] = model.(pe).(cc).setupExternalCoupling(phi, E);
            
            state.(pe).(cc).jExternal = jExternal;
            state.(pe).(cc).jFaceExternal = jFaceExternal;
            
        end

        function state = updateDerivativeControlValues(model, state, state0, dt)

            ctrl = 'Control';
            
            dEdt = (state.(ctrl).E - state0.(ctrl).E)/dt;
            dIdt = (state.(ctrl).I - state0.(ctrl).I)/dt;
            
            state.(ctrl).dEdt = value(dEdt);
            state.(ctrl).dIdt = value(dIdt);
            state.(ctrl).prevCtrlType = state0.(ctrl).ctrlType;
            
        end
        
        function state = setupEIEquation(model, state)
            
            pe   = 'PositiveElectrode';
            cc   = 'CurrentCollector';
            ctrl = 'Control';
            
            I = state.(ctrl).I;
            E = state.(ctrl).E;
            phi = state.(pe).(cc).phi;
            
            coupterm = model.(pe).(cc).couplingTerm;
            faces    = coupterm.couplingfaces;
            cond_pcc = model.(pe).(cc).EffectiveElectricalConductivity;
            [trans_pcc, cells] = model.(pe).(cc).operators.harmFaceBC(cond_pcc, faces);
            
            state.Control.EIequation = sum(trans_pcc.*(state.(pe).(cc).phi(cells) - E)) - I;

        end
        
        function state = initStateAD(model, state)

            [pnames, extras]  = model.getPrimaryVariables();
            vars = cell(numel(pnames),1);
            for i=1:numel(pnames)
                vars{i} = model.getProp(state,pnames{i});
            end
            % Get the AD state for this model           
            [vars{:}] = model.AutoDiffBackend.initVariablesAD(vars{:});
            newstate =struct();
            for i=1:numel(pnames)
               newstate = model.setNewProp(newstate, pnames{i}, vars{i});
            end
            
            for i = 1 : numel(extras)
                var = model.getProp(state,extras{i});
                assert(isnumeric(var) | ischar(var));
                newstate = model.setNewProp(newstate, extras{i}, var);
            end
            
            time = state.time;
            state = newstate;
            state.time = time;
            
        end 
        
        function [p, extra] = getPrimaryVariables(model)
            
            bat = model;
            elyte   = 'Electrolyte';
            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            am      = 'ActiveMaterial';
            sd      = 'SolidDiffusion';
            cc      = 'CurrentCollector';
            ctrl    = 'Control';
            thermal = 'ThermalModel';
            
            %% chaged order of c cElectrode to get easier reduction
            p = {{elyte, 'c'}             , ...
                 {elyte, 'phi'}           , ...   
                 {ne, am, sd, 'cSurface'} , ...    
                 {ne, am, 'phi'}          , ...   
                 {ne, am, 'c'}            , ...
                 {pe, am, sd, 'cSurface'} , ...    
                 {pe, am, 'phi'}          , ...   
                 {pe, am, 'c'}            , ...
                 {ne, cc, 'phi'}          , ...    
                 {pe, cc, 'phi'}          , ...
                 {thermal, 'T'}           , ...
                 {ctrl, 'E'}              , ...
                 {ctrl, 'I'}};          
            
            neq = numel(p);
            keep = model.getVariablesToUses(neq);
            extra = {p{find(not(keep))}};
            if (not(all(keep)))
                p = {p{find(keep)}};
            end
            extra{end + 1} = {ctrl, 'ctrlType'};
            if strcmp(model.(ctrl).controlPolicy, 'CCCV')
                extra{end + 1} = {ctrl, 'nextCtrlType'};
            end
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
              otherwise
                error('Error controlPolicy not recognized');
            end
            
        end

        function model = validateModel(model, varargin)
            mnames = {{'Electrolyte'}, ...
                      {'PositiveElectrode','ActiveMaterial'}, ...
                      {'NegativeElectrode','ActiveMaterial'}, ...
                      {'NegativeElectrode','CurrentCollector'}, ...
                      {'PositiveElectrode','CurrentCollector'}};
            model.Electrolyte.AutoDiffBackend=model.AutoDiffBackend;
            model.Electrolyte=model.Electrolyte.validateModel(varargin{:});
            model.PositiveElectrode.ActiveMaterial.AutoDiffBackend= model.AutoDiffBackend;
            model.PositiveElectrode.ActiveMaterial = model.PositiveElectrode.ActiveMaterial.validateModel(varargin{:});
            model.NegativeElectrode.ActiveMaterial.AutoDiffBackend= model.AutoDiffBackend;
            model.NegativeElectrode.ActiveMaterial = model.NegativeElectrode.ActiveMaterial.validateModel(varargin{:});
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
            am    = 'ActiveMaterial';
            itf   = 'Interface';

            cmin = model.cmin;
            
            state.(elyte).c = max(cmin, state.(elyte).c);
            
            eldes = {ne, pe};
            for ind = 1 : numel(eldes)
                elde = eldes{ind};
                state.(elde).(am).c = max(cmin, state.(elde).(am).c);
                cmax = model.(elde).(am).(itf).cmax;
                state.(elde).(am).c = min(cmax, state.(elde).(am).c);
            end
            
            report = [];
            
        end

        function cleanState = addVariable(model, cleanState, state, state0)
            
            cleanState = addVariable@BaseModel(model, cleanState, state, state0);

            ctrl = 'Control';            
            cleanState.(ctrl).ctrlType = state.(ctrl).ctrlType;
            
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
        
        
        function outputvars = extractGlobalVariables(model, states)
            ns = numel(states);
            ws = cell(ns, 1);
            for i = 1 : ns
                E    = states{i}.Control.E;
                I    = states{i}.Control.I;
                T    = states{i}.ThermalModel.T;
                time = states{i}.time;
                
                Tmax = max(T);
                outputvars{i} = struct('E'   , E   , ...
                                       'I'   , I   , ...
                                       'Tmax', Tmax, ...
                                       'time', time);
            end
        end

    end
    
end



%{
Copyright 2009-2021 SINTEF Industry, Sustainable Energy Technology
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
