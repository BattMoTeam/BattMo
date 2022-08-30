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
        
        % flag that decide the model setup
        use_solid_diffusion
        use_thermal
        include_current_collectors
        
        primaryVariableNames
        addedVariableNames
        selectedEquationNames
        selectedEquationTypes
        
    end
    
    methods
        
        function model = Battery(paramobj)
            
            model = model@BaseModel();
            
            % All the submodels should have same backend (this is not assigned automaticallly for the moment)
            model.AutoDiffBackend = SparseAutoDiffBackend('useBlocks', false);
            
            %% Setup the model using the input parameters
            fdnames = {'G'            , ...
                       'couplingTerms', ...
                       'initT'        , ...
                       'use_thermal'  , ...
                       'include_current_collectors' , ...
                       'use_thermal'               , ...
                       'use_solid_diffusion'       , ...
                       'SOC'};
            
            model = dispatchParams(model, paramobj, fdnames);

            % Assign the components : Electrolyte, NegativeElectrode, PositiveElectrode
            params.include_current_collectors = model.include_current_collectors;
            params.use_thermal = model.use_thermal;

            model.NegativeElectrode = model.setupElectrode(paramobj.NegativeElectrode, params);
            model.PositiveElectrode = model.setupElectrode(paramobj.PositiveElectrode, params);
            model.Electrolyte       = model.setupElectrolyte(paramobj.Electrolyte);
            if model.use_thermal
                model.ThermalModel = ThermalComponent(paramobj.ThermalModel);
            end
            model.Control = model.setupControl(paramobj.Control);
            
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

            if model.use_thermal
                % setup Thermal Model by assigning the effective heat capacity and conductivity, which is computed from the sub-models.
                model = model.setupThermalModel();
            end
            
            % setup couplingNames
            model.couplingNames = cellfun(@(x) x.name, model.couplingTerms, 'uniformoutput', false);
            
            % setup equations and variable names selected in the model
            model = model.setupSelectedModel();
            
            % setup some mappings (mappings from electrodes to electrolyte)
            model = model.setupMappings();
            
            % setup capping
            cmax_ne = model.(ne).(am).(itf).cmax;
            cmax_pe = model.(pe).(am).(itf).cmax;
            model.cmin = 1e-5*max(cmax_ne, cmax_pe);
            
        end

        
        function model = setupSelectedModel(model)

            elyte   = 'Electrolyte';
            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            am      = 'ActiveMaterial';
            sd      = 'SolidDiffusion';
            cc      = 'CurrentCollector';
            ctrl    = 'Control';
            thermal = 'ThermalModel';
            
            % IMPORTANT : In the iterative solver, we remove parts of the linear systems. There, we assume a certain
            % order of the equations. The order is as indicated in the following cell structure: There is a given match
            % between the variables names (first column) and the equation names (second column). NOTE : this order MUST
            % be respected in the function getEquations where the equations are collected; nothing is imposed there.
            
            allNames = { ...
                {elyte, 'c'}            , 'elyte_massCons'      , 'cell' ; ... 
                {elyte, 'phi'}          , 'elyte_chargeCons'    , 'cell' ; ... 
                {ne, am, 'phi'}         , 'ne_am_chargeCons'    , 'cell' ; ... 
                {pe, am, 'phi'}         , 'pe_am_chargeCons'    , 'cell' ; ... 
                {ne, am, sd, 'cSurface'}, 'ne_am_massCons'      , 'scell' ; ... 
                {ne, am, sd, 'c'}       , 'ne_am_sd_massCons'   , 'sdiff'; ... 
                {ne, am, 'c'}           , 'ne_am_sd_soliddiffeq', 'sdiff'; ... 
                {pe, am, sd, 'cSurface'}, 'pe_am_massCons'      , 'scell' ; ...
                {pe, am, sd, 'c'}       , 'pe_am_sd_massCons'   , 'sdiff'; ... 
                {pe, am, 'c'}           , 'pe_am_sd_soliddiffeq', 'sdiff'; ... 
                {ne, cc, 'phi'}         , 'ne_cc_chargeCons'    , 'cell' ; ... 
                {pe, cc, 'phi'}         , 'pe_cc_chargeCons'    , 'cell' ; ... 
                {thermal, 'T'}          , 'energyCons'          , 'cell' ; ...
                {ctrl, 'E'}             , 'EIeq'                , 'ctrl' ; ...
                {ctrl, 'I'}             , 'controlEq'           , 'ctrl' };           
                
            allPrimaryVarNames = allNames(:, 1);
            allEquationNames   = allNames(:, 2);
            allEquationTypes   = allNames(:, 3);
                
            addedVariableNames = {};
            
            selectedvarinds = false(numel(allPrimaryVarNames), 1);
            selectedEqnInds = false(numel(allEquationNames), 1);
            
            % default set of primary variables (enter in all setups)
            pickVarInd = @(varnames) Battery.pickVarInd(varnames, allPrimaryVarNames);
            inds = pickVarInd({{elyte, 'c'}     , ...
                                {elyte, 'phi'}  , ...   
                                {ne, am, 'phi'} , ...   
                                {pe, am, 'phi'} , ...   
                                {ctrl, 'E'}     , ...
                                {ctrl, 'I'}});
            selectedVarInds(inds) = true;
            
            % default set of equations (enter in all setups)
            [isok, inds] = ismember({'elyte_massCons'   , ...
                                     'elyte_chargeCons' , ...
                                     'ne_am_chargeCons' , ...
                                     'pe_am_chargeCons' , ...
                                     'EIeq'             , ...
                                     'controlEq'}, allEquationNames);
            selectedEqnInds(inds) = true;
            
            if model.use_thermal
                selectedVarInds(pickVarInd({thermal, 'T'})) = true;
                [isok, inds] = ismember('energyCons', allEquationNames);
                selectedEqnInds(inds) = true;
            else
                addedVariableNames{end + 1} = {thermal, 'T'};
            end

            if model.use_solid_diffusion
                selectedVarInds(pickVarInd({{ne, am, sd, 'cSurface'}, {pe, am, sd, 'cSurface'}})) = true;
                [isok, inds] = ismember({'ne_am_sd_soliddiffeq', 'pe_am_sd_soliddiffeq'}, allEquationNames);
                selectedEqnInds(inds) = true;
                if model.(ne).(am).useSimplifiedDiffusionModel
                    selectedVarInds(pickVarInd({ne, am, 'c'})) = true;
                    [isok, inds] = ismember('ne_am_massCons', allEquationNames);
                    selectedEqnInds(inds) = true;                    
                else
                    selectedVarInds(pickVarInd({ne, am, sd, 'c'})) = true;
                    [isok, inds] = ismember('ne_am_sd_massCons', allEquationNames);
                    selectedEqnInds(inds) = true;                    
                end
                if model.(pe).(am).useSimplifiedDiffusionModel
                    selectedVarInds(pickVarInd({pe, am, 'c'})) = true;
                    [isok, inds] = ismember('pe_am_massCons', allEquationNames);
                    selectedEqnInds(inds) = true;                    
                else
                    selectedVarInds(pickVarInd({pe, am, sd, 'c'})) = true;
                    [isok, inds] = ismember('pe_am_sd_massCons', allEquationNames);
                    selectedEqnInds(inds) = true;                    
                end
            end

            if model.include_current_collectors
                selectedVarInds(pickVarInd({{ne, cc, 'phi'}, {pe, cc, 'phi'}})) = true;
                [isok, inds] = ismember({'ne_cc_chargeCons', 'pe_cc_chargeCons'}, allEquationNames);
                selectedEqnInds(inds) = true;
            end
            
            addedVariableNames{end + 1} = {ctrl, 'ctrlType'};
            
            if strcmp(model.(ctrl).controlPolicy, 'CCCV')
                addedVariableNames{end + 1} = {ctrl, 'nextCtrlType'};
            end
            
            model.primaryVariableNames  = allPrimaryVarNames(selectedVarInds);
            model.selectedEquationNames = allEquationNames(selectedEqnInds);
            model.selectedEquationTypes = allEquationTypes(selectedEqnInds);
            model.addedVariableNames    = addedVariableNames;
            
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
            
            % Dispatch electrolyte concentration and potential in the electrodes
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
            inputnames = {{pe, cc, 'E'}, ...
                          {pe, cc, 'I'}, ...
                          {pe, cc, 'phi'}, ...
                         };
            model = model.registerPropFunction({{'controlEq'}, fn, inputnames});

            % Function that update Thermal accumulation terms
            fn = @Battery.updateThermalAccumTerms;
            inputnames = {{thermal, 'T'}};
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
            
            %% Function that updates Thermal Reaction Terms
            fn = @Battery.updateThermalReactionSourceTerms;
            inputnames = {{ne, am, itf, 'R'}  , ...
                          {ne, am, itf, 'eta'}, ...
                          {pe, am, itf, 'R'}  , ...
                          {pe, am, itf, 'eta'}};
            model = model.registerPropFunction({{thermal, 'jHeatReactionSource'}, fn, inputnames});
                                                    
            %% Functions that setup external coupling at positive and negative electrodes
            
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
                
                if model.include_current_collectors
                    
                    % The effecive and intrinsic thermal parameters for the current collector are the same.
                    cc_map = model.(elde).(cc).G.mappings.cellmap;
                    cc_hcond = model.(elde).(cc).thermalConductivity;
                    cc_hcap = model.(elde).(cc).heatCapacity;

                    hcap(cc_map) = hcap(cc_map) + cc_hcap;
                    hcond(cc_map) = hcond(cc_map) + cc_hcond;
                    
                end
                
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

            if model.use_thermal
                model.ThermalModel.EffectiveHeatCapacity = hcap;
                model.ThermalModel.EffectiveThermalConductivity = hcond;
            end
            
        end
        
        
        function electrode = setupElectrode(model, paramobj, params)
        % Setup the electrode models (both :attr:`NegativeElectrode` and :attr:`PositiveElectrode`). Here, :code:`paramobj`
        % is instance of :class:`ElectrodeInputParams <Electrochemistry.Electrodes.ElectrodeInputParams>`
            electrode = Electrode(paramobj, params);
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
            model.(elyte).volumeFraction = subsasgnAD(model.(elyte).volumeFraction,elyte_cells(model.(ne).(am).G.mappings.cellmap), model.(ne).(am).porosity);
            model.(elyte).volumeFraction = subsasgnAD(model.(elyte).volumeFraction,elyte_cells(model.(pe).(am).G.mappings.cellmap), model.(pe).(am).porosity);
            sep_cells = elyte_cells(model.(elyte).(sep).G.mappings.cellmap); 
            model.(elyte).volumeFraction = subsasgnAD(model.(elyte).volumeFraction,sep_cells, model.(elyte).(sep).porosity);

            if model.use_thermal
                model.(elyte).EffectiveThermalConductivity = model.(elyte).volumeFraction.*model.(elyte).thermalConductivity;
            end

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
                
                m     = (1 ./ (elde_itf.theta100 - elde_itf.theta0));
                b     = -m .* elde_itf.theta0;
                theta = (SOC - b) ./ m;
                c     = theta .* elde_itf.cmax;
                nc    = elde_itf.G.cells.num;

                if model.(elde).(am).useSimplifiedDiffusionModel
                    initstate.(elde).(am).(sd).cSurface = c*ones(nc, 1);
                    initstate.(elde).(am).c = c*ones(nc, 1);
                else
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

            if model.(ne).include_current_collectors
                OCP = initstate.(ne).(am).(itf).OCP;
                OCP = OCP(1) .* ones(bat.(ne).(cc).G.cells.num, 1);
                initstate.(ne).(cc).phi = OCP - ref;
            end
            
            if model.(pe).include_current_collectors
                OCP = initstate.(pe).(am).(itf).OCP;
                OCP = OCP(1) .* ones(bat.(pe).(cc).G.cells.num, 1);
                initstate.(pe).(cc).phi = OCP - ref;
            end
            
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
                    state.(elde).(am).(itf).cElectrodeSurface = state.(elde).(am).c;
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
                state.(elde).(am).(itf) = battery.(elde).(am).(itf).updateEta(state.(elde).(am).(itf));
                state.(elde).(am).(itf) = battery.(elde).(am).(itf).updateReactionRate(state.(elde).(am).(itf));
            end

            %% Update Electrodes -> Electrolyte  coupling

            state = battery.updateElectrolyteCoupling(state);
            
            %% Update coupling within electrodes and external coupling
            
            for ind = 1 : numel(electrodes)
                elde = electrodes{ind};
                state.(elde).(am) = model.(elde).(am).updateConductivity(state.(elde).(am));
                if model.include_current_collectors
                    state.(elde).(cc) = model.(elde).(cc).updateConductivity(state.(elde).(cc));
                end
                state.(elde) = battery.(elde).updateCoupling(state.(elde));
                switch elde
                  case ne
                    state = model.setupExternalCouplingNegativeElectrode(state);
                  case pe
                    state = model.setupExternalCouplingPositiveElectrode(state);
                end
                if model.include_current_collectors
                    state.(elde).(cc) = battery.(elde).(cc).updatejBcSource(state.(elde).(cc));
                    state.(elde).(am) = battery.(elde).(am).updatejBcSource(state.(elde).(am));
                else
                    state.(elde).(am) = battery.(elde).(am).updatejBcSourceNoCurrentCollector(state.(elde).(am));
                end
                
            end
            
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

                if model.include_current_collectors
                    %% Electrodes charge conservation - current collector part
                    state.(elde).(cc) = battery.(elde).(cc).updateCurrent(state.(elde).(cc));
                    state.(elde).(cc) = battery.(elde).(cc).updateChargeConservation(state.(elde).(cc));
                end
                
            end
            
            %% elyte mass conservation

            state.(elyte) = battery.(elyte).updateDiffusionCoefficient(state.(elyte));
            state.(elyte) = battery.(elyte).updateMassFlux(state.(elyte));
            state.(elyte) = battery.(elyte).updateMassConservation(state.(elyte));

            %% update solid diffustion equations
            for ind = 1 : numel(electrodes)
                elde = electrodes{ind};
                state.(elde).(am).(sd) = battery.(elde).(am).(sd).updateDiffusionCoefficient(state.(elde).(am).(sd));
                state.(elde).(am) = battery.(elde).(am).dispatchRate(state.(elde).(am));
                if model.use_solid_diffusion
                    if model.(elde).(am).useSimplifiedDiffusionModel
                        state.(elde).(am) = battery.(elde).(am).assembleAccumTerm(state.(elde).(am), state0.(elde).(am), dt);
                        state.(elde).(am) = battery.(elde).(am).updateMassSource(state.(elde).(am));
                        state.(elde).(am) = battery.(elde).(am).updateMassConservation(state.(elde).(am));
                    else
                        state.(elde).(am).(sd) = battery.(elde).(am).(sd).updateMassSource(state.(elde).(am).(sd));
                        state.(elde).(am).(sd) = battery.(elde).(am).(sd).updateFlux(state.(elde).(am).(sd));
                        state.(elde).(am).(sd) = battery.(elde).(am).(sd).updateAccumTerm(state.(elde).(am).(sd), state0.(elde).(am).(sd), dt);
                        state.(elde).(am).(sd) = battery.(elde).(am).(sd).updateMassConservation(state.(elde).(am).(sd));
                    end
                    state.(elde).(am).(sd) = battery.(elde).(am).(sd).assembleSolidDiffusionEquation(state.(elde).(am).(sd));
                end
            end

            if model.use_thermal
                
                %% update Face fluxes
                for ind = 1 : numel(electrodes)
                    elde = electrodes{ind};
                    state.(elde).(am) = battery.(elde).(am).updateFaceCurrent(state.(elde).(am));
                    if model.include_current_collectors
                        state.(elde).(cc) = battery.(elde).(cc).updateFaceCurrent(state.(elde).(cc));
                    end
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
                
            end
            
            %% Setup relation between E and I at positive current collectror
            
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
            
            % Equation name : 'elyte_massCons';
            eqs{end + 1} = state.(elyte).massCons*massConsScaling;

            % Equation name : 'elyte_chargeCons';
            eqs{end + 1} = state.(elyte).chargeCons;
            
            % Equation name : 'ne_am_chargeCons';
            eqs{end + 1} = state.(ne).(am).chargeCons;

            % Equation name : 'pe_am_chargeCons';
            eqs{end + 1} = state.(pe).(am).chargeCons;
            
            if model.use_solid_diffusion
                
                if model.(ne).(am).useSimplifiedDiffusionModel
                    % Equation name : 'ne_am_massCons';
                    eqs{end + 1} = state.(ne).(am).massCons*massConsScaling;
                else
                    % Equation name : 'ne_am_sd_massCons';
                    n    = model.(ne).(am).(itf).n; % number of electron transfer (equal to 1 for Lithium)
                    F    = model.con.F;
                    vol  = model.(ne).(am).operators.pv;
                    rp   = model.(ne).(am).(sd).rp;
                    vsf  = model.(ne).(am).(sd).volumetricSurfaceArea;
                    surfp = 4*pi*rp^2;
                    
                    scalingcoef = (vsf*vol(1)*n*F)/surfp;
                    eqs{end + 1} = scalingcoef*state.(ne).(am).(sd).massCons;
                end
                
                % Equation name : 'ne_am_sd_soliddiffeq';
                eqs{end + 1} = state.(ne).(am).(sd).solidDiffusionEq.*massConsScaling.*battery.(ne).(am).(itf).G.cells.volumes/dt;
                
                if model.(pe).(am).useSimplifiedDiffusionModel
                    % Equation name : 'pe_am_massCons';
                    eqs{end + 1} = state.(pe).(am).massCons*massConsScaling;
                else
                    % Equation name : 'pe_am_sd_massCons';
                    n    = model.(ne).(am).(itf).n; % number of electron transfer (equal to 1 for Lithium)
                    F    = model.con.F;
                    vol  = model.(pe).(am).operators.pv;
                    rp   = model.(pe).(am).(sd).rp;
                    vsf  = model.(pe).(am).(sd).volumetricSurfaceArea;
                    surfp = 4*pi*rp^2;
                    
                    scalingcoef = (vsf*vol(1)*n*F)/surfp;
                    eqs{end + 1} = scalingcoef*state.(pe).(am).(sd).massCons;
                end
                
                % Equation name : 'pe_am_sd_soliddiffeq';
                eqs{end + 1} = state.(pe).(am).(sd).solidDiffusionEq.*massConsScaling.*battery.(pe).(am).(itf).G.cells.volumes/dt;
                
            end
            
            % Equation name : 'ne_cc_chargeCons';
            if model.(ne).include_current_collectors
                eqs{end + 1} = state.(ne).(cc).chargeCons;
            end
            
            % Equation name : 'pe_cc_chargeCons';
            if model.(pe).include_current_collectors
                eqs{end + 1} = state.(pe).(cc).chargeCons;
            end

            % Equation name : 'energyCons';
            if model.use_thermal
                eqs{end + 1} = state.(thermal).energyCons;
            end
            
            % Equation name : 'EIeq';
            eqs{end + 1} = - state.(ctrl).EIequation;
            
            % Equation name : 'controlEq'                                    
            eqs{end + 1} = state.(ctrl).controlEquation;
            
            eqs{1} = eqs{1} - model.Electrolyte.sp.t(1)*eqs{2};
            
            names = model.selectedEquationNames;
            types = model.selectedEquationTypes;
            
            %% The equations are reordered in a way that is consitent with the linear iterative solver 
            % (the order of the equation does not matter if we do not use an iterative solver)
            ctrltype = state.Control.ctrlType;
            switch ctrltype
              case {'constantCurrent', 'CC_charge', 'CC_discharge'}
                types{end - 1} = 'cell';   
              case {'constantVoltage', 'CV_charge'}
                neqs  = numel(types);
                order = [1 : neqs - 2, neqs, neqs - 1];
                types = {types{order}};
                eqs   = {eqs{order}};
                names = {names{order}};
              case 'CV_discharge'
                % stil I which is variable
                types{end - 1} = 'cell'; 
              otherwise 
                error()
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
            
            electrodes = {ne, pe};
            for i = 1 : numel(electrodes)
                elde = electrodes{i}; % electrode name
                if model.(elde).(am).useSimplifiedDiffusionModel | ~model.use_solid_diffusion
                    state.(elde).(am) = model.(elde).(am).assembleAccumTerm(state.(elde).(am), state0.(elde).(am), dt);
                end
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
            elyte_model = model.(elyte);
            elyte_map   = elyte_model.G.mappings.cellmap;
            elyte_vf    = elyte_model.volumeFraction;
            elyte_j     = state.(elyte).jFace;
            elyte_cond  = state.(elyte).conductivity;
            elyte_econd = elyte_cond.*elyte_vf.^1.5;
            elyte_vols  = elyte_model.G.cells.volumes;
            elyte_jsq   = computeCellFluxNorm(elyte_model, elyte_j);
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
            
            if model.(ne).include_current_collectors
                
                cc = 'CurrentCollector';

                phi   = state.(ne).(cc).phi;
                sigma = state.(ne).(cc).conductivity;

                [jExternal, jFaceExternal] = setupExternalCoupling(model.(ne).(cc), phi, 0, sigma);
                
                state.(ne).(cc).jExternal = jExternal;
                state.(ne).(cc).jFaceExternal = jFaceExternal;
                
            else
                
                am = 'ActiveMaterial';
                
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
            
            E   = state.(ctrl).E;

            if model.(pe).include_current_collectors
                
                cc   = 'CurrentCollector';
                
                phi   = state.(pe).(cc).phi;
                sigma = state.(pe).(cc).conductivity;
                
                [jExternal, jFaceExternal] = setupExternalCoupling(model.(pe).(cc), phi, E, sigma);
                
                state.(pe).(cc).jExternal = jExternal;
                state.(pe).(cc).jFaceExternal = jFaceExternal;
                
            else
                
                am = 'ActiveMaterial';
                
                phi   = state.(pe).(am).phi;
                sigma = state.(pe).(am).conductivity;
                
                [jExternal, jFaceExternal] = setupExternalCoupling(model.(pe).(am), phi, E, sigma);
                
                state.(pe).(am).jExternal = jExternal;
                state.(pe).(am).jFaceExternal = jFaceExternal;
                
            end
            
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
              otherwise
                error('Error controlPolicy not recognized');
            end
            %NB hack to get thing go
            forces.Imax = [];%model.Control.Imax;
            
        end

        function model = validateModel(model, varargin)

            model.PositiveElectrode.ActiveMaterial = model.PositiveElectrode.ActiveMaterial.setupDependentProperties();
            model.NegativeElectrode.ActiveMaterial = model.NegativeElectrode.ActiveMaterial.setupDependentProperties();
            model.Electrolyte.Separator = model.Electrolyte.Separator.setupDependentProperties();
            
            model = model.setupElectrolyteModel();
                
            model.PositiveElectrode.ActiveMaterial.AutoDiffBackend= model.AutoDiffBackend;
            model.PositiveElectrode.ActiveMaterial = model.PositiveElectrode.ActiveMaterial.validateModel(varargin{:});
            model.NegativeElectrode.ActiveMaterial.AutoDiffBackend= model.AutoDiffBackend;
            model.NegativeElectrode.ActiveMaterial = model.NegativeElectrode.ActiveMaterial.validateModel(varargin{:});
            model.Electrolyte.AutoDiffBackend=model.AutoDiffBackend;
            model.Electrolyte=model.Electrolyte.validateModel(varargin{:});
            if model.NegativeElectrode.include_current_collectors
                model.NegativeElectrode.CurrentCollector.AutoDiffBackend= model.AutoDiffBackend;
                model.NegativeElectrode.CurrentCollector= model.NegativeElectrode.CurrentCollector.validateModel(varargin{:});
            end
            
            if model.PositiveElectrode.include_current_collectors
                model.PositiveElectrode.CurrentCollector.AutoDiffBackend=model.AutoDiffBackend;
                model.PositiveElectrode.CurrentCollector= model.PositiveElectrode.CurrentCollector.validateModel(varargin{:});
            end
            
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
                if model.(elde).(am).useSimplifiedDiffusionModel | ~model.use_solid_diffusion
                    state.(elde).(am).c = max(cmin, state.(elde).(am).c);
                    cmax = model.(elde).(am).(itf).cmax;
                    state.(elde).(am).c = min(cmax, state.(elde).(am).c);
                else
                    state.(elde).(am).(sd).c = max(cmin, state.(elde).(am).(sd).c);
                    cmax = model.(elde).(am).(itf).cmax;
                    state.(elde).(am).(sd).c = min(cmax, state.(elde).(am).(sd).c);
                end
            end
            
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
        
        function ind = pickVarInd(varnames, givenvarnames)
            for ivarnames = 1 : numel(givenvarnames)
                givenvarnames{ivarnames} = strjoin(givenvarnames{ivarnames}, '_');
            end
            if ischar(varnames{1})
                varnames = {varnames};
            end
            for ivarnames = 1 : numel(varnames)
                varnames{ivarnames} = strjoin(varnames{ivarnames}, '_');
            end
            [isok, ind] = ismember(varnames, givenvarnames);
            
            assert(all(isok), 'Some variable names are not been found');
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
