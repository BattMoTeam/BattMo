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

        NegativeElectrode % Negative Electrode Model, instance of :class:`Electrode <Electrochemistry.Electrodes.Electrode>`
        PositiveElectrode % Positive Electrode Model, instance of :class:`Electrode <Electrochemistry.Electrodes.Electrode>`
        Electrolyte       % Electrolyte model, instance of :class:`Electrolyte <Electrochemistry.Electrodes.Electrolyte>`
        Separator         % Separator model, instance of :class:`Separator <Electrochemistry.Electrodes.Separator>`
        ThermalModel      % Thermal model, instance of :class:`ThermalComponent <Electrochemistry.ThermalComponent>`
        Control           % Control Model

        SOC % State Of Charge

        initT % Initial temperature

        couplingTerms % Coupling terms
        cmin % mininum concentration used in capping

        couplingNames

        mappings

        % flag that decide the model setup
        use_thermal
        include_current_collectors

        funcCallList

        primaryVarNames
        equationVarNames
        equationNames
        equationTypes
        equationIndices

    end

    methods

        function model = Battery(paramobj)

            model = model@BaseModel();

            paramobj.validateInputParams();

            % All the submodels should have same backend (this is not assigned automaticallly for the moment)
            model.AutoDiffBackend = SparseAutoDiffBackend('useBlocks', false);

            %% Setup the model using the input parameters
            fdnames = {'G'                         , ...
                       'couplingTerms'             , ...
                       'initT'                     , ...
                       'use_thermal'               , ...
                       'include_current_collectors', ...
                       'use_thermal'               , ...
                       'SOC'};

            model = dispatchParams(model, paramobj, fdnames);

            model.NegativeElectrode = Electrode(paramobj.NegativeElectrode);
            model.PositiveElectrode = Electrode(paramobj.PositiveElectrode);
            model.Separator         = Separator(paramobj.Separator);

            % We setup the electrolyte model (in particular we compute the volume fraction from the other components)
            model = model.setupElectrolyteModel(paramobj);

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
            model = model.setupCapping();

            % setup computational graph
            model = model.setupComputationalGraph();
            model.funcCallList = model.computationalGraph.getOrderedFunctionCallList();
            
        end

        function model = setupSelectedModel(model, varargin)
        % The system of equation should fullfill a special structure to fit into the iterative linear solver with
        % preconditioner. We create this structure here.

            opt = struct('reduction', []);
            opt = merge_options(opt, varargin{:});
            % For the reduction structure format, see battmodDir()/Utilities/JsonSchemas/linearsolver.schema.json and the reduction property

            elyte   = 'Electrolyte';
            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            co      = 'Coating';
            am      = 'ActiveMaterial';
            am1     = 'ActiveMaterial1';
            am2     = 'ActiveMaterial2';
            sd      = 'SolidDiffusion';
            cc      = 'CurrentCollector';
            ctrl    = 'Control';
            thermal = 'ThermalModel';

            addedVarNames = {};

            varEqTypes ={{elyte, 'c'}   , {elyte,  'massCons'}     , 'cell'; ...
                         {elyte, 'phi'} , {elyte,  'chargeCons'}   , 'cell'; ...
                         {ne, co, 'phi'}, {ne, co, 'chargeCons'}   , 'cell'; ...
                         {pe, co, 'phi'}, {pe, co, 'chargeCons'}   , 'cell'; ...
                         {ctrl, 'E'}    , {ctrl, 'EIequation'}     , 'ctrl'; ...
                         {ctrl, 'I'}    , {ctrl, 'controlEquation'}, 'ctrl'};

            eldes = {ne, pe};

            
            for ielde = 1 : numel(eldes)

                elde = eldes{ielde};

                switch model.(elde).(co).active_material_type
                  case 'default'
                    ams = {am};
                  case 'composite'
                    ams = {am1, am2};
                  otherwise
                    error('active_material_type not recognized');
                end
                
                for iam = 1 : numel(ams)
                    amc = ams{iam};
                    switch model.(elde).(co).(amc).diffusionModelType
                      case 'simple'
                        newentries = {{elde, co, amc, sd, 'cAverage'}, {elde, co, amc, sd, 'massCons'}        , 'cell'; ...
                                      {elde, co, amc, sd, 'cSurface'}, {elde, co, amc, sd, 'solidDiffusionEq'}, 'cell'};
                      case 'full'
                        newentries = {{elde, co, amc, sd, 'c'}       , {elde, co, amc, sd, 'massCons'}        , 'cell'; ...
                                      {elde, co, amc, sd, 'cSurface'}, {elde, co, amc, sd, 'solidDiffusionEq'}, 'cell'};
                      otherwise
                        error('diffusionModelType not recognized');
                    end
                    
                    varEqTypes = vertcat(varEqTypes, newentries);

                end
                
                if model.include_current_collectors

                    newentries = {{elde, cc, 'phi'}, {elde, cc, 'chargeCons'}, 'cell'};
                    varEqTypes = vertcat(varEqTypes, newentries);

                end

            end

            if model.use_thermal

                newentries = {{thermal, 'T'}, {thermal, 'energyCons'}, 'cell'};
                varEqTypes = vertcat(varEqTypes, newentries);

            end

            primaryVarNames = varEqTypes(:, 1);
            equationTypes   = varEqTypes(:, 3);

            % The variable and equation lists are not a priori ordered (in the sense that we have 'cell' types first and
            % the other types after). It is a requirement in some setup of the linear solver.
            % Note : if you use a direct solver, this is not used.

            variableReordered = false;

            if ~isempty(opt.reduction)
                reduc = opt.reduction;
                % We set the type of the variable to be reduced as 'reduced' or 'specialReduced' (anything but 'cell')
                % and we move them at the end of list respecting the order in which they have been given.
                if ~isempty(reduc) && reduc.doReduction
                    neq = numel(equationTypes);
                    equationTypes = cell(neq, 1);
                    for ieqtype = 1 : neq
                        equationTypes{ieqtype} = 'cell';
                    end

                    variables = reduc.variables;
                    if isstruct(variables)
                        variables = num2cell(variables);
                    end
                    einds = nan(numel(variables),1);
                    for ivar = 1 : numel(variables)
                        var = variables{ivar};
                        [found, ind] = Battery.getVarIndex(var.name, primaryVarNames);
                        if ~found
                            error('variable to be reduce has not been found');
                        end
                        equationTypes{ind} = 'reduced';
                        if isfield(var, "special") && var.special
                            equationTypes{ind} = 'specialReduced';
                        end
                        einds(ivar) = ind;
                        order(ivar) = var.order;
                    end

                    [~, ind] = sort(order);
                    einds = einds(ind);

                    inds = (1 : neq)';
                    inds(einds) = [];
                    inds = [inds; einds];
                    variableReordered = true;

                end
            end

            if ~variableReordered
                % We reorder to get first the 'cells' type (required for reduction in iterative solver)
                iscell = ismember(equationTypes, {'cell'});
                inds = [find(iscell); find(~iscell)];
            end


            primaryVarNames  = varEqTypes(inds, 1);
            equationVarNames = varEqTypes(inds, 2);
            equationTypes    = equationTypes(inds);

            % We use shortened names for easier visualisation and because Matlab also has a limitation on the lenght of
            % a field name in a structure.
            function str = setupName(varname)
                shortvarname = cellfun(@(elt) Battery.shortenName(elt), varname, 'uniformoutput', false);
                str = Battery.varToStr(shortvarname);
            end
            equationNames = cellfun(@(varname) setupName(varname), equationVarNames, 'uniformoutput', false);

            equationIndices = struct();
            for ieq = 1 : numel(equationNames)
                equationIndices.(equationNames{ieq}) = ieq;
            end

            model.primaryVarNames  = primaryVarNames;
            model.equationVarNames = equationVarNames;
            model.equationNames    = equationNames;
            model.equationTypes    = equationTypes;
            model.equationIndices  = equationIndices;

        end

        function model = registerVarAndPropfuncNames(model)

            %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)

            model = registerVarAndPropfuncNames@BaseModel(model);

            % defines shorthands for the submodels
            elyte   = 'Electrolyte';
            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            co      = 'Coating';
            am      = 'ActiveMaterial';
            am1     = 'ActiveMaterial1';
            am2     = 'ActiveMaterial2';
            cc      = 'CurrentCollector';
            itf     = 'Interface';
            sd      = 'SolidDiffusion';
            thermal = 'ThermalModel';
            ctrl    = 'Control';

            if ~model.use_thermal
                % we register the temperature variable, as it is not done by the ThermalModel which is empty in this case
                model = model.registerVarName({thermal, 'T'});
            end

            eldes = {ne, pe};


            %% Temperature dispatch functions
            fn = @Battery.updateTemperature;

            inputnames = {{thermal, 'T'}};
            model = model.registerPropFunction({{elyte, 'T'}  , fn, inputnames});
            for ielde = 1 : numel(eldes)
                elde = eldes{ielde};
                model = model.registerPropFunction({{elde, co, 'T'} , fn, inputnames});
            end

            %% Coupling functions

            % Dispatch electrolyte concentration and potential in the electrodes
            for ielde = 1 : numel(eldes)
                
                elde = eldes{ielde};

                switch model.(elde).(co).active_material_type
                  case 'default'
                    ams = {am};
                  case 'composite'
                    ams = {am1, am2};
                  otherwise
                    error('active_material_type not recognized');
                end

                for iam = 1 : numel(ams)

                    amc = ams{iam};
                    
                    fn = @Battery.updateElectrodeCoupling;
                    inputnames = {{elyte, 'c'}, ...
                                  {elyte, 'phi'}};
                    
                    model = model.registerPropFunction({{elde, co, amc, itf, 'phiElectrolyte'}, fn, inputnames});
                    model = model.registerPropFunction({{elde, co, amc, itf, 'cElectrolyte'}  , fn, inputnames});


                end

            end

            fn = @Battery.updateElectrolyteCoupling;
            inputnames = {{ne, co, 'eSource'}, ...
                          {pe, co, 'eSource'}};
            model = model.registerPropFunction({{elyte, 'massSource'}, fn, inputnames});
            model = model.registerPropFunction({{elyte, 'eSource'}, fn, inputnames});
                    

            % Function that assemble the control equation
            fn = @Battery.setupEIEquation;
            inputnames = {{ctrl, 'E'}, ...
                          {ctrl, 'I'}};
            if model.include_current_collectors
                inputnames{end + 1} = {pe, cc, 'phi'};
            else
                inputnames{end + 1} = {pe, co, 'phi'};
            end
            model = model.registerPropFunction({{ctrl, 'EIequation'}, fn, inputnames});


            inputnames = {};
            fn = @Battery.updateControl;
            fn = {fn, @(propfunction) PropFunction.drivingForceFuncCallSetupFn(propfunction)};
            model = model.registerPropFunction({{ctrl, 'ctrlVal'}, fn, inputnames});
            model = model.registerPropFunction({{ctrl, 'ctrlType'}, fn, inputnames});


            %% Function that update the Thermal Ohmic Terms

            if model.use_thermal

                fn = @Battery.updateThermalOhmicSourceTerms;
                inputnames = {{elyte, 'jFace'}        , ...
                              {ne, co, 'jFace'}       , ...
                              {pe, co, 'jFace'}       , ...
                              {elyte, 'conductivity'} , ...
                              {ne, co, 'conductivity'}, ...
                              {pe, co, 'conductivity'}};

                if model.include_current_collectors
                    varnames ={{ne, cc, 'jFace'}       , ...
                               {pe, cc, 'jFace'}       , ...
                               {ne, cc, 'conductivity'}, ...
                               {pe, cc, 'conductivity'}};
                    inputnames = horzcat(inputnames, varnames);
                end

                model = model.registerPropFunction({{thermal, 'jHeatOhmSource'}, fn, inputnames});

                %% Function that updates the Thermal Chemical Terms
                fn = @Battery.updateThermalChemicalSourceTerms;
                inputnames = {{elyte, 'diffFlux'}, ...
                              {elyte, 'D'}       , ...
                              VarName({elyte}, 'dmudcs', 2)};
                model = model.registerPropFunction({{thermal, 'jHeatChemicalSource'}, fn, inputnames});

                %% Function that updates Thermal Reaction Terms
                fn = @Battery.updateThermalReactionSourceTerms;
                inputnames = {{thermal, 'T'}           , ...
                              {ne, co, am, sd, 'Rvol'} , ...
                              {ne, co, am, itf, 'eta'} , ...
                              {ne, co, am, itf, 'dUdT'}, ...
                              {pe, co, am, sd, 'Rvol'} , ...
                              {pe, co, am, itf, 'eta'} , ...
                              {pe, co, am, itf, 'dUdT'}};
                model = model.registerPropFunction({{thermal, 'jHeatReactionSource'}, fn, inputnames});

            else
                model = model.removeVarName({elyte, 'diffFlux'});
                model = model.removeVarName({ne, co, am, itf, 'dUdT'});
                model = model.removeVarName({pe, co, am, itf, 'dUdT'});
            end

            %% Functions that setup external  coupling for negative electrode


            fns{1} = @Battery.setupExternalCouplingNegativeElectrode;
            fns{2} = @Battery.setupExternalCouplingPositiveElectrode;

            for ielde = 1 : numel(eldes)

                elde = eldes{ielde};
                fn   = fns{ielde};

                if model.include_current_collectors

                    inputnames = {{elde, cc, 'phi'}, ...
                                  {elde, cc, 'conductivity'}};
                    model = model.registerPropFunction({{elde, cc, 'jExternal'}, fn, inputnames});

                    if model.use_thermal
                        model = model.registerPropFunction({{elde, cc, 'jFaceExternal'}, fn, inputnames});
                    end

                    inputnames = {{elde, co, 'phi'}, ...
                                  {elde, co, 'conductivity'}};
                    model = model.registerPropFunction({{elde, co, 'jExternal'}, fn, inputnames});

                else

                    inputnames = {{elde, co, 'phi'}, ...
                                  {elde, co, 'conductivity'}};
                    model = model.registerPropFunction({{elde, co, 'jExternal'}, fn, inputnames});

                    if model.use_thermal
                        model = model.registerPropFunction({{elde, co, 'jFaceExternal'}, fn, inputnames});
                    end

                end

            end

            %% Declare the "static" variables
            varnames = {};
            if ~model.use_thermal
                varnames{end + 1} = {thermal, 'T'};
            end
            model = model.registerStaticVarNames(varnames);


        end

        function control = setupControl(model, paramobj)

            C = computeCellCapacity(model);

            switch paramobj.controlPolicy
                
              case "CCDischarge"

                control = CCDischargeControlModel(paramobj);
                CRate = control.CRate;
                control.Imax = (C/hour)*CRate;
                
              case 'CCCharge'

                control = CCChargeControlModel(paramobj);
                CRate = control.CRate;
                control.Imax = (C/hour)*CRate;
                
              case "CCCV"
                
                control = CcCvControlModel(paramobj);
                CRate = control.CRate;
                control.Imax = (C/hour)*CRate;
                
              case "powerControl"
                
                control = PowerControlModel(paramobj);
                
              case "CC"
                
                control = CcControlModel(paramobj);
                CRate = control.CRate;
                control.Imax = (C/hour)*CRate;
                
              otherwise
                
                error('Error controlPolicy not recognized');
            end

        end

        function model = setupThermalModel(model, paramobj)
        % Setup the thermal model :attr:`ThermalModel`. Here, :code:`paramobj` is instance of
        % :class:`ThermalComponentInputParams <Electrochemistry.ThermalComponentInputParams>`

            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            co      = 'Coating';
            am      = 'ActiveMaterial';
            cc      = 'CurrentCollector';
            elyte   = 'Electrolyte';
            sep     = 'Separator';
            thermal = 'ThermalModel';

            eldes = {ne, pe}; % electrodes

            G = model.G;
            nc = G.cells.num;

            vhcap = zeros(nc, 1); % effective heat capacity
            hcond = zeros(nc, 1); % effective heat conductivity

            for ind = 1 : numel(eldes)

                elde = eldes{ind};

                if model.include_current_collectors

                    % The effecive and intrinsic thermal parameters for the current collector are the same.
                    cc_map = model.(elde).(cc).G.mappings.cellmap;

                    cc_hcond = model.(elde).(cc).effectiveThermalConductivity;
                    cc_vhcap = model.(elde).(cc).effectiveVolumetricHeatCapacity;

                    hcond(cc_map) = hcond(cc_map) + cc_hcond;
                    vhcap(cc_map) = vhcap(cc_map) + cc_vhcap;

                end

                % Effective parameters from the Electrode Active Component region.
                co_map       = model.(elde).(co).G.mappings.cellmap;

                co_hcond     = model.(elde).(co).effectiveThermalConductivity;
                co_vhcap     = model.(elde).(co).effectiveVolumetricHeatCapacity;

                hcond(co_map) = hcond(co_map) + co_hcond;
                vhcap(co_map) = vhcap(co_map) + co_vhcap;

            end

            % Electrolyte

            elyte_map   = model.(elyte).G.mappings.cellmap;

            elyte_hcond = model.(elyte).effectiveThermalConductivity;
            elyte_vhcap = model.(elyte).effectiveVolumetricHeatCapacity;

            vhcap(elyte_map) = vhcap(elyte_map) + elyte_vhcap;
            hcond(elyte_map) = hcond(elyte_map) + elyte_hcond;

            % Separator

            sep_map   = model.(sep).G.mappings.cellmap;

            sep_hcond = model.(sep).effectiveThermalConductivity;
            sep_vhcap = model.(sep).effectiveVolumetricHeatCapacity;

            vhcap(sep_map) = vhcap(sep_map) + sep_vhcap;
            hcond(sep_map) = hcond(sep_map) + sep_hcond;

            % Assign values

            model.(thermal).effectiveVolumetricHeatCapacity = vhcap;
            model.(thermal).effectiveThermalConductivity = hcond;


        end


        function model = setupMappings(model)

            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            co    = 'Coating';
            cc    = 'CurrentCollector';
            elyte = 'Electrolyte';

            G_elyte = model.(elyte).G;
            elytecelltbl.cells = (1 : G_elyte.cells.num)';
            elytecelltbl.globalcells = G_elyte.mappings.cellmap;
            elytecelltbl = IndexArray(elytecelltbl);

            eldes = {ne, pe};

            for ind = 1 : numel(eldes)

                elde = eldes{ind};
                G_elde  = model.(elde).(co).G;
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

        function model = setupElectrolyteModel(model, paramobj)
        % Assign the electrolyte volume fractions in the different regions

            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            co    = 'Coating';
            sep   = 'Separator';

            elyte_cells = zeros(model.G.cells.num, 1);
            elyte_cells(paramobj.(elyte).G.mappings.cellmap) = (1 : paramobj.(elyte).G.cells.num)';

            paramobj.(elyte).volumeFraction = ones(paramobj.(elyte).G.cells.num, 1);
            paramobj.(elyte).volumeFraction = subsasgnAD(paramobj.(elyte).volumeFraction, elyte_cells(model.(ne).(co).G.mappings.cellmap), 1 - model.(ne).(co).volumeFraction);
            paramobj.(elyte).volumeFraction = subsasgnAD(paramobj.(elyte).volumeFraction, elyte_cells(model.(pe).(co).G.mappings.cellmap), 1 - model.(pe).(co).volumeFraction);
            paramobj.(elyte).volumeFraction = subsasgnAD(paramobj.(elyte).volumeFraction, elyte_cells(model.(sep).G.mappings.cellmap), model.(sep).porosity);

            model.(elyte) = Electrolyte(paramobj.(elyte));

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
            co      = 'Coating';
            am      = 'ActiveMaterial';
            am1     = 'ActiveMaterial1';
            am2     = 'ActiveMaterial2';
            itf     = 'Interface';
            sd      = 'SolidDiffusion';
            cc      = 'CurrentCollector';
            thermal = 'ThermalModel';
            ctrl    = 'Control';

            initstate.(thermal).T = T*ones(nc, 1);

            %% Synchronize temperatures
            initstate = model.updateTemperature(initstate);

            %% Setup initial state for electrodes

            eldes = {ne, pe};

            for ielde = 1 : numel(eldes)

                elde = eldes{ielde};
                
                switch model.(elde).(co).active_material_type
                  case 'default'
                    ams = {am};
                  case 'composite'
                    ams = {am1, am2};
                  otherwise
                    error('active_material_type not recognized');
                end

                for iam = 1 : numel(ams)

                    amc = ams{iam};
                    
                    elde_itf = bat.(elde).(co).(amc).(itf);

                    theta = SOC*(elde_itf.guestStoichiometry100 - elde_itf.guestStoichiometry0) + elde_itf.guestStoichiometry0;
                    c     = theta*elde_itf.saturationConcentration;
                    nc    = model.(elde).(co).G.cells.num;

                    switch model.(elde).(co).(amc).diffusionModelType
                      case 'simple'
                        initstate.(elde).(co).(amc).(sd).cSurface = c*ones(nc, 1);
                        initstate.(elde).(co).(amc).(sd).cAverage = c*ones(nc, 1);
                      case 'full'
                        initstate.(elde).(co).(amc).(sd).cSurface = c*ones(nc, 1);
                        N = model.(elde).(co).(amc).(sd).N;
                        np = model.(elde).(co).(amc).(sd).np; % Note : we have by construction np = nc
                        initstate.(elde).(co).(amc).(sd).c = c*ones(N*np, 1);
                      otherwise
                        error('diffusionModelType not recognized')
                    end

                end

                % In case of two material, we choose first material for the initial OCP
                amc = ams{1};
                
                initstate = model.evalVarName(initstate, {elde, co, amc, itf, 'OCP'});

                OCP = initstate.(elde).(co).(amc).(itf).OCP;
                
                if ielde == 1
                    % The value in the first cell is used as reference.
                    ref = OCP(1);
                end

                initstate.(elde).(co).phi = OCP - ref;

            end

            %% Setup initial Electrolyte state

            initstate.(elyte).phi = zeros(bat.(elyte).G.cells.num, 1) - ref;
            initstate.(elyte).c   = 1*mol/litre*ones(bat.(elyte).G.cells.num, 1);

            %% Setup initial Current collectors state

            if model.(ne).include_current_collectors
                OCP = initstate.(ne).(co).(amc).(itf).OCP;
                OCP = OCP(1) .* ones(bat.(ne).(cc).G.cells.num, 1);
                initstate.(ne).(cc).phi = OCP - ref;
            end

            if model.(pe).include_current_collectors
                OCP = initstate.(pe).(co).(amc).(itf).OCP;
                OCP = OCP(1) .* ones(bat.(pe).(cc).G.cells.num, 1);
                initstate.(pe).(cc).phi = OCP - ref;
            end

            initstate.(ctrl).E = OCP(1) - ref;

            switch model.(ctrl).controlPolicy
                
              case 'CCDischarge'
                
                initstate.(ctrl).ctrlType = 'constantCurrent';
                initstate.(ctrl).I = model.(ctrl).Imax;
                
              case 'CCCharge'
                
                initstate.(ctrl).ctrlType = 'constantCurrent';
                initstate.(ctrl).I = -model.(ctrl).Imax;
                
              case 'CCCV'
                
                switch model.(ctrl).initialControl
                  case 'discharging'
                    initstate.(ctrl).ctrlType     = 'CC_discharge1';
                    initstate.(ctrl).nextCtrlType = 'CC_discharge1';
                    initstate.(ctrl).I            = model.(ctrl).Imax;
                  case 'charging'
                    initstate.(ctrl).ctrlType     = 'CC_charge1';
                    initstate.(ctrl).nextCtrlType = 'CC_charge1';
                    initstate.(ctrl).I            = - model.(ctrl).Imax;
                  otherwise
                    error('initialControl not recognized');
                end
                
              case 'powerControl'
                
                switch model.(ctrl).initialControl
                  case 'discharging'
                    error('to implement (should be easy...)')
                  case 'charging'
                    initstate.(ctrl).ctrlType = 'charge';
                    E = initstate.(ctrl).E;
                    P = model.(ctrl).chargingPower;
                    initstate.(ctrl).I = -P/E;
                  otherwise
                    error('initialControl not recognized');
                end
                
              case 'CC'
                
                % this value will be overwritten after first iteration
                initstate.(ctrl).I = 0;
                switch model.(ctrl).initialControl
                  case 'discharging'
                    initstate.(ctrl).ctrlType = 'discharge';
                  case 'charging'
                    initstate.(ctrl).ctrlType = 'charge';
                  otherwise
                    error('initialControl not recognized');
                end
              otherwise
                error('control policy not recognized');
            end

            initstate.time = 0;

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
            paramobj = ControlModelInputParams([]);
            model.Control = ControlModel(paramobj);
            model.Control.controlPolicy = 'None';
            drivingForces = model.getValidDrivingForces();

            % We call getEquations to update state

            [~, state] = getEquations(model, state0, state, dt, drivingForces, 'ResOnly', true);

            % We set to empty the fields we know that are not meaningfull

            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            am      = 'ActiveMaterial';
            co      = 'Coating';
            cc      = 'CurrentCollector';
            elyte   = 'Electrolyte';
            am1     = 'ActiveMaterial1';
            am2     = 'ActiveMaterial2';
            itf     = 'Interface';
            sd      = "SolidDiffusion";
            thermal = 'ThermalModel';
            ctrl    = 'Control';

            state.(elyte).massAccum = [];
            state.(elyte).massCons = [];
            
            eldes = {ne, pe};
            
            for ielde = 1 : numel(eldes)

                elde = eldes{ielde};
                
                switch model.(elde).(co).active_material_type
                  case 'default'
                    ams = {am};
                  case 'composite'
                    ams = {am1, am2};
                  otherwise
                    error('active_material_type not recognized');
                end

                for iam = 1 : numel(ams)
                    
                    amc = ams{iam};
                    
                    state.(elde).(co).(amc).(sd).massAccum = [];
                    state.(elde).(co).(amc).(sd).massCons = [];

                end

                state = model.evalVarName(state, {elde, co, 'SOC'});
                
            end

            if model.use_thermal
                state = model.updateThermalIrreversibleReactionSourceTerms(state);
                state = model.updateThermalReversibleReactionSourceTerms(state);
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
               dispif(mrstVerbose, 'No AD initialization in equation old style')
               state0 = model.initStateAD(state0);
            else
                assert(opts.ResOnly);
            end


            %% We call the assembly equations ordered from the graph

            funcCallList = model.funcCallList;

            for ifunc = 1 : numel(funcCallList)
                eval(funcCallList{ifunc});
            end

            %% We apply some scaling

            % Shorthands used in this function
            battery = model;
            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            co      = 'Coating';
            am      = 'ActiveMaterial';
            am1     = 'ActiveMaterial1';
            am2     = 'ActiveMaterial2';
            elyte   = 'Electrolyte';
            am      = 'ActiveMaterial';
            itf     = 'Interface';
            sd      = "SolidDiffusion";
            thermal = 'ThermalModel';
            ctrl    = 'Control';

            eldes = {ne, pe};

            massConsScaling = model.con.F;

            state.(elyte).massCons = state.(elyte).massCons*massConsScaling;

            for ielde = 1 : numel(eldes)

                elde = eldes{ielde};

                switch model.(elde).(co).active_material_type
                  case 'default'
                    ams = {am};
                  case 'composite'
                    ams = {am1, am2};
                  otherwise
                    error('active_material_type not recognized');
                end

                for iam = 1 : numel(ams)
                    
                    amc = ams{iam};
                    
                    switch model.(elde).(co).(amc).diffusionModelType

                      case 'simple'

                        state.(elde).(co).(amc).(sd).massCons         = massConsScaling*state.(elde).(co).(amc).(sd).massCons;
                        state.(elde).(co).(amc).(sd).solidDiffusionEq = massConsScaling.*battery.(elde).(co).G.cells.volumes/dt.*state.(elde).(co).(amc).(sd).solidDiffusionEq;

                      case 'full'

                        n    = model.(elde).(co).(amc).(itf).numberOfElectronsTransferred;
                        F    = model.con.F;
                        vol  = model.(elde).(co).operators.pv;
                        rp   = model.(elde).(co).(amc).(sd).particleRadius;
                        vsf  = model.(elde).(co).(amc).(sd).volumetricSurfaceArea;

                        surfp = 4*pi*rp^2;

                        scalingcoef = (vsf*vol(1)*n*F)/surfp;

                        state.(elde).(co).(amc).(sd).massCons         = scalingcoef*state.(elde).(co).(amc).(sd).massCons;
                        state.(elde).(co).(amc).(sd).solidDiffusionEq = scalingcoef*state.(elde).(co).(amc).(sd).solidDiffusionEq;

                      otherwise

                        error('diffusionModelType not recognized');

                    end
                end
            end

            for ieq = 1 : numel(model.equationVarNames)
                eqs{ieq} = model.getProp(state, model.equationVarNames{ieq});
            end

            ei = model.equationIndices;

            % By doing this linear transformation, we remove the direct dependency of the mass conservation equation with
            % respect to the potential gradien. Only the concentration gradient remains in the equation. This is a special
            % property when the transferance t is a constant.

            mieq = ei.(Battery.varToStr({'elyte', 'massCons'}));
            cieq = ei.(Battery.varToStr({'elyte', 'chargeCons'}));

            eqs{mieq} = eqs{mieq} - model.(elyte).sp.t(1)*eqs{cieq};

            names       = model.equationNames;
            types       = model.equationTypes;
            primaryVars = model.primaryVarNames;

            %% The equations are reordered in a way that is consitent with the linear iterative solver
            % (the order of the equation does not matter if we do not use an iterative solver)
            ctrltype = state.Control.ctrlType;
            switch ctrltype

              case {'constantCurrent', 'CC_discharge1', 'CC_discharge2', 'CC_charge1', 'charge', 'discharge'}

                eqname = Battery.varToStr({'ctrl', 'EIequation'});
                types{ei.(eqname)} = 'cell';

              case {'constantVoltage', 'CV_charge2'}

                eieqname = Battery.varToStr({'ctrl', 'EIequation'});
                cteqname = Battery.varToStr({'ctrl', 'controlEquation'});

                eqs([ei.(eieqname), ei.(cteqname)]) = eqs([ei.(cteqname), ei.(eieqname)]);

              otherwise

                error('control type not recognized');

            end


            %% Setup LinearizedProblem that can be processed by MRST Newton API
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

        end

        function state = updateTemperature(model, state)
        % Dispatch the temperature in all the submodels

            elyte   = 'Electrolyte';
            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            co      = 'Coating';
            am      = 'ActiveMaterial';
            cc      = 'CurrentCollector';
            thermal = 'ThermalModel';

            % (here we assume that the ThermalModel has the "parent" grid)
            state.(elyte).T   = state.(thermal).T(model.(elyte).G.mappings.cellmap);

            eldes = {ne, pe};

            for ielde = 1 : numel(eldes)
                elde = eldes{ielde};
                state.(elde).(co).T = state.(thermal).T(model.(elde).(co).G.mappings.cellmap);
                if model.include_current_collectors
                    state.(elde).(cc).T = state.(thermal).T(model.(elde).(cc).G.mappings.cellmap);
                end

                % Update temperature in the active materials of the electrodes.
                state.(elde).(co) = model.(elde).(co).dispatchTemperature(state.(elde).(co));
            end
        end

        function state = updateElectrolyteCoupling(model, state)
        % Assemble the electrolyte coupling by adding the ion sources from the electrodes

            battery = model;

            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            co    = 'Coating';

            F = battery.con.F;

            couplingterms = battery.couplingTerms;

            elyte_e_source = zeros(battery.(elyte).G.cells.num, 1);

            % setup AD
            phi = state.(elyte).phi;
            if isa(phi, 'ADI')
                adsample = getSampleAD(phi);
                adbackend = model.AutoDiffBackend;
                elyte_e_source = adbackend.convertToAD(elyte_e_source, adsample);
            end

            coupnames = model.couplingNames;

            ne_esource = state.(ne).(co).eSource;
            if isa(ne_esource, 'ADI') && ~isa(elyte_e_source, 'ADI')
                adsample = getSampleAD(ne_esource);
                adbackend = model.AutoDiffBackend;
                elyte_e_source = adbackend.convertToAD(elyte_e_source, adsample);
            end
            
            coupterm = getCoupTerm(couplingterms, 'NegativeElectrode-Electrolyte', coupnames);
            elytecells = coupterm.couplingcells(:, 2);
            elyte_e_source(elytecells) = - ne_esource;

            pe_esource = state.(pe).(co).eSource;
            if isa(pe_esource, 'ADI') && ~isa(elyte_e_source, 'ADI')
                adsample = getSampleAD(pe_esource);
                adbackend = model.AutoDiffBackend;
                elyte_e_source = adbackend.convertToAD(elyte_e_source, adsample);
            end

            coupterm = getCoupTerm(couplingterms, 'PositiveElectrode-Electrolyte', coupnames);
            elytecells = coupterm.couplingcells(:, 2);
            elyte_e_source(elytecells) = -pe_esource;

            elyte_c_source = elyte_e_source./(battery.(elyte).sp.z(1)*F);

            state.Electrolyte.eSource    = elyte_e_source;
            state.Electrolyte.massSource = elyte_c_source;

        end


        function state = updateControl(model, state, drivingForces)

            ctrl = "Control";

            switch model.(ctrl).controlPolicy

              case {'CCCV', 'powerControl'}

                % nothing to do here

              case {'CCDischarge', 'CCCharge'}

                E    = state.(ctrl).E;
                I    = state.(ctrl).I;
                time = state.time;

                if model.(ctrl).useCVswitch
                    
                    [ctrlVal, ctrltype] = drivingForces.src(time, value(I), value(E));
                    state.(ctrl).ctrlType = ctrltype;
                    
                else
                    ctrlVal = drivingForces.src(time);
                end
                
                state.(ctrl).ctrlVal  = ctrlVal;

              case 'CC'

                time = state.time;
                ctrlVal = drivingForces.src(time);
                state.(ctrl).ctrlVal  = ctrlVal;

              case 'None'

                % nothing done here. This case is only used for the addVariables method

              otherwise

                error('control type not recognized');

            end


        end

        function state = updateThermalOhmicSourceTerms(model, state)
        % Assemble the ohmic source term :code:`state.jHeatOhmSource`, see :cite:t:`Latz2016`

            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            co      = 'Coating';
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
                    cc_econd = cc_model.effectiveElectronicConductivity;
                    cc_vols  = cc_model.G.cells.volumes;
                    cc_jsq   = computeCellFluxNorm(cc_model, cc_j);
                    state.(elde).(cc).jsq = cc_jsq;  %store square of current density
                    src = subsetPlus(src,cc_vols.*cc_jsq./cc_econd, cc_map);
                    %                    src(cc_map) = src(cc_map) + cc_vols.*cc_jsq./cc_econd;
                end

                co_model = model.(elde).(co);
                co_map   = co_model.G.mappings.cellmap;
                co_j     = state.(elde).(co).jFace;
                co_econd = co_model.effectiveElectronicConductivity;
                co_vols  = co_model.G.cells.volumes;
                co_jsq   = computeCellFluxNorm(co_model, co_j);
                state.(elde).(co).jsq = co_jsq;

                %src(co_map) = src(co_map) + co_vols.*co_jsq./co_econd;
                src = subsetPlus(src, co_vols.*co_jsq./co_econd, co_map);
            end

            % Electrolyte
            elyte_model    = model.(elyte);
            elyte_map      = elyte_model.G.mappings.cellmap;
            elyte_vf       = elyte_model.volumeFraction;
            elyte_j        = state.(elyte).jFace;
            elyte_bruggman = elyte_model.bruggemanCoefficient;
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

            elyte   = 'Electrolyte';
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

        function state = updateThermalIrreversibleReactionSourceTerms(model, state)
        % not used in assembly, just as added variable (in addVariables method). Method updateThermalReactionSourceTerms is used in assembly

            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            co      = 'Coating';
            am      = 'ActiveMaterial';
            sd      = 'SolidDiffusion';
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

                F      = model.(elde).(co).(am).(itf).constants.F;
                n      = model.(elde).(co).(am).(itf).numberOfElectronsTransferred;
                co_map = model.(elde).(co).G.mappings.cellmap;
                vsa    = model.(elde).(co).(am).(itf).volumetricSurfaceArea;
                vols   = model.(elde).(co).G.cells.volumes;

                Rvol = locstate.(elde).(co).(am).(sd).Rvol;
                dUdT = locstate.(elde).(co).(am).(itf).dUdT;
                eta  = locstate.(elde).(co).(am).(itf).eta;

                itf_src = n*F*vols.*Rvol.*eta;

                src(co_map) = src(co_map) + itf_src;


            end

            state.(thermal).jHeatIrrevReactionSource = src;

        end

        function state = updateThermalReversibleReactionSourceTerms(model, state)
        % not used in assembly, just as added variable (in addVariables method). Method updateThermalReactionSourceTerms is used in assembly

            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            co      = 'Coating';
            am      = 'ActiveMaterial';
            itf     = 'Interface';
            sd      = 'SolidDiffusion';
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

                F      = model.(elde).(co).(am).(itf).constants.F;
                n      = model.(elde).(co).(am).(itf).numberOfElectronsTransferred;
                co_map = model.(elde).(co).G.mappings.cellmap;
                vsa    = model.(elde).(co).(am).(itf).volumetricSurfaceArea;
                vols   = model.(elde).(co).G.cells.volumes;

                Rvol = locstate.(elde).(co).(am).(sd).Rvol;
                dUdT = locstate.(elde).(co).(am).(itf).dUdT;
                eta  = locstate.(elde).(co).(am).(itf).eta;

                itf_src = n*F*vols.*Rvol.*T(co_map).*dUdT;

                src(co_map) = src(co_map) + itf_src;

            end

            state.(thermal).jHeatRevReactionSource = src;

        end

        function state = updateThermalReactionSourceTerms(model, state)
        % Assemble the source term from chemical reaction :code:`state.jHeatReactionSource`, see :cite:t:`Latz2016`

            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            co      = 'Coating';
            am      = 'ActiveMaterial';
            itf     = 'Interface';
            sd      = 'SolidDiffusion';
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

                F      = model.(elde).(co).(am).(itf).constants.F;
                n      = model.(elde).(co).(am).(itf).numberOfElectronsTransferred;
                co_map = model.(elde).(co).G.mappings.cellmap;
                vsa    = model.(elde).(co).(am).(itf).volumetricSurfaceArea;
                vols   = model.(elde).(co).G.cells.volumes;

                Rvol = locstate.(elde).(co).(am).(sd).Rvol;
                dUdT = locstate.(elde).(co).(am).(itf).dUdT;
                eta  = locstate.(elde).(co).(am).(itf).eta;

                itf_src = n*F*vols.*Rvol.*(eta + T(co_map).*dUdT);

                src(co_map) = src(co_map) + itf_src;

            end

            state.(thermal).jHeatReactionSource = src;

        end


        function state = updateElectrodeCoupling(model, state)
        % Setup the electrode coupling by updating the potential and concentration of the electrolyte in the active
        % component of the electrodes. There, those quantities are considered as input and used to compute the reaction
        % rate.
        %

            bat = model;
            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            co    = 'Coating';
            am    = 'ActiveMaterial';
            am1   = 'ActiveMaterial1';
            am2   = 'ActiveMaterial2';
            itf   = 'Interface';
            cc    = 'CurrentCollector';

            eldes = {ne, pe};
            phi_elyte = state.(elyte).phi;
            c_elyte = state.(elyte).c;

            elyte_cells = zeros(model.G.cells.num, 1);
            elyte_cells(bat.(elyte).G.mappings.cellmap) = (1 : bat.(elyte).G.cells.num)';


            
            for ind = 1 : numel(eldes)
                elde = eldes{ind};

                switch model.(elde).(co).active_material_type
                  case 'default'
                    ams = {am};
                  case 'composite'
                    ams = {am1, am2};
                  otherwise
                    error('active_material_type not recognized');
                end

                for iam = 1 : numel(ams)
                    
                    amc = ams{iam};
                    state.(elde).(co).(amc).(itf).phiElectrolyte = phi_elyte(elyte_cells(bat.(elde).(co).G.mappings.cellmap));
                    state.(elde).(co).(amc).(itf).cElectrolyte   = c_elyte(elyte_cells(bat.(elde).(co).G.mappings.cellmap));
                    
                end

            end

        end

        function state = setupExternalCouplingNegativeElectrode(model, state)
        %
        % Setup external electronic coupling of the negative electrode at the current collector
        %
            ne = 'NegativeElectrode';
            co = 'Coating';

            if model.(ne).include_current_collectors

                cc = 'CurrentCollector';

                phi   = state.(ne).(cc).phi;
                sigma = state.(ne).(cc).conductivity;

                [jExternal, jFaceExternal] = setupExternalCoupling(model.(ne).(cc), phi, 0, sigma);

                state.(ne).(cc).jExternal = jExternal;
                state.(ne).(cc).jFaceExternal = jFaceExternal;
                state.(ne).(co).jExternal     = 0;
                state.(ne).(co).jFaceExternal = 0;

            else

                phi   = state.(ne).(co).phi;
                sigma = state.(ne).(co).conductivity;

                [jExternal, jFaceExternal] = setupExternalCoupling(model.(ne).(co), phi, 0, sigma);

                state.(ne).(co).jExternal = jExternal;
                state.(ne).(co).jFaceExternal = jFaceExternal;

            end


        end

        function state = setupExternalCouplingPositiveElectrode(model, state)
        %
        % Setup external electronic coupling of the positive electrode at the current collector
        %
            pe   = 'PositiveElectrode';
            ctrl = 'Control';
            co   = 'Coating';

            E   = state.(ctrl).E;

            if model.(pe).include_current_collectors

                cc   = 'CurrentCollector';

                phi   = state.(pe).(cc).phi;
                sigma = state.(pe).(cc).conductivity;

                [jExternal, jFaceExternal] = setupExternalCoupling(model.(pe).(cc), phi, E, sigma);

                state.(pe).(cc).jExternal     = jExternal;
                state.(pe).(cc).jFaceExternal = jFaceExternal;
                state.(pe).(co).jExternal     = 0;
                state.(pe).(co).jFaceExternal = 0;
            else

                phi   = state.(pe).(co).phi;
                sigma = state.(pe).(co).conductivity;

                [jExternal, jFaceExternal] = setupExternalCoupling(model.(pe).(co), phi, E, sigma);

                state.(pe).(co).jExternal     = jExternal;
                state.(pe).(co).jFaceExternal = jFaceExternal;

            end

        end

        function state = setupEIEquation(model, state)

            pe   = 'PositiveElectrode';
            ctrl = 'Control';

            I = state.(ctrl).I;
            E = state.(ctrl).E;

            if model.include_current_collectors

                mat   = 'CurrentCollector';

            else

                mat   = 'Coating';

            end

            phi = state.(pe).(mat).phi;

            coupterm = model.(pe).(mat).externalCouplingTerm;
            faces    = coupterm.couplingfaces;
            cond_pcc = model.(pe).(mat).effectiveElectronicConductivity;
            [trans_pcc, cells] = model.(pe).(mat).operators.transFaceBC(faces);

            state.Control.EIequation = sum(cond_pcc.*trans_pcc.*(phi(cells) - E)) - I;


        end

        function primaryvarnames = getPrimaryVariableNames(model)

            primaryvarnames = model.primaryVarNames;

        end

        function forces = getValidDrivingForces(model)

            forces = getValidDrivingForces@PhysicalModel(model);

            ctrl = 'Control';
            switch model.(ctrl).controlPolicy
              case 'CCCV'
                forces.CCCV = true;
              case 'CCDischarge'
                forces.CCDischarge = true;
                forces.src = [];
              case 'CCCharge'
                forces.CCCharge = true;
                forces.src = [];
              case 'powerControl'
                forces.powerControl = true;
                forces.src = [];
              case 'CC'
                forces.CC = true;
                forces.src = [];
              case 'None'
                % used only in addVariables
              otherwise
                error('Error controlPolicy not recognized');
            end
            % TODO this is a hack to get thing go
            forces.Imax = [];

        end

        function model = validateModel(model, varargin)

            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            co    = 'Coating';
            cc    = 'CurrentCollector';

            model.(elyte).AutoDiffBackend   = model.AutoDiffBackend;
            model.(elyte)                   = model.(elyte).validateModel(varargin{:});

            eldes = {ne, pe};

            for ind = 1 : numel(eldes)

                elde = eldes{ind};

                model.(elde).(co).AutoDiffBackend = model.AutoDiffBackend;
                model.(elde).(co)                 = model.(elde).(co).validateModel(varargin{:});

                if model.(elde).include_current_collectors

                    model.(elde).(cc).AutoDiffBackend = model.AutoDiffBackend;
                    model.(elde).(cc)                 = model.(elde).(cc).validateModel(varargin{:});

                end

            end

            if isempty(model.computationalGraph)
                model = model.setupComputationalGraph();
            end

            cgt = model.computationalGraph;
            

        end


        function [state, report] = updateState(model, state, problem, dx, drivingForces)

            [state, report] = updateState@BaseModel(model, state, problem, dx, drivingForces);

            %% cap concentrations
            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            co    = 'Coating';
            am    = 'ActiveMaterial';
            am1   = 'ActiveMaterial1';
            am2   = 'ActiveMaterial2';
            sd    = 'SolidDiffusion';
            itf   = 'Interface';

            cmin = model.cmin;

            state.(elyte).c = max(cmin, state.(elyte).c);

            eldes = {ne, pe};
            for ind = 1 : numel(eldes)
                elde = eldes{ind};

                switch model.(elde).(co).active_material_type
                  case 'default'
                    ams = {am};
                  case 'composite'
                    ams = {am1, am2};
                  otherwise
                    error('active_material_type not recognized');
                end
                
                for iam = 1 : numel(ams)
                    
                    amc = ams{iam};
                    cmax = model.(elde).(co).(amc).(itf).saturationConcentration;
                    switch model.(elde).(co).(amc).diffusionModelType
                      case 'simple'
                        state.(elde).(co).(amc).(sd).cAverage = max(cmin, state.(elde).(co).(amc).(sd).cAverage);
                        state.(elde).(co).(amc).(sd).cAverage = min(cmax, state.(elde).(co).(amc).(sd).cAverage);
                      case 'full'
                        state.(elde).(co).(amc).(sd).c = max(cmin, state.(elde).(co).(amc).(sd).c);
                        state.(elde).(co).(amc).(sd).c = min(cmax, state.(elde).(co).(amc).(sd).c);
                      otherwise
                        error('diffusionModelType not recognized')

                    end

                end
            end

            ctrl = 'Control';
            state.(ctrl) = model.(ctrl).updateControlState(state.(ctrl));

            report = [];

        end

        function cleanState = addStaticVariables(model, cleanState, state)
        % Variables that are no AD initiated (but should be "carried over")

            cleanState = addStaticVariables@BaseModel(model, cleanState, state);


            cleanState.time = state.time;

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

            if strcmp(model.(ctrl).controlPolicy, 'powerControl')
                state.(ctrl).time = state.time;
            end

            state.(ctrl) = model.(ctrl).prepareStepControl(state.(ctrl), state0.(ctrl), dt, drivingForces);

        end

        function model = setupCapping(model)

            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            co      = 'Coating';
            itf     = 'Interface';

            eldes = {pe, ne};

            for ielde = 1 : numel(eldes)
                
                elde = eldes{ielde};

                switch model.(elde).(co).active_material_type

                  case 'default'

                    am = 'ActiveMaterial';
                    
                    cmaxs{ielde} = model.(elde).(co).(am).(itf).saturationConcentration;
                    
                  case 'composite'
                    
                    am1 = 'ActiveMaterial1';
                    am2 = 'ActiveMaterial2';

                    ams = {am1, am2};

                    cmax = 0;
                    
                    for iam = 1 : numel(ams)
                        
                        amc = ams{iam};
                        cmax = max(cmax, model.(elde).(co).(amc).(itf).saturationConcentration);
                        
                    end
                    
                    cmaxs{ielde} = cmax;
                    
                  otherwise
                    
                    error('active_material_type not recognized');
                    
                end
                

            end
            
            model.cmin = 1e-5*max(cmaxs{1}, cmaxs{2});
            

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

        function [found, varind] = getVarIndex(varname, pvarnames)

            varname   = Battery.varToStr(varname);
            pvarnames = cellfun(@(name) Battery.varToStr(name), pvarnames, 'uniformoutput', false);

            [found, varind] = ismember(varname, pvarnames);

        end

        function str = varToStr(varname)

            str = strjoin(varname, '_');

        end

        function str = shortenName(name)

            namemapping = {'NegativeElectrode', 'ne'      ; ...
                           'PositiveElectrode', 'pe'      ; ...
                           'Coating'          , 'co'      ; ...
                           'ActiveMaterial'   , 'am'      ; ...
                           'ActiveMaterial1'  , 'am1'     ; ...
                           'ActiveMaterial2'  , 'am2'     ; ...
                           'CurrentCollector' , 'cc'      ; ...
                           'Electrolyte'      , 'elyte'   ; ...
                           'Interface'        , 'itf'     ; ...
                           'SolidDiffusion'   , 'sd'      ; ...
                           'ThermalModel'     , 'thermal' ; ...
                           'Control'          , 'ctrl'};


            [found, ind] = ismember(name, namemapping(:, 1));

            if found
                str = namemapping{ind, 2};
            else
                str = name;
            end

        end

    end

end



%{
Copyright 2021-2023 SINTEF Industry, Sustainable Energy Technology
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
