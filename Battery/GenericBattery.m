classdef GenericBattery < BaseModel
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
        include_swelling % true if swelling is included
        
        % We write here the jsonstruct that has been process when constructing inputparams
        jsonstruct

    end

    methods

        function model = GenericBattery(inputparams)

            model = model@BaseModel();

            % All the submodels should have same backend (this is not assigned automaticallly for the moment)
            model.AutoDiffBackend = SparseAutoDiffBackend('useBlocks', true);

            %% Setup the model using the input parameters
            fdnames = {'G'                         , ...
                       'couplingTerms'             , ...
                       'initT'                     , ...
                       'use_thermal'               , ...
                       'include_current_collectors', ...
                       'use_thermal'               , ...
                       'SOC'};

            model = dispatchParams(model, inputparams, fdnames);

            model.jsonstruct = inputparams.buildJsonStruct();

            model.NegativeElectrode = Electrode(inputparams.NegativeElectrode);
            model.PositiveElectrode = Electrode(inputparams.PositiveElectrode);
            model.Separator         = Separator(inputparams.Separator);

            % We setup the electrolyte model (in particular we compute the volume fraction from the other components)
            model = model.setupElectrolyteModel(inputparams);

            if model.use_thermal
                model.ThermalModel = ThermalComponent(inputparams.ThermalModel);
            end

            model.Control = model.setupControl(inputparams.Control);


            if model.use_thermal
                % setup Thermal Model by assigning the effective heat capacity and conductivity, which is computed from the sub-models.
                model = model.setupThermalModel();
            end

            % setup couplingNames
            model.couplingNames = cellfun(@(x) x.name, model.couplingTerms, 'uniformoutput', false);

            % setup some mappings (mappings from electrodes to electrolyte)
            model = model.setupMappings();

            % setup capping
            model = model.setupCapping();

            % setup scalings
            model = model.setupScalings();

            % equip model for computation
            %
            % Define short-names
            %
            shortNames =  {'Electrolyte'            ,'elyte'  ; ...
                           'NegativeElectrode'      ,'ne'     ; ...
                           'PositiveElectrode'      ,'pe'     ; ...
                           'Coating'                ,'co'     ; ...
                           'ActiveMaterial'         ,'am'     ; ...
                           'ActiveMaterial1'        ,'am1'    ; ...
                           'ActiveMaterial2'        ,'am2'    ; ...
                           'CurrentCollector'       ,'cc'     ; ...
                           'Interface'              ,'itf'    ; ...
                           'SolidDiffusion'         ,'sd'     ; ...
                           'ThermalModel'           ,'thermal'; ...
                           'Control'                ,'ctrl'   ; ...
                           'SideReaction'           ,'sr'     ; ...
                           'SolidElectrodeInterface','sei'};

            model = model.equipModelForComputation('shortNames', shortNames);

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
            sr      = 'SideReaction';
            sei     = 'SolidElectrodeInterface';

            if ~model.use_thermal
                % we register the temperature variable, as it is not done by the ThermalModel which is empty in this case
                model = model.registerVarName({thermal, 'T'});
            end

            eldes = {ne, pe};


            %% Temperature dispatch functions
            fn = @GenericBattery.updateTemperature;

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

                if  model.(elde).(co).activeMaterialModelSetup.composite
                    ams = {am1, am2};
                else
                    ams = {am};
                end

                for iam = 1 : numel(ams)

                    amc = ams{iam};

                    fn = @GenericBattery.updateElectrodeCoupling;
                    inputnames = {{elyte, 'c'}, ...
                                  {elyte, 'phi'}};

                    model = model.registerPropFunction({{elde, co, amc, itf, 'phiElectrolyte'}, fn, inputnames});
                    model = model.registerPropFunction({{elde, co, amc, itf, 'cElectrolyte'}  , fn, inputnames});

                    switch model.(elde).(co).activeMaterialModelSetup.SEImodel
                      case {'none', 'Bolay'}
                        % nothing more to add
                      case 'Safari'
                        model = model.registerPropFunction({{elde, co, amc, sr, 'phiElectrolyte'}, fn, inputnames});
                      otherwise
                        error('SEI model not recognized');
                    end

                end

            end

            fn = @GenericBattery.updateElectrolyteCoupling;
            inputnames = {{ne, co, 'eSource'}, ...
                          {pe, co, 'eSource'}};
            model = model.registerPropFunction({{elyte, 'massSource'}, fn, inputnames});
            model = model.registerPropFunction({{elyte, 'eSource'}, fn, inputnames});

            % Function that assemble the control equation
            fn = @GenericBattery.setupEIEquation;
            inputnames = {{ctrl, 'E'}, ...
                          {ctrl, 'I'}};
            if model.include_current_collectors
                inputnames{end + 1} = {pe, cc, 'phi'};
            else
                inputnames{end + 1} = {pe, co, 'phi'};
            end
            model = model.registerPropFunction({{ctrl, 'EIequation'}, fn, inputnames});

            if model.include_current_collectors && model.(pe).use_normed_current_collector
                fn = @GenericBattery.updateCurrentCollectorPhiRef;
                inputnames = {{ctrl, 'E'}};
                model = model.registerPropFunction({{pe, cc, 'phiRef'}, fn, inputnames});
            end

            inputnames = {};
            fn = @GenericBattery.updateControl;
            fn = {fn, @(propfunction) PropFunction.drivingForceFuncCallSetupFn(propfunction)};
            switch model.(ctrl).controlPolicy
              case {'CCDischarge', 'CCCharge', 'CC', 'timeControl'}
                model = model.registerPropFunction({{ctrl, 'ctrlVal'}, fn, inputnames});
              case {'CCCV', 'Impedance'}
                % do nothing
              otherwise
                error('controlPolicy not recognized');
            end
            model = model.registerPropFunction({{ctrl, 'ctrlType'}, fn, inputnames});

            if model.include_swelling

                fn = @GenericBattery.updateSwellingElectrolyteAccumTerm;
                fn = {fn, @(propfunction) PropFunction.accumFuncCallSetupFn(propfunction)};
                inputNames = {{elyte,'c'}, {elyte, 'volumeFraction'}};
                model = model.registerPropFunction({{elyte, 'massAccum'}, fn, inputNames});

                fn = @GenericBattery.updateSwellingElectrolyteVolumeFraction;
                inputNames = {};
                for ielde = 1 : numel(eldes)
                    elde = eldes{ielde};
                    if model.(elde).coatingModelSetup.swelling
                        inputNames = horzcat(inputNames, {{elde, co, 'volumeFraction'}});                        
                    end
                end                
                model = model.registerPropFunction({{elyte,'volumeFraction'}, fn, inputNames});

                fn = @GenericBattery.updateSwellingElectrolyteConvFlux;
                inputNames = {{elyte, 'j'}, {elyte, 'c'}};
                for ielde = 1 : numel(eldes)
                    elde = eldes{ielde};
                    if model.(elde).coatingModelSetup.swelling
                        inputNames = horzcat(inputNames, {{elde, co, am, sd, 'x'}});
                    end
                end                
                model = model.registerPropFunction({{elyte,'convFlux'}, fn, inputNames});
                
            end
            
            %% Function that update the Thermal Ohmic Terms

            if model.use_thermal

                fn = @GenericBattery.updateThermalOhmicSourceTerms;
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
                fn = @GenericBattery.updateThermalChemicalSourceTerms;
                inputnames = {{elyte, 'diffFlux'}, ...
                              {elyte, 'D'}       , ...
                              VarName({elyte}, 'dmudcs', 2)};
                model = model.registerPropFunction({{thermal, 'jHeatChemicalSource'}, fn, inputnames});

                %% Function that updates Thermal Reaction Terms
                fn = @GenericBattery.updateThermalReactionSourceTerms;
                inputnames = {{thermal, 'T'}           , ...
                              {ne, co, am, sd, 'Rvol'} , ...
                              {ne, co, am, itf, 'eta'} , ...
                              {pe, co, am, sd, 'Rvol'} , ...
                              {pe, co, am, itf, 'eta'}};
                for ielde = 1 : numel(eldes)
                    elde = eldes{ielde};
                    if model.(elde).(co).(am).(itf).includeEntropyChange
                        inputnames{end + 1} = {elde, co, am, itf, 'dUdT'};
                    end
                end
                model = model.registerPropFunction({{thermal, 'jHeatReactionSource'}, fn, inputnames});

            else
                model = model.removeVarName({elyte, 'diffFlux'});
            end

            %% Functions that setup external  coupling for negative electrode


            fns{1} = @GenericBattery.setupExternalCouplingNegativeElectrode;
            fns{2} = @GenericBattery.setupExternalCouplingPositiveElectrode;

            for ielde = 1 : numel(eldes)

                elde = eldes{ielde};
                fn   = fns{ielde};

                if model.include_current_collectors

                    inputnames = {{elde, cc, 'phi'}         , ...
                                  {elde, cc, 'conductivity'}, ...
                                  {elde, co, 'phi'}         , ...
                                  {elde, co, 'conductivity'}};

                    model = model.registerPropFunction({{elde, cc, 'jExternal'}, fn, inputnames});
                    model = model.registerPropFunction({{elde, co, 'jExternal'}, fn, inputnames});

                    if model.use_thermal
                        model = model.registerPropFunction({{elde, cc, 'jFaceExternal'}, fn, inputnames});
                    end


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
            model = model.setAsStaticVarNames(varnames);


        end

        function model = setupScalings(model)

            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            co      = 'Coating';
            cc      = 'CurrentCollector';
            am      = 'ActiveMaterial';
            am1     = 'ActiveMaterial1';
            am2     = 'ActiveMaterial2';
            elyte   = 'Electrolyte';
            am      = 'ActiveMaterial';
            itf     = 'Interface';
            sd      = 'SolidDiffusion';
            thermal = 'ThermalModel';
            ctrl    = 'Control';
            sr      = 'SideReaction';
            sei     = 'SolidElectrodeInterface';

            eldes = {ne, pe};

            % We compute a scaling for the reaction rate coefficient j0 (see Interface model) and Volumetric reaction rate Rvol
            % in mol/(s*m^3) (see SolidDiffusion model)

            for ielde = 1 : numel(eldes)

                elde = eldes{ielde};

                if  model.(elde).(co).activeMaterialModelSetup.composite
                    ams = {am1, am2};
                else
                    ams = {am};
                end

                % for simplicity, we use only the first material in case of composite;

                amc = ams{1};

                state.T            = model.initT;
                state.cElectrolyte = model.(elyte).species.nominalConcentration;
                % We use half max concentration
                state.cElectrodeSurface = 0.5*model.(elde).(co).(amc).(itf).saturationConcentration;

                state = model.(elde).(co).(amc).(itf).updateReactionRateCoefficient(state);

                j0s(ielde) = state.j0;

                vsa = model.(elde).(co).(amc).(sd).volumetricSurfaceArea;
                F   = model.con.F;

                Rvols(ielde) = j0s(ielde)*vsa/F;

                % reference for the conductivity
                condRef.(elde).(co) = model.(elde).(co).effectiveElectronicConductivity;
                if model.(elde).include_current_collectors
                    condRef.(elde).(cc) = model.(elde).(cc).effectiveElectronicConductivity;
                end

            end

            j0Ref   = mean(j0s);
            RvolRef = mean(Rvols);

            volRef.(ne).(co) = mean(model.(ne).(co).G.getVolumes());
            volRef.(pe).(co) = mean(model.(pe).(co).G.getVolumes());

            if model.include_current_collectors
                volRef.(ne).(cc) = mean(model.(ne).(cc).G.getVolumes());
                volRef.(pe).(cc) = mean(model.(pe).(cc).G.getVolumes());
            end

            volRef.(elyte) = mean(model.(elyte).G.getVolumes());

            F = model.con.F;

            scalings = {};

            scalings{end + 1} = {{elyte, 'chargeCons'}, F*volRef.(elyte)*RvolRef};
            scalings{end + 1} = {{elyte, 'massCons'}, volRef.(elyte)*RvolRef};
            
            for ielde = 1 : numel(eldes)

                elde = eldes{ielde};
                
                scalings{end + 1} = {{elde, co, 'chargeCons'}, F*volRef.(elde).(co)*RvolRef};

                if model.include_current_collectors
                    % We use the same scaling as for the coating multiplied by the conductivity ration
                    coef = condRef.(elde).(cc)/condRef.(elde).(co);
                    scalings{end + 1} = {{elde, cc, 'chargeCons'}, coef*F*volRef.(elde).(cc)*RvolRef};
                end

                if  model.(elde).(co).activeMaterialModelSetup.composite
                    ams = {am1, am2};
                else
                    ams = {am};
                end

                for iam = 1 : numel(ams)

                    amc = ams{iam};

                    switch model.(elde).(co).(amc).diffusionModelType

                      case 'simple'


                        scalings{end + 1} = {{elde, co, amc, sd, 'massCons'}, RvolRef};

                        rp  = model.(elde).(co).(amc).(sd).particleRadius;
                        vsa = model.(elde).(co).(amc).(sd).volumetricSurfaceArea;
                        D0  = model.(elde).(co).(amc).(sd).referenceDiffusionCoefficient;
                        scalings{end + 1} = {{elde, co, amc, sd, 'solidDiffusionEq'}, rp*RvolRef/(5*vsa*D0)};

                      case {'full', 'swelling'}

                        rp   = model.(elde).(co).(amc).(sd).particleRadius;

                        volp = 4/3*pi*rp^3;

                        coef = RvolRef*volp;
                        scalings{end + 1} = {{elde, co, amc, sd, 'massCons'}, coef};
                        scalings{end + 1} = {{elde, co, amc, sd, 'solidDiffusionEq'}, coef};

                      otherwise

                        error('diffusionModelType not recognized');

                    end


                    switch model.(elde).(co).activeMaterialModelSetup.SEImodel

                      case 'none'

                        % nothing more to add

                      case 'Bolay'

                        F   = model.(elde).(co).(amc).(itf).constants.F;
                        vsa = model.(elde).(co).(amc).(sd).volumetricSurfaceArea;
                        L   = model.(elde).(co).(amc).(itf).SEIlengthInitial;
                        k   = model.(elde).(co).(amc).(itf).SEIionicConductivity;

                        SEIvoltageDropRef = F*RvolRef/vsa*L/k;
                        scalings{end + 1} = {{elde, co, amc, itf, 'SEIvoltageDropEquation'}, SEIvoltageDropRef};
                        model.(elde).(co).(amc).(itf).SEIvoltageDropRef = SEIvoltageDropRef;

                        De = model.(elde).(co).(amc).(itf).SEIelectronicDiffusionCoefficient;
                        ce = model.(elde).(co).(amc).(itf).SEIintersticialConcentration;

                        scalings{end + 1} = {{elde, co, amc, itf, 'SEImassCons'}, De*ce/L};

                      case 'Safari'
                        % we use the scaling given for {elde, co, amc, sd, 'massCons'} and use the same ratio as in the
                        % scaling for the standalone model SEIActiveMaterial for the SEI parts.

                        rp  = model.(elde).(co).(am).(sd).particleRadius;
                        vsa = model.(elde).(co).(am).(sd).volumetricSurfaceArea;
                        Mw  = model.(elde).(co).(am).(sei).molecularWeight;
                        rho = model.(elde).(co).(am).(sei).density;

                        scalings{end + 1} = {{elde, co, am, sei, 'widthEq'}, RvolRef/vsa*Mw/rho};

                        deltaref = 1*nano*meter;

                        scalings{end + 1} = {{elde, co, am, sei, 'interfaceBoundaryEq'}, deltaref*RvolRef/vsa};
                        scalings{end + 1} = {{elde, co, am, sei, 'massCons'}, deltaref*RvolRef/vsa};

                      otherwise

                        error('SEI model not recognized');
                    end

                end

            end

            model.scalings = scalings;

        end

        function printJsonStruct(model)

            fjv = flattenJsonStruct(model.jsonstruct);
            fjv.print();

        end


        function control = setupControl(model, inputparams)

            C = computeCellCapacity(model);

            switch inputparams.controlPolicy

              case 'timeControl'

                control = TimeControlModel(inputparams);

              case "Impedance"

                control = ImpedanceControlModel(inputparams);

              case "CCDischarge"

                control = CCDischargeControlModel(inputparams);
                rate = control.DRate;
                control.Imax = (C/hour)*rate;

              case 'CCCharge'

                control = CCChargeControlModel(inputparams);
                if isempty(control.Imax)
                    rate = control.CRate;
                    control.Imax = (C/hour)*rate;
                end

              case "CCCV"

                control = CcCvControlModel(inputparams);
                CRate = control.CRate;
                DRate = control.DRate;
                control.ImaxCharge    = (C/hour)*CRate;
                control.ImaxDischarge = (C/hour)*DRate;

              case "powerControl"

                control = PowerControlModel(inputparams);

              case "CC"

                control = CCcontrolModel(inputparams);

              otherwise

                error('Error controlPolicy not recognized');
            end

        end

        function model = setupThermalModel(model, inputparams)
        % Setup the thermal model :attr:`ThermalModel`. Here, :code:`inputparams` is instance of
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
            nc = G.getNumberOfCells();

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
            elytecelltbl.cells = (1 : G_elyte.getNumberOfCells())';
            elytecelltbl.globalcells = G_elyte.mappings.cellmap;
            elytecelltbl = IndexArray(elytecelltbl);

            eldes = {ne, pe};

            for ind = 1 : numel(eldes)

                elde = eldes{ind};
                G_elde  = model.(elde).(co).G;
                clear eldecelltbl;
                eldecelltbl.cells = (1 : G_elde.getNumberOfCells())';
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

        function model = setupElectrolyteModel(model, inputparams)
        % Assign the electrolyte volume fractions in the different regions

            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            co    = 'Coating';
            sep   = 'Separator';

            eldes = {ne, pe};
            
            elyte_cells = zeros(model.G.getNumberOfCells(), 1);
            elyte_cells(inputparams.(elyte).G.mappings.cellmap) = (1 : inputparams.(elyte).G.getNumberOfCells())';

            inputparams.(elyte).volumeFraction = ones(inputparams.(elyte).G.getNumberOfCells(), 1);
            inputparams.(elyte).volumeFraction = subsasgnAD(inputparams.(elyte).volumeFraction, elyte_cells(model.(ne).(co).G.mappings.cellmap), 1 - model.(ne).(co).volumeFraction);
            inputparams.(elyte).volumeFraction = subsasgnAD(inputparams.(elyte).volumeFraction, elyte_cells(model.(pe).(co).G.mappings.cellmap), 1 - model.(pe).(co).volumeFraction);
            inputparams.(elyte).volumeFraction = subsasgnAD(inputparams.(elyte).volumeFraction, elyte_cells(model.(sep).G.mappings.cellmap), model.(sep).porosity);

            if inputparams.(ne).coatingModelSetup.swelling || inputparams.(pe).coatingModelSetup.swelling
                model.include_swelling = true;
                model.(elyte) = ElectrolyteSwelling(inputparams.(elyte));
            else
                model.include_swelling = false;
                model.(elyte) = Electrolyte(inputparams.(elyte));
            end

        end

        function initstate = setupInitialState(model, jsonstruct)
        % Setup the values of the primary variables at initial state
        %
        % The jsonstruct structure (standard matlab struct type) contains some parameters for the initialization. For
        % the moment, it includes only the initial electrolyte concentration
        %

            if nargin < 2 | isempty(jsonstruct)
                jsonstruct.Electrolyte.initialConcentration = model.Electrolyte.species.nominalConcentration;
                if ~isempty(model.Electrolyte.nominalEthyleneCarbonateConcentration)
                    jsonstruct.Electrolyte.initialEthyleneCarbonateConcentration = model.Electrolyte.nominalEthyleneCarbonateConcentration;
                end
            end

            nc = model.G.getNumberOfCells();

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
            sr      = 'SideReaction';
            sei     = 'SolidElectrodeInterface';

            initstate.(thermal).T = T*ones(nc, 1);

            %% Synchronize temperatures
            initstate = model.updateTemperature(initstate);

            %% Setup initial state for electrodes

            eldes = {ne, pe};

            for ielde = 1 : numel(eldes)

                elde = eldes{ielde};

                if  model.(elde).(co).activeMaterialModelSetup.composite
                    ams = {am1, am2};
                else
                    ams = {am};
                end

                for iam = 1 : numel(ams)

                    amc = ams{iam};

                    elde_itf = bat.(elde).(co).(amc).(itf);

                    theta = SOC*(elde_itf.guestStoichiometry100 - elde_itf.guestStoichiometry0) + elde_itf.guestStoichiometry0;
                    c     = theta*elde_itf.saturationConcentration;
                    nc    = model.(elde).(co).G.getNumberOfCells();

                    switch model.(elde).(co).(amc).diffusionModelType
                      case 'simple'
                        initstate.(elde).(co).(amc).(sd).cSurface = c*ones(nc, 1);
                        initstate.(elde).(co).(amc).(sd).cAverage = c*ones(nc, 1);
                      case {'full'}
                        initstate.(elde).(co).(amc).(sd).cSurface = c*ones(nc, 1);
                        N = model.(elde).(co).(amc).(sd).N;
                        np = model.(elde).(co).(amc).(sd).np; % Note : we have by construction np = nc
                        initstate.(elde).(co).(amc).(sd).c = c*ones(N*np, 1);
                      case  'swelling'
                        % theta is interpretated as a fill-in level
                        
                        compmodel = model.(elde).(co);
                        compmodel = compmodel.registerVarAndPropfuncNames();
                        compmodel = compmodel.removePropFunction({amc, sd, 'x'});
                        compmodel = compmodel.setupComputationalGraph();

                        clear state
                        state.(amc).(sd).x = theta;
                        state = compmodel.evalVarName(state, {'volumeFraction'});
                        state = compmodel.evalVarName(state, {amc, sd, 'radiusElongation'});
                        
                        re = state.(am).(sd).radiusElongation;
                        vf = state.volumeFraction;

                        N  = model.(elde).(co).(amc).(sd).N;
                        np = model.(elde).(co).(amc).(sd).np; % Note : we have by construction np = nc
                        
                        c = model.(elde).(co).maximumTotalConcentration/vf*theta;

                        initstate.(elde).(co).volumeFraction              = vf .* ones(nc, 1);
                        initstate.(elde).(co).(amc).(sd).cSurface         = c*ones(nc, 1);
                        initstate.(elde).(co).(amc).(sd).c                = c*ones(N*np, 1);
                        initstate.(elde).(co).(amc).(sd).radiusElongation = re*ones(nc, 1);
                        
                      otherwise
                        
                        error('diffusionModelType not recognized')
                        
                    end

                    if elde_itf.useDoubleLayerCapacity
                        initstate.(elde).(co).(amc).(itf).capacityR = zeros(nc, 1);
                    end

                    switch model.(elde).(co).activeMaterialModelSetup.SEImodel

                      case 'none'
                        %  nothing to do

                      case 'Safari'

                        N = model.(elde).(co).(amc).(sei).N;
                        np = model.(elde).(co).(amc).(sei).np; % Note : we have by construction np = nc

                        cExternal = jsonstruct.Electrolyte.initialEthyleneCarbonateConcentration;

                        initstate.(elde).(co).(amc).(sei).cExternal   = cExternal;
                        initstate.(elde).(co).(amc).(sei).c           = cExternal*ones(N*np, 1);
                        initstate.(elde).(co).(amc).(sei).cInterface  = cExternal*ones(np, 1);
                        initstate.(elde).(co).(amc).(sei).delta       = 5*nano*meter*ones(np, 1);
                        initstate.(elde).(co).(amc).totalReactionRate = zeros(np, 1);

                      case 'Bolay'

                        l = model.(elde).(co).(amc).(itf).SEIlengthInitial/model.(elde).(co).(amc).(itf).SEIlengthRef;
                        initstate.(elde).(co).(amc).(itf).normalizedSEIlength      = l*ones(np, 1);
                        initstate.(elde).(co).(amc).(itf).normalizedSEIvoltageDrop = zeros(np, 1);

                        initstate = model.evalVarName(initstate, {ne, co, amc, itf, 'SEIlength'});

                      otherwise

                        error('SEI model not recognized');

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

            initstate.(elyte).phi = zeros(bat.(elyte).G.getNumberOfCells(), 1) - ref;
            initstate.(elyte).c   = jsonstruct.Electrolyte.initialConcentration*ones(bat.(elyte).G.getNumberOfCells(), 1);


            %% Setup initial variables needed in case of double layer

            for ielde = 1 : numel(eldes)

                elde = eldes{ielde};

                if  model.(elde).(co).activeMaterialModelSetup.composite
                    ams = {am1, am2};
                else
                    ams = {am};
                end

                for iam = 1 : numel(ams)

                    amc = ams{iam};

                    if model.(elde).(co).(amc).(itf).useDoubleLayerCapacity
                        initstate = model.evalVarName(initstate, {elde, co, amc, itf, 'cElectrolyte'});
                        initstate = model.evalVarName(initstate, {elde, co, amc, itf, 'phiElectrode'});
                        initstate = model.evalVarName(initstate, {elde, co, amc, itf, 'phiElectrolyte'});
                    end

                end

            end

            %% Setup initial Current collectors state

            if model.(ne).include_current_collectors
                OCP = initstate.(ne).(co).(amc).(itf).OCP;
                OCP = OCP(1) .* ones(bat.(ne).(cc).G.getNumberOfCells(), 1);
                initstate.(ne).(cc).phi = OCP - ref;
            end

            if model.(pe).include_current_collectors
                OCP = initstate.(pe).(co).(amc).(itf).OCP;
                OCP = OCP(1) .* ones(bat.(pe).(cc).G.getNumberOfCells(), 1);
                if model.(pe).use_normed_current_collector
                    initstate.(pe).(cc).scaledDeltaPhi = 0*OCP;
                else
                    initstate.(pe).(cc).phi = OCP - ref;
                end
            end

            initstate.(ctrl).E = OCP(1) - ref;

            switch model.(ctrl).controlPolicy

              case {'Impedance'}

                initstate.(ctrl).I = 0;

              case {'timeControl'}

                %  We initiate to some values, but they should be overriden as the simulation starts
                initstate.(ctrl).I        = 0;
                initstate.(ctrl).ctrlType = 'constantCurrent';

              case {'CCDischarge', 'CCCharge'}

                initstate.(ctrl).ctrlType = 'constantCurrent';
                initstate.(ctrl).I = model.(ctrl).Imax;

              case 'CC'

                initstate.(ctrl).ctrlType = 'constantCurrent';
                initstate.(ctrl).I = 0;

              case 'CCCV'

                initstate.(ctrl).numberOfCycles = 0;

                switch model.(ctrl).initialControl
                  case 'discharging'
                    initstate.(ctrl).ctrlType = 'CC_discharge1';
                    initstate.(ctrl).I        = model.(ctrl).ImaxDischarge;
                  case 'charging'
                    initstate.(ctrl).ctrlType = 'CC_charge1';
                    initstate.(ctrl).I        = - model.(ctrl).ImaxCharge;
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
            inputparams = ControlModelInputParams([]);
            model.Control = ControlModel(inputparams);
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
            sd      = 'SolidDiffusion';
            thermal = 'ThermalModel';
            ctrl    = 'Control';

            state.(elyte).massAccum = [];
            state.(elyte).massCons = [];

            eldes = {ne, pe};

            for ielde = 1 : numel(eldes)

                elde = eldes{ielde};

                if  model.(elde).(co).activeMaterialModelSetup.composite
                    ams = {am1, am2};
                else
                    ams = {am};
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

        function jsonstruct = exportParams(model, varargin)

            opt = struct('removeEmpty', false);
            opt = merge_options(opt, varargin{:});

            jsonstruct = exportParams@BaseModel(model);

            fdnames = {'SOC'        , ...
                       'initT'      , ...
                       'use_thermal', ...
                       'include_current_collectors'};

            for ifd = 1 : numel(fdnames)
                fdname = fdnames{ifd};
                jsonstruct.(fdname) = model.(fdname);
            end

            jsonstruct = removeJsonStructEmptyField(jsonstruct);

        end

        function state = updateCurrentCollectorPhiRef(model, state);

            pe   = 'PositiveElectrode';
            cc   = 'CurrentCollector';
            ctrl = 'Control';

            state.(pe).(cc).phiRef = state.(ctrl).E;

        end

        %% Update at each step the electrolyte volume fractions in the different regions (neg_elde, elyte, pos_elde)
        function state = updateSwellingElectrolyteVolumeFraction(model, state)

            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            co    = 'Coating';
            am    = 'ActiveMaterial';
            sep   = 'Separator';

            elyte_cells = zeros(model.G.getNumberOfCells(), 1);
            elyte_cells(model.(elyte).G.mappings.cellmap) = (1 : model.(elyte).G.getNumberOfCells)';
            
            % Define the volumeFraction in the electrodes
            eldes = {ne, pe};
            % Initialisation of AD for the porosity of the elyte
            state.(elyte).volumeFraction = 0 * state.(elyte).c;

            % Define the porosity in the separator
            sep_cells = elyte_cells(model.(sep).G.mappings.cellmap);
            state.(elyte).volumeFraction = subsasgnAD(state.(elyte).volumeFraction, sep_cells, model.(sep).porosity);
            
            for ielde = 1 : numel(eldes)
                elde = eldes{ielde};
                if model.(elde).coatingModelSetup.swelling
                    porosity = 1 - state.(elde).(co).volumeFraction;
                else
                    porosity = 1 - model.(elde).(co).volumeFraction;
                end
                state.(elyte).volumeFraction = subsasgnAD(state.(elyte).volumeFraction                     , ...
                                                          elyte_cells(model.(elde).(co).G.mappings.cellmap), ...
                                                          porosity);                    
            end

            if model.use_thermal
                warning('swelling model not in sync with thermal simulation yet...')
                model.(elyte).EffectiveThermalConductivity = model.(elyte).volumeFraction.*model.(elyte).thermalConductivity;
            end

        end

        
        %% Definition of the accumulation term (dc/dt)
        function state = updateSwellingElectrolyteAccumTerm(model, state, state0, dt)

            elyte   = 'Electrolyte';
            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            co      = 'Coating';
            am      = 'ActiveMaterial';

            eldes = {ne, pe};
            
            c         = state.(elyte).c;
            vf        = state.(elyte).volumeFraction;
            c0        = state0.(elyte).c;

            elyte_cells = zeros(model.G.getNumberOfCells(), 1);
            elyte_cells(model.(elyte).G.mappings.cellmap) = (1 : model.(elyte).G.getNumberOfCells())';

            vf0 = model.(elyte).volumeFraction;
            vols = model.(elyte).G.getVolumes();

            for ielde = 1 : numel(eldes)

                elde = eldes{ielde};

                if model.(elde).coatingModelSetup.swelling
                    
                    elde_cells = elyte_cells(model.(elde).(co).G.mappings.cellmap);
                    
                    porosity0 = 1 - state0.(elde).(co).volumeFraction;
                    vf0(elde_cells) = porosity0;
                    
                end
                
            end
            
            state.(elyte).massAccum  = vols.*(vf.*c - vf0.*c0)/dt;
            
        end

        function state = updateSwellingElectrolyteConvFlux(model, state)

            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            co    = 'Coating';
            am    = 'ActiveMaterial';
            itf   = 'Interface';
            sd    = 'SolidDiffusion';

            eldes = {ne, pe};

            % ad-hoc AD compatible initialization
            j = state.(elyte).j;
            state.(elyte).convFlux = 0 .* j;

            for ielde = 1 : numel(eldes)

                elde = eldes{ielde};

                if model.(elde).coatingModelSetup.swelling

                    G  = model.(elyte).G;

                    F = model.(elde).(co).(am).(itf).constants.F;
                    s = -1;
                    n = model.(elde).(co).(am).(itf).numberOfElectronsTransferred;

                    x = state.(elde).(co).(am).(sd).x;
                    
                    molarVolumeLithiated   = model.(elde).(co).computeMolarVolumeLithiated(x);
                    molarVolumeDelithiated = model.(elde).(co).computeMolarVolumeLithiated(0);

                    elyte_cells = zeros(model.G.getNumberOfCells(), 1);
                    elyte_cells(model.(elyte).G.mappings.cellmap) = (1 : model.(elyte).G.getNumberOfCells)';
                    elyte_cells_elde = elyte_cells(model.(elde).G.mappings.cellmap);
                    
                    j = state.(elyte).j;
                    c = state.(elyte).c;

                    flux = c(elyte_cells_elde).*(s./(n.*F)).*(molarVolumeLithiated - molarVolumeDelithiated).*j(elyte_cells_elde);

                    state.(elyte).convFlux = subsasgnAD(state.(elyte).convFlux, elyte_cells_elde, flux);

                end
                
            end
            
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

            sp = battery.Electrolyte.species;
            z = sp.chargeNumber;

            couplingterms = battery.couplingTerms;

            elyte_e_source = zeros(battery.(elyte).G.getNumberOfCells(), 1);

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

            elyte_c_source = elyte_e_source./(z*F);

            state.Electrolyte.eSource    = elyte_e_source;
            state.Electrolyte.massSource = elyte_c_source;

        end


        function state = updateControl(model, state, drivingForces)

            ctrl = "Control";

            switch model.(ctrl).controlPolicy

              case {'CCCV', 'powerControl'}

                % nothing to do here

              case {'timeControl'}

                [ctrlVal, ctrlType] = model.(ctrl).computeInput(state.time);

                state.(ctrl).ctrlVal  = ctrlVal;
                state.(ctrl).ctrlType = ctrlType;

              case {'CCDischarge', 'CCCharge'}

                Imax = model.(ctrl).Imax;

                E    = state.(ctrl).E;
                I    = state.(ctrl).I;
                time = state.time;

                ctrlType = 'constantCurrent';

                if model.(ctrl).useCVswitch

                    [ctrlVal, ctrlType] = drivingForces.src(time, value(I), value(E), Imax);

                else
                    
                    ctrlVal = drivingForces.src(time, Imax);
                    
                end

                state.(ctrl).ctrlVal  = ctrlVal;
                state.(ctrl).ctrlType = ctrlType;
                
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

            nc = model.G.getNumberOfCells();

            src = zeros(nc, 1);

            for ind = 1 : numel(eldes)

                elde = eldes{ind};

                if model.include_current_collectors

                    cc_model = model.(elde).(cc);
                    cc_map   = cc_model.G.mappings.cellmap;
                    cc_j     = state.(elde).(cc).jFace;
                    cc_econd = cc_model.effectiveElectronicConductivity;
                    cc_vols  = cc_model.G.getVolumes();
                    cc_jsq   = cc_model.G.getCellFluxNorm(cc_j);

                    src = subsetPlus(src, cc_vols .* cc_jsq ./ cc_econd, cc_map);
                    
                end

                co_model = model.(elde).(co);
                co_map   = co_model.G.mappings.cellmap;
                co_j     = state.(elde).(co).jFace;
                co_econd = co_model.effectiveElectronicConductivity;
                co_vols  = co_model.G.getVolumes();
                co_jsq   = co_model.G.getCellFluxNorm(co_j);

                src = subsetPlus(src, co_vols.*co_jsq./co_econd, co_map);

            end

            % Electrolyte
            elyte_model = model.(elyte);
            elyte_map   = elyte_model.G.mappings.cellmap;
            elyte_j     = state.(elyte).jFace;
            elyte_econd = state.(elyte).conductivity; % effective conductivity
            elyte_vols  = elyte_model.G.getVolumes();
            elyte_jsq   = elyte_model.G.getCellFluxNorm(elyte_j);

            src = subsetPlus(src, elyte_vols.*elyte_jsq./elyte_econd, elyte_map);

            state.(thermal).jHeatOhmSource = src;

        end

        function state = updateThermalChemicalSourceTerms(model, state)
        % Assemble the thermal source term from transport :code:`state.jHeatChemicalSource`, see :cite:t:`Latz2016`

            elyte   = 'Electrolyte';
            thermal = 'ThermalModel';

            % prepare term
            nc       = model.G.getNumberOfCells();
            T        = state.(thermal).T;
            nf       = model.(elyte).G.getNumberOfFaces();
            intfaces = model.(elyte).G.getIntFaces;

            phi = state.(elyte).phi;

            zeroFace = model.AutoDiffBackend.convertToAD(zeros(nf, 1), phi);

            src = zeros(nc, 1);

            % Compute chemical heat source in electrolyte
            dmudcs = state.(elyte).dmudcs;   % Derivative of chemical potential with respect to concentration
            D      = state.(elyte).D;        % Effective diffusion coefficient
            Dgradc = state.(elyte).diffFlux; % Diffusion flux (-D*grad(c))
            DFaceGradc = zeroFace;
            DFaceGradc = subsetPlus(DFaceGradc, Dgradc, intfaces);

            % Compute norm of square norm of diffusion flux
            elyte_model   = model.(elyte);
            elyte_map     = elyte_model.G.mappings.cellmap;
            elyte_vols    = elyte_model.G.getVolumes();
            elyte_jchemsq = elyte_model.G.getCellFluxNorm(DFaceGradc);
            elyte_src     = elyte_vols.*elyte_jchemsq./D;

            % This is a bit hacky for the moment (we should any way consider all the species)
            elyte_src = dmudcs{1}.*elyte_src;

            % map to source term at battery level
            src = subsetPlus(src, elyte_src, elyte_map);

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

            nc = model.G.getNumberOfCells();

            src = zeros(nc, 1);

            T = state.(thermal).T;

            for ind = 1 : numel(eldes)

                elde = eldes{ind};

                F      = model.(elde).(co).(am).(itf).constants.F;
                n      = model.(elde).(co).(am).(itf).numberOfElectronsTransferred;
                co_map = model.(elde).(co).G.mappings.cellmap;
                vsa    = model.(elde).(co).(am).(itf).volumetricSurfaceArea;
                vols   = model.(elde).(co).G.getVolumes();

                Rvol = state.(elde).(co).(am).(sd).Rvol;
                dUdT = state.(elde).(co).(am).(itf).dUdT;
                eta  = state.(elde).(co).(am).(itf).eta;

                itf_src = n*F*vols.*Rvol.*eta;

                src = subsetPlus(src, itf_src, co_map);

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

            nc = model.G.getNumberOfCells();

            T = state.(thermal).T;
            
            src = zeros(nc, 1);

            for ind = 1 : numel(eldes)

                elde = eldes{ind};

                F      = model.(elde).(co).(am).(itf).constants.F;
                n      = model.(elde).(co).(am).(itf).numberOfElectronsTransferred;
                co_map = model.(elde).(co).G.mappings.cellmap;
                vsa    = model.(elde).(co).(am).(itf).volumetricSurfaceArea;
                vols   = model.(elde).(co).G.getVolumes();

                Rvol = state.(elde).(co).(am).(sd).Rvol;
                dUdT = state.(elde).(co).(am).(itf).dUdT;
                eta  = state.(elde).(co).(am).(itf).eta;

                itf_src = n*F*vols.*Rvol.*T(co_map).*dUdT;

                src = subsetPlus(src, itf_src, co_map);

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

            nc = model.G.getNumberOfCells();

            src = zeros(nc, 1);

            T = state.(thermal).T;

            for ind = 1 : numel(eldes)

                elde = eldes{ind};

                F      = model.(elde).(co).(am).(itf).constants.F;
                n      = model.(elde).(co).(am).(itf).numberOfElectronsTransferred;
                co_map = model.(elde).(co).G.mappings.cellmap;
                vsa    = model.(elde).(co).(am).(itf).volumetricSurfaceArea;
                vols   = model.(elde).(co).G.getVolumes();

                Rvol = state.(elde).(co).(am).(sd).Rvol;
                eta  = state.(elde).(co).(am).(itf).eta;

                if model.(elde).(co).(am).(itf).includeEntropyChange
                    dUdT = state.(elde).(co).(am).(itf).dUdT;
                    itf_src = n*F*vols.*Rvol.*(eta + T(co_map).*dUdT);
                else
                    itf_src = n*F*vols.*Rvol.*eta;
                end
                
                src = subsetPlus(src, itf_src, co_map);

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
            sr    = 'SideReaction';
            sei   = 'SolidElectrodeInterface';

            eldes = {ne, pe};
            phi_elyte = state.(elyte).phi;
            c_elyte = state.(elyte).c;

            elyte_cells = zeros(model.G.getNumberOfCells(), 1);
            elyte_cells(bat.(elyte).G.mappings.cellmap) = (1 : bat.(elyte).G.getNumberOfCells())';

            for ind = 1 : numel(eldes)

                elde = eldes{ind};

                if  model.(elde).(co).activeMaterialModelSetup.composite
                    ams = {am1, am2};
                else
                    ams = {am};
                end

                for iam = 1 : numel(ams)

                    amc = ams{iam};
                    state.(elde).(co).(amc).(itf).phiElectrolyte = phi_elyte(elyte_cells(bat.(elde).(co).G.mappings.cellmap));
                    state.(elde).(co).(amc).(itf).cElectrolyte   = c_elyte(elyte_cells(bat.(elde).(co).G.mappings.cellmap));

                    switch model.(elde).(co).activeMaterialModelSetup.SEImodel
                      case {'none', 'Bolay'}
                        % nothing more to do
                      case 'Safari'
                        state.(elde).(co).(amc).(sr).phiElectrolyte = phi_elyte(elyte_cells(bat.(elde).(co).G.mappings.cellmap));
                      otherwise
                        error('SEI model not recognized');
                    end

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

                coupterm = model.(ne).(cc).externalCouplingTerm;
                bcfaces = coupterm.couplingfaces;
                [jExternal, jFaceExternal] = assembleBoundarySource(model.(ne).(cc), phi, 0, sigma, bcfaces);

                state.(ne).(cc).jExternal     = jExternal;
                state.(ne).(cc).jFaceExternal = jFaceExternal;
                state.(ne).(co).jExternal     = 0;
                state.(ne).(co).jFaceExternal = 0;

            else

                phi   = state.(ne).(co).phi;
                sigma = state.(ne).(co).conductivity;

                coupterm = model.(ne).(co).externalCouplingTerm;
                bcfaces = coupterm.couplingfaces;
                [jExternal, jFaceExternal] = assembleBoundarySource(model.(ne).(co), phi, 0, sigma, bcfaces);

                state.(ne).(co).jExternal     = jExternal;
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

                coupterm = model.(pe).(cc).externalCouplingTerm;
                bcfaces = coupterm.couplingfaces;
                [jExternal, jFaceExternal] = assembleBoundarySource(model.(pe).(cc), phi, E, sigma, bcfaces);

                state.(pe).(cc).jExternal     = jExternal;
                state.(pe).(cc).jFaceExternal = jFaceExternal;
                state.(pe).(co).jExternal     = 0;
                state.(pe).(co).jFaceExternal = 0;
            else

                phi   = state.(pe).(co).phi;
                sigma = state.(pe).(co).conductivity;

                coupterm = model.(pe).(co).externalCouplingTerm;
                bcfaces = coupterm.couplingfaces;
                [jExternal, jFaceExternal] = assembleBoundarySource(model.(pe).(co), phi, E, sigma, bcfaces);

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
            [trans_pcc, cells] = model.(pe).(mat).G.getBcTrans(faces);

            state.Control.EIequation = sum(cond_pcc.*trans_pcc.*(phi(cells) - E)) - I;


        end


        function forces = getValidDrivingForces(model)

            forces = getValidDrivingForces@PhysicalModel(model);

            forces.src = [];

            ctrl = 'Control';
            ctrlpol = model.(ctrl).controlPolicy;
            forces.(ctrlpol) = true;

            % TODO this is a hack to get thing go
            forces.Imax = [];

        end

        function model = setTPFVgeometry(model, tPFVgeometry)
        % tPFVgeometry should be instance of TwoPointFiniteVolumeGeometry

            model.G.parentGrid.tPFVgeometry = tPFVgeometry;

            elyte   = 'Electrolyte';
            sep     = 'Separator';
            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            thermal = 'ThermalModel';

            model.(elyte) = model.(elyte).setTPFVgeometry(tPFVgeometry);
            model.(sep)   = model.(sep).setTPFVgeometry(tPFVgeometry);
            model.(ne)    = model.(ne).setTPFVgeometry(tPFVgeometry);
            model.(pe)    = model.(pe).setTPFVgeometry(tPFVgeometry);

            if model.use_thermal
                model.(thermal) = model.(thermal).setTPFVgeometry(tPFVgeometry);
            end

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

                if  model.(elde).(co).activeMaterialModelSetup.composite
                    ams = {am1, am2};
                else
                    ams = {am};
                end

                for iam = 1 : numel(ams)

                    amc = ams{iam};
                    cmax = model.(elde).(co).(amc).(itf).saturationConcentration;
                    switch model.(elde).(co).(amc).diffusionModelType
                      case 'simple'
                        state.(elde).(co).(amc).(sd).cAverage = max(cmin, state.(elde).(co).(amc).(sd).cAverage);
                        state.(elde).(co).(amc).(sd).cAverage = min(cmax, state.(elde).(co).(amc).(sd).cAverage);
                      case {'full', 'swelling'}
                        state.(elde).(co).(amc).(sd).c = max(cmin, state.(elde).(co).(amc).(sd).c);
                        state.(elde).(co).(amc).(sd).c = min(cmax, state.(elde).(co).(amc).(sd).c);
                      otherwise
                        error('diffusionModelType not recognized')
                    end

                end
                
            end
            
            report = [];

        end

        function cleanState = addStaticVariables(model, cleanState, state)
        % Variables that are no AD initiated (but should be "carried over")

            cleanState = addStaticVariables@BaseModel(model, cleanState, state);

            cleanState.time = state.time;

            thermal = 'ThermalModel';

            if ~model.use_thermal
                thermal = 'ThermalModel';
                cleanState.(thermal).T = state.(thermal).T;
            end


        end

        function [model, state] = prepareTimestep(model, state, state0, dt, drivingForces)

            [model, state] = prepareTimestep@BaseModel(model, state, state0, dt, drivingForces);

            ctrl = 'Control';

            state.(ctrl) = model.(ctrl).prepareStepControl(state.(ctrl), state0.(ctrl), dt);
            
            if strcmp(model.(ctrl).controlPolicy, 'powerControl')
                state.(ctrl).time = state.time;
            end

        end

        function model = setupCapping(model)

            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            co      = 'Coating';
            itf     = 'Interface';

            eldes = {pe, ne};

            for ielde = 1 : numel(eldes)

                elde = eldes{ielde};

                if model.(elde).(co).activeMaterialModelSetup.composite

                    am1 = 'ActiveMaterial1';
                    am2 = 'ActiveMaterial2';

                    ams = {am1, am2};

                    cmax = 0;

                    for iam = 1 : numel(ams)

                        amc = ams{iam};
                        cmax = max(cmax, model.(elde).(co).(amc).(itf).saturationConcentration);

                    end

                    cmaxs{ielde} = cmax;

                else

                    am = 'ActiveMaterial';

                    cmaxs{ielde} = model.(elde).(co).(am).(itf).saturationConcentration;


                end


            end

            model.cmin = 0;


        end

        function [state, report] = updateAfterConvergence(model, state0, state, dt, drivingForces)

             [state, report] = updateAfterConvergence@BaseModel(model, state0, state, dt, drivingForces);

             ctrl = 'Control';
             state.(ctrl) = model.(ctrl).updateControlAfterConvergence(state.(ctrl), state0.(ctrl), dt);
        end


        function [state, report] = stepFunction(model, state, state0, dt, drivingForces, linsolver, nonlinsolver, iteration, varargin)

            [state, report] = stepFunction@BaseModel(model, state, state0, dt, drivingForces, linsolver, nonlinsolver, iteration, varargin{:});

            if ~report.Failure
                ctrl = 'Control';
                state.(ctrl) = model.(ctrl).updateControlState(state.(ctrl), state0.(ctrl), dt);
            end
            
            if report.Converged

                if ismember(model.(ctrl).controlPolicy, {'CCCV'})
                    % we check for the constraints

                    [arefulfilled, state.(ctrl)] = model.(ctrl).checkConstraints(state.(ctrl), state0.(ctrl), dt);

                    if ~arefulfilled
                        report.Converged = false;
                        report.Failure   = true;
                        report.FailureMsg = 'The current time step converged but for a control that is no longer valid.';
                    end

                end

            end

        end

        function outputvars = extractGlobalVariables(model, states)

            ns = numel(states);

            if ns == 0
                % This happens when simulation fail to converge at first step
                outputvars = [];
            else
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

    end

    methods(Static)

        function [found, varind] = getVarIndex(varname, pvarnames)

            varname   = GenericBattery.varToStr(varname);
            pvarnames = cellfun(@(name) GenericBattery.varToStr(name), pvarnames, 'uniformoutput', false);

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
Copyright 2021-2024 SINTEF Industry, Sustainable Energy Technology
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
