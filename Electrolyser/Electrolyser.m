classdef Electrolyser < PhysicalModel
    
    properties
        
        con = PhysicalConstants();

        % Components
        IonomereMembrane
        HydrogenElectrode
        HydrogenCatalyser
        OxygenElectrode
        OxygenCatalyser
        
        couplingTerms
        couplingNames
        
    end
    
    methods
        
        function model = Electrolyser(paramobj)

            
            model = model@PhysicalModel([]);
            
            % OBS : All the submodels should have same backend (this is not assigned automaticallly for the moment)
            model.AutoDiffBackend = SparseAutoDiffBackend('useBlocks', false);
            
            %% Setup the model using the input parameters
            fdnames = {'G'            , ...
                       'couplingTerms', ...
                       'initT'};
            
            model = dispatchParams(model, paramobj, fdnames);
            
            % Assign the components : Electrolyte, NegativeElectrode, PositiveElectrode
            model.HydrogenElectrode = OpenElectrode(paramobj.HydrogenElectrode);
            model.OxygenElectrode   = OpenElectrode(paramobj.OxygenElectrode);
            model.HydrogenCatalyser = Catalyser(paramobj.HydrogenCatalyser);
            model.OxygenCatalyser   = Catalyser(paramobj.OxygenCatalyser);
            model.IonomerMembrane   = IonomerMembrane(paramobj.IonomerMembrane);
            
            % setup couplingNames
            model.couplingNames = cellfun(@(x) x.name, model.couplingTerms, 'uniformoutput', false);
            
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
        end

        function state = dispatchTemperature(model, state)
            mea   = 'IonomereMembrane';
            her   = 'HydrogenElectrode';
            hcatl = 'HydrogenCatalyser';
            oer   = 'OxygenElectrode';
            ocatl = 'OxygenCatalyser';
            
            submodelnames = {mea, her, hcatl, oer, ocatl};
            for ind = 1 : numel(submodelnames)
                mname = submodelnames{ind};
                state.(mname).T = state.T(model.(mname).G.cellmap);
            end
        end
        
        function state = dispatchToCatalyser(model, state)
            
            state = model.updateHydrogenPotentialDiffs(state); 
            state = model.updateOxygenPotentialDiffs(state);
            state = model.dispatchGasPartialPressureHydrogenCatalyser(state);
            state = model.dispatchGasPartialPressureOxygenCatalyser(state);
            state = model.dispatchConcentrationsHydrogenCatalyser(state);
            state = model.dispatchConcentrationsOxygenCatalyser(state);
            state = model.dispatchH2OactivityHydrogenCatalyser(state);
            state = model.dispatchH2OactivityOxygenCatalyser(state);
            
        end
        
        function state = dispatchFromCatalyser(model, state)
            
            state = model.dispatchCatalyserToHydrogenElectrode(state);
            state = model.dispatchCatalyserToOxygenElectrode(state);
            state = model.dipatchCatalysersToIonomer(state)
            
        end
        

        function state = updateHydrogenPotentialDiffs(model, state)
            mea   = 'IonomereMembrane';
            her   = 'HydrogenElectrode';
            hcatl = 'HydrogenCatalyser';
            
            coupterms = model.couplingTerms;
            coupnames = model.couplingNames;
            
            coupterm = getCoupTerm(coupterms, {her, hcatl}, coupnames);
            coupcells = coupterm.couplingcells;
            state.(hcatl).phiElyte(coupcells(:, 2)) = state.(her).phi(coupcells(:, 1));
            
            coupterm = getCoupTerm(coupterms, {mea, hcatl}, coupnames);
            coupcells = coupterm.couplingcells;
            state.(hcatl).phiInmr(coupcells(:, 2)) = state.(mea).phi(coupcells(:, 1));
        end
        
        
        function state = updateOxygenPotentialDiffs(model, state)
            mea   = 'IonomereMembrane';
            oer   = 'OxygenElectrode';
            ocatl = 'OxygenCatalyser';
            
            coupterms = model.couplingTerms;
            coupnames = model.couplingNames;
            
            coupterm = getCoupTerm(coupterms, {oer, ocatl}, coupnames);
            coupcells = coupterm.couplingcells;
            state.(ocatl).phiElyte(coupcelsl(:, 2)) = state.(oer).phi(coupcells(:, 1));
            
            coupterm = getCoupTerm(coupterms, {mea, ocatl}, coupnames);
            coupcells = coupterm.couplingcells;
            state.(ocatl).phiInmr(coupcells(:, 2)) = state.(mea).phi(coupcells(:, 1));
        end

        function state = dispatchGasPartialPressureHydrogenCatalyser(model, state)
            her   = 'HydrogenElectrode';
            hcatl = 'HydrogenCatalyser';
            
            gind = model.(her).phaseInd.gas;
            m  = state.(her).compGasMasses;
            T  = state.(her).T;
            vg = state.(her).volumeFraction{gind};
            
            gInd = model.(her).gasInd;
            
            % compute mol amount for each gas component
            coupterm = getCoupTerm(coupterms, {oer, ocatl}, coupnames);
            coupcells = coupterm.couplingcells;
            
            mH2O = m{gInd.H2Ogas};
            mH2 = m{gInd.H2}
            
            mH2O = mH2O(coupcells(:, 1));
            mH2  = mH2(coupcells(:, 1));
            vg   = vg(coupcells(:, 1));
            T    = T(coupcells(:, 1));
            
            p{gInd.H2Ogas} = mH2O./(model.(her).H2O.MW).*R.*T./vg;
            p{gInd.H2} = mH2./(model.(her).H2.MW).*R.*T./vg;
            
            state.(hcatl).gasPressuresElyte(coupcells(:, 2)) = p;
        end
        
        function state = dispatchGasPartialPressureOxygenCatalyser(model, state)
            oer   = 'OxygenElectrode';
            ocatl = 'OxygenCatalyser';
            
            gind = model.(oer).phaseInd.gas;
            m  = state.(oer).compGasMasses;
            T  = state.(oer).T;
            vg = state.(oer).volumeFraction{gind};
            
            gInd = model.(oer).gasInd;
            
            % compute mol amount for each gas component
            coupterm = getCoupTerm(coupterms, {oer, ocatl}, coupnames);
            coupcells = coupterm.couplingcells;
            
            mH2O = m{gInd.H2Ogas};
            mO2 = m{gInd.O2}
            
            mH2O = mH2O(coupcells(:, 1));
            mO2  = mO2(coupcells(:, 1));
            vg   = vg(coupcells(:, 1));
            T    = T(coupcells(:, 1));
            
            p{gInd.H2Ogas} = mH2O./(model.(oer).H2O.MW).*R.*T./vg;
            p{gInd.O2} = mO2./(model.(oer).O2.MW).*R.*T./vg;
            
            state.(ocatl).gasPressuresElyte(coupcells(:, 2)) = p;
        end
        
        function state = dispatchConcentrationsHydrogenCatalyser(model, state)
            her   = 'HydrogenElectrode';
            hcatl = 'HydrogenCatalyser';
            
            lind = model.(her).liquidInd;
            
            cH2O = state.(her).concentrations{lind.H2Oliquid};
            cOH = state.(her).concentrations{lind.OH};
            
            coupterms = model.couplingTerms;
            coupnames = model.couplingNames;
            
            coupterm = getCoupTerm(coupterms, {her, hcatl}, coupnames);
            coupcells = coupterm.couplingcells;
            cOH = cOH(coupcells(:, 1));
            cH = cH(coupcells(:, 1));

            state.(hcatl).cOHElyte(coupcells(:, 2)) = cOH;
            state.(hcatl).cHElyte(coupcells(:, 2)) = cH;
        end
        
        function state = dispatchConcentrationsOxygenCatalyser(model, state)
            oer   = 'OxygenElectrode';
            ocatl = 'OxygenCatalyser';
            
            lind = model.(oer).liquidInd;
            
            cH2O = state.(oer).concentrations{lind.H2Oliquid};
            cOH = state.(oer).concentrations{lind.OH};
            
            coupterms = model.couplingTerms;
            coupnames = model.couplingNames;
            
            coupterm = getCoupTerm(coupterms, {oer, ocatl}, coupnames);
            coupcells = coupterm.couplingcells;
            cOH = cOH(coupcells(:, 1));
            cH = cH(coupcells(:, 1));

            state.(ocatl).cOHElyte(coupcells(:, 2)) = cOH;
            state.(ocatl).cHElyte(coupcells(:, 2)) = cH;
        end
        
        function state = dispatchH2OactivityHydrogenCatalyser(model, state)
            her   = 'HydrogenElectrode';
            hcatl = 'HydrogenCatalyser';
            mea   = 'IonomerMembrane';
            
            coupterms = model.couplingTerms;
            coupnames = model.couplingNames;
            
            coupterm = getCoupTerm(coupterms, {her, hcatl}, coupnames);
            coupcells = coupterm.couplingcells;
            state.(hcatl).H2OaElyte(coupcells(:, 2)) = state.(her).H2Oactivity(coupcells(:, 1));
            
            coupterm = getCoupTerm(coupterms, {mea, hcatl}, coupnames);
            coupcells = coupterm.couplingcells;
            state.(hcatl).H2OaInmr(coupcells(:, 2)) = state.(mea).H2Oactivity(coupcells(:, 1));
        end
                    
        
        function state = dispatchH2OactivityOxygenCatalyser(model, state)
            oer   = 'OxygenElectrode';
            ocatl = 'OxygenCatalyser';
            mea   = 'IonomerMembrane';
            
            coupterms = model.couplingTerms;
            coupnames = model.couplingNames;
            
            coupterm = getCoupTerm(coupterms, {oer, ocatl}, coupnames);
            coupcells = coupterm.couplingcells;
            state.(ocatl).H2OaElyte(coupcells(:, 2)) = state.(oer).H2Oactivity(coupcells(:, 1));
            
            coupterm = getCoupTerm(coupterms, {mea, ocatl}, coupnames);
            coupcells = coupterm.couplingcells;
            state.(ocatl).H2OaInmr(coupcells(:, 2)) = state.(mea).H2Oactivity(coupcells(:, 1));
        end
        
        
        function state = dispatchCatalyserToHydrogenElectrode(model, state)
            
            her   = 'HydrogenElectrode';
            hcatl = 'HydrogenCatalyser';
            mea   = 'IonomerMembrane';
            
            leps = model.(mea).liquidVolumeFraction;
            gInd = model.(her).gasInd;
            
            coupterms = model.couplingTerms;
            coupnames = model.couplingNames;
            
            coupterm = getCoupTerm(coupterms, {her, hcatl}, coupnames);
            coupcells = coupterm.couplingcells;
            
            jElyte    = state.(hcatl).jElyte(coupcells(:, 2));
            j         = state.(hcatl).j(coupcells(:, 2));
            Rxch      = state.(hcatl).Rxch(coupcells(:, 2));
            Rsorption = state.(hcatl).Rsorption(coupcells(:, 2));
            
            OHSource        = 2*jElyte./(n*F) + Rxch.*leps;
            H2OliquidSource = -Rsorption;
            compH2GasSource = -j/(n*F).*H2.MW;
            
            state.(her).OHSource(coupcells(:, 1)) = OHSource;
            state.(her).H2OliquidSource(coupcells(:, 1)) = H2OliquidSource;
            state.(her).compGasSources{gInd.H2}(coupcells(:, 1)) = compH2GasSource;
        end
        
        function state = dispatchCatalyserToOxygenElectrode(model, state)
            
            oer   = 'OxygenElectrode';
            ocatl = 'OxygenCatalyser';
            mea   = 'IonomerMembrane';
            
            leps = model.(mea).liquidVolumeFraction;
            gInd = model.(oer).gasInd;
            sp   = model.(oer).sp;
            
            coupterms = model.couplingTerms;
            coupnames = model.couplingNames;
            coupterm = getCoupTerm(coupterms, {oer, ocatl}, coupnames);
            coupcells = coupterm.couplingcells;
            
            jElyte    = state.(ocatl).jElyte(coupcells(:, 2));
            j         = state.(ocatl).j(coupcells(:, 2));
            Rxch      = state.(ocatl).Rxch(coupcells(:, 2));
            Rsorption = state.(ocatl).Rsorption(coupcells(:, 2));
            
            OHSource        = 2*jElyte./(n*F) + Rxch.*leps;
            H2OliquidSource = -Rsorption;
            compO2GasSource = -0.5.*j./(n*F).*sp.O2.MW;
            
            state.(oer).OHSource(coupcells(:, 1)) = OHSource;
            state.(oer).H2OliquidSource(coupcells(:, 1)) = H2OliquidSource;
            state.(oer).compGasSources{gInd.O2}(coupcells(:, 1)) = compO2GasSource;
        end
        
        function state = dipatchCatalysersToIonomer(model, state)
            
            % initialize sources
            OHSource = 0*state.(mea).phi;
            H2OSource = 0*state.(mea).phi;

            
            % Dispatch values from Hydrogen Catalyser
            coupterms = model.couplingTerms;
            coupnames = model.couplingNames;
            coupterm = getCoupTerm(coupterms, {mea, hcatl}, coupnames);
            coupcells = coupterm.couplingcells;
            jInmr = state.(hcatl).jInmr(coupcells(:, 2));
            Rxch  = state.(hcatl).Rxch(coupcells(:, 2));
            leps  = model.(mea).liquidVolumeFraction(coupcells(:, 2));
            
            OHSource(coupcells(:, 1))  = -2.*jInmr./(n.*F) - Rxch.*leps;
            H2OSource(coupcells(:, 1)) = -jInmr./(n.*F) - Rsorption;
            
            % Dispatch values from Oxygen Catalyser
            coupterms = model.couplingTerms;
            coupnames = model.couplingNames;
            coupterm = getCoupTerm(coupterms, {mea, ocatl}, coupnames);
            coupcells = coupterm.couplingcells;
            jInmr = state.(ocatl).jInmr(coupcells(:, 2));
            Rxch  = state.(ocatl).Rxch(coupcells(:, 2));
            leps  = model.(mea).liquidVolumeFraction(coupcells(:, 2));

            OHSource(coupcells(:, 1))  = -2.*jInmr./(n.*F) - Rxch.*leps;
            H2OSource(coupcells(:, 1)) = -jInmr./(n.*F) - Rsorption;
            
            state.(mea).OHSource  = OHSource;
            state.(mea).H2OSource = H2OSource;
        end
        
        
        function [problem, state] = getEquations(model, state0, state,dt, drivingForces, varargin)
            
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
            thermal = 'ThermalModel';
            
            electrodes = {ne, pe};
            electrodecomponents = {eac, cc};

            %% Synchronization across components

            % temperature
            state = battery.updateTemperature(state);
            
            for ind = 1 : numel(electrodes)
                elde = electrodes{ind};
                % potential and concentration between active material and electode active component
                state.(elde).(eac) = battery.(elde).(eac).updatePhi(state.(elde).(eac));
                state.(elde).(eac) = battery.(elde).(eac).updateChargeCarrier(state.(elde).(eac));
                
            end
            
            %% Accumulation terms

            state = battery.updateAccumTerms(state, state0, dt);

            %% Update Electrolyte -> Electrodes coupling 
            
            state = battery.updateElectrodeCoupling(state); 

            %% Update reaction rates in both electrodes

            for ind = 1 : numel(electrodes)
                elde = electrodes{ind};
                state.(elde).(eac).(am) = battery.(elde).(eac).(am).updateDiffusionConductivityCoefficients(state.(elde).(eac).(am));
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
            
            
            %% Set up the governing equations
            
            eqs = {};
            
            %% We collect mass and charge conservation equations for the electrolyte and the electrodes

            eqs{end + 1} = state.(elyte).massCons;
            eqs{end + 1} = state.(elyte).chargeCons;
            
            eqs{end + 1} = state.(ne).(eac).massCons;
            eqs{end + 1} = state.(ne).(eac).chargeCons;
            
            eqs{end + 1} = state.(pe).(eac).massCons;
            eqs{end + 1} = state.(pe).(eac).chargeCons;
            
            eqs{end + 1} = state.(ne).(cc).chargeCons;
            eqs{end + 1} = state.(pe).(cc).chargeCons;
            
            eqs{end + 1} = state.(thermal).energyCons;
            
            %% We setup and add the control equation (fixed total current at PositiveCurrentCollector)
            
            src = drivingForces.src(time);
            coupterm = battery.(pe).(cc).couplingTerm;
            faces = coupterm.couplingfaces;
            bcval = state.(pe).(cc).E;
            cond_pcc = battery.(pe).(cc).EffectiveElectricalConductivity;
            [trans_pcc, cells] = battery.(pe).(cc).operators.harmFaceBC(cond_pcc, faces);
            control = sum(trans_pcc.*(state.(pe).(cc).phi(cells) - bcval)) - src;
            
            eqs{end + 1} = - control;

            %% Give type and names to equations and names of the primary variables (for book-keeping)
            
            types = {'cell','cell','cell','cell',...
                     'cell','cell','cell','cell','cell','cell'};
            
            names = {'elyte_massCons'   , ...
                     'elyte_chargeCons' , ...
                     'ne_eac_massCons'  , ...
                     'ne_eac_chargeCons', ...
                     'pe_eac_massCons'  , ...
                     'pe_eac_chargeCons', ...
                     'ne_cc_chargeCons' , ...
                     'pe_cc_chargeCons' , ...
                     'energyCons'       , ...
                     'control'};
            
            primaryVars = model.getPrimaryVariables();

            %% setup LinearizedProblem that can be processed by MRST Newton API
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
        
        end

        function state = updateTemperature(model, state)
            
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
        
        
        function state = updateElectrolyteCoupling(model, state)
        % Setup the electrolyte coupling by adding ion sources from the electrodes
        % shortcuts:
        % c_source : Source term for charge carrier.
                        
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
            
            elyte_e_source = elyte_c_source.*battery.(elyte).sp.z(1); % we divide with F later
            
            state.Electrolyte.(ccSourceName) = elyte_c_source/F; 
            state.Electrolyte.eSource = elyte_e_source;
            
        end
        
        function state = updateAccumTerms(model, state, state0, dt)
                    
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

        
        function state = updateThermalAccumTerms(model, state, state0, dt)
                    
            thermal = 'ThermalModel';
         
            hcap = model.(thermal).EffectiveHeatCapacity;
            
            T = state.(thermal).T;
            T0 = state0.(thermal).T;

            % (here we assume that the ThermalModel has the "parent" grid)
            vols = model.G.cells.volumes;
            
            state.(thermal).accumHeat = hcap.*vols.*(T - T0)/dt;
            
        end


        function state = updateThermalOhmicSourceTerms(model, state)
        % reference Latz et al (ref1 in reference list)

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
            end
            
            for ind = 1 : numel(eldes)

                elde = eldes{ind};
                
                cc_model = model.(elde).(cc);
                cc_map   = cc_model.G.mappings.cellmap;
                cc_j     = state.(elde).(cc).j;
                cc_econd = cc_model.EffectiveElectricalConductivity;
                
                cc_src = computeCellFluxNorm(cc_model, cc_j, 1, cc_econd); % note that we send volumeFraction = 1
                src(cc_map) = src(cc_map) + cc_src;

                eac_model = model.(elde).(eac);
                eac_map   = eac_model.G.mappings.cellmap;
                eac_j     = state.(elde).(eac).j;
                eac_econd = eac_model.EffectiveElectricalConductivity;
                eac_vf    = eac_model.volumeFraction;
                
                eac_src = computeCellFluxNorm(eac_model, eac_j, eac_vf, eac_econd);
                src(eac_map) = src(eac_map) + eac_src;
                
            end

            % Electrolyte
            elyte_model = model.(elyte);
            elyte_map   = elyte_model.G.mappings.cellmap;
            elyte_j     = state.(elyte).j;
            elyte_econd = state.(elyte).conductivity;
            elyte_vf    = elyte_model.volumeFraction;            
            
            elyte_src = computeCellFluxNorm(elyte_model, elyte_j, elyte_vf, elyte_econd);
            src(elyte_map) = src(elyte_map) + elyte_src;
            
            state.(thermal).jHeatOhmSource = src;
            
            % For simplicity, we set the boundary heat transfer term equal to zero here.
            state.(thermal).jHeatBcSource = 0;
            
        end
        
        function state = updateThermalChemicalSourceTerms(model, state)
        % reference Latz et al (ref1 in reference list)            
            
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
            end            

            % Compute chemical heat source in electrolyte
            dmudcs = state.(elyte).dmudcs; % Derivative of chemical potential with respect to concentration
            D = state.(elyte).D;           % Diffusion coefficient 
            Ddc = state.(elyte).diffFlux;  % Diffusion flux (-D*grad(c))
            
            elyte_map = model.(elyte).G.mappings.cellmap;
            vf = model.(elyte).volumeFraction;
            
            elyte_src = computeCellFluxNorm(model.(elyte), Ddc, vf, D);
            % This is a bit hacky for the moment (we should any way consider all the species)
            elyte_src = dmudcs{1}.*elyte_src;
            
            % map to source term at battery level
            src(elyte_map) = src(elyte_map) + elyte_src;
            
            state.(thermal).jHeatChemicalSource = src;
            
        end
        
        
        function state = updateThermalReactionSourceTerms(model, state)
        % reference Latz et al (ref1 in reference list)            
            
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
            end
            
            for ind = 1 : numel(eldes)

                elde = eldes{ind};
                am_model = model.(elde).(eac).(am);
                am_map   = am_model.G.mappings.cellmap;
                R = state.(elde).(eac).(am).R;
                eta = state.(elde).(eac).(am).eta;
                vols = model.(elde).(eac).G.cells.volumes;
                am_src = vols.*R.*eta;
                
                src(am_map) = src(am_map) + am_src;
                
            end

            % We multiply by volumes
            src = model.G.cells.volumes.*src;
            
            state.(thermal).jHeatReactionSource = src;

        end
        
        
        function state = updateElectrodeCoupling(model, state)
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
            
            ne = 'NegativeElectrode';
            cc = 'CurrentCollector';
           
            phi = state.(ne).(cc).phi;

            jExternal = model.(ne).(cc).setupExternalCoupling(phi, 0);
            
            state.(ne).(cc).jExternal = jExternal;
            
        end
        
        function state = setupExternalCouplingPositiveElectrode(model, state)
            
            pe = 'PositiveElectrode';
            cc = 'CurrentCollector';
           
            phi = state.(pe).(cc).phi;
            E = state.(pe).(cc).E;
            
            jExternal = model.(pe).(cc).setupExternalCoupling(phi, E);
            
            state.(pe).(cc).jExternal = jExternal;
            
        end

        
        function state = initStateAD(model,state)
            
            bat = model;
            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            am    = 'ActiveMaterial';
            eac   = 'ElectrodeActiveComponent';
            cc    = 'CurrentCollector';
            thermal = 'ThermalModel';
            
            adbackend = model.AutoDiffBackend();
            useMex=false;
            if(isprop(adbackend,'useMex'))
               useMex = adbackend.useMex; 
            end
            opts=struct('types',[1,1,2,2,3,3,4,5,6,7],'useMex',useMex);
            [state.(elyte).cs{1}  , ...
             state.(elyte).phi    , ...   
             state.(ne).(eac).c   , ...   
             state.(ne).(eac).phi , ...   
             state.(pe).(eac).c   , ...    
             state.(pe).(eac).phi , ...   
             state.(ne).(cc).phi  , ...    
             state.(pe).(cc).phi  , ...    
             state.(thermal).T    , ...
             state.(pe).(cc).E] = ...
                adbackend.initVariablesAD(...
                    state.(elyte).cs{1}  , ...
                    state.(elyte).phi    , ...   
                    state.(ne).(eac).c   , ...    
                    state.(ne).(eac).phi , ...   
                    state.(pe).(eac).c   , ...    
                    state.(pe).(eac).phi , ...   
                    state.(ne).(cc).phi  , ...    
                    state.(pe).(cc).phi  , ...    
                    state.(thermal).T    , ...
                    state.(pe).(cc).E, ...
                    opts); 
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
            thermal = 'ThermalModel';
            
            p = {{elyte, 'cs', 1} , ...
                 {elyte, 'phi'}   , ...   
                 {ne, eac, 'c'}   , ...    
                 {ne, eac, 'phi'} , ...   
                 {pe, eac, 'c'}   , ...    
                 {pe, eac, 'phi'} , ...   
                 {ne, cc, 'phi'}  , ...    
                 {pe, cc, 'phi'}  , ...
                 {thermal, 'T'}   , ...
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

end

