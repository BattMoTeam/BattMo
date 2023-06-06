classdef BatterySwelling < Battery
% 
% Used when at least one of the two electrodes is a swelling Material.
% It functions as the Battery class but implements a new equation necessary because of the particle swelling :
% the volume Conservation equation 
%
    properties

    end
    
    methods
        
        function model = BatterySwelling(paramobj)
            model = model@Battery(paramobj)
        end

    %% Same setup as in the BatteryClass but adding the volumeConservationEquation
        function model = setupSelectedModel(model, varargin)
            

            opt = struct('reduction', []);
            opt = merge_options(opt, varargin{:});
            % for the reduction structrure format see battmodDir()/Utilities/JsonSchemas/linearsolver.schema.json and /reduction property
            
            elyte   = 'Electrolyte';
            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            am      = 'ActiveMaterial';
            sd      = 'SolidDiffusion';
            cc      = 'CurrentCollector';
            ctrl    = 'Control';
            thermal = 'ThermalModel';
            
            addedVariableNames = {};
            
            varEqTypes ={{elyte, 'c'}   , 'elyte_massCons'  , 'cell'; ...  
                         {elyte, 'phi'} , 'elyte_chargeCons', 'cell'; ...    
                         {ne, am, 'phi'}, 'ne_am_chargeCons', 'cell'; ...    
                         {pe, am, 'phi'}, 'pe_am_chargeCons', 'cell'; ...
                         {ctrl, 'E'}    , 'EIeq'            , 'ctrl'; ...  
                         {ctrl, 'I'}    , 'controlEq'       , 'ctrl'};

            %Specific to swelling materials
            if model.(pe).(am).isSwellingMaterial
                newentry = {{pe, am, 'porosity'}, 'pe_am_volumeCons', 'cell'};
                varEqTypes = vertcat(varEqTypes, newentry);
            end
            if model.(ne).(am).isSwellingMaterial
                newentry = {{ne, am, 'porosity'}, 'ne_am_volumeCons', 'cell'};
                varEqTypes = vertcat(varEqTypes, newentry);
            end
            
            if model.use_thermal
                newentries = {{thermal, 'T'}, 'energyCons', 'cell'};
                varEqTypes = vertcat(varEqTypes, newentries);
            else
                addedVariableNames{end + 1} = {thermal, 'T'};
            end

            switch model.(ne).(am).diffusionModelType
              case 'simple'
                newentries = {{ne, am, 'c'}, 'ne_am_massCons', 'scell';
                              {ne, am, sd, 'cSurface'}, 'ne_am_sd_soliddiffeq', 'scell'};
              case 'full'
                newentries = {{ne, am, sd, 'c'}, 'ne_am_sd_massCons', 'cell';
                              {ne, am, sd, 'cSurface'}, 'ne_am_sd_soliddiffeq', 'cell'};
              case 'interParticleOnly'
                newentries = {{ne, am, 'c'}, 'ne_am_massCons', 'cell'};
              otherwise
                error('diffusionModelType not recognized');
            end
            varEqTypes = vertcat(varEqTypes, newentries);

            
            switch model.(pe).(am).diffusionModelType
              case 'simple'
                newentries = {{pe, am, 'c'}, 'pe_am_massCons', 'scell';
                              {pe, am, sd, 'cSurface'}, 'pe_am_sd_soliddiffeq', 'scell'};
              case 'full'
                newentries = {{pe, am, sd, 'c'}, 'pe_am_sd_massCons', 'cell';
                              {pe, am, sd, 'cSurface'}, 'pe_am_sd_soliddiffeq', 'cell'};
              case 'interParticleOnly'
                newentries = {{pe, am, 'c'}, 'pe_am_massCons', 'cell'};
              otherwise
                error('diffusionModelType not recognized');
            end
            varEqTypes = vertcat(varEqTypes, newentries);

            
            if model.include_current_collectors
                
                newentries = {{ne, cc, 'phi'}, 'ne_cc_chargeCons', 'cell'; ...
                              {pe, cc, 'phi'}, 'pe_cc_chargeCons', 'cell'};
                         
                varEqTypes = vertcat(varEqTypes, newentries);
                
            end

            switch model.(ctrl).controlPolicy
              case 'IEswitch'
                addedVariableNames{end + 1} = {ctrl, 'ctrlType'};
              case 'CCCV'
                addedVariableNames{end + 1} = {ctrl, 'ctrlType'};
                addedVariableNames{end + 1} = {ctrl, 'nextCtrlType'};
              case 'CC'
                addedVariableNames{end + 1} = {ctrl, 'ctrlType'};
              otherwise
                error('controlPolicy not recognized');
            end
            
            primaryVariableNames = varEqTypes(:, 1);
            equationTypes = varEqTypes(:, 3);

            % The variable and equation lists are not a priori ordered (in the sense that we have 'cell' types first and
            % the other types after). It is a requirement in some setup of the linear solver.
            % Note : if you use a direct solver, this is not used.
            
            variableReordered = false; 
            
            if ~isempty(opt.reduction)
                reduc = opt.reduction;
                % We set the type of the variable to be reduced as 'other' and we move them at the end of list
                % respecting the order they have been given.
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
                        [found, ind] = Battery.getVarIndex(var.name, primaryVariableNames);
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


            primaryVariableNames = varEqTypes(inds, 1);
            equationNames = varEqTypes(inds, 2);
            equationTypes = equationTypes(inds);
            
            equationIndices = struct();
            for ieq = 1 : numel(equationNames)
                equationIndices.(equationNames{ieq}) = ieq;
            end
            
            model.addedVariableNames   = addedVariableNames; 
            model.primaryVariableNames = primaryVariableNames;
            model.equationNames        = equationNames; 
            model.equationTypes        = equationTypes;
            model.equationIndices      = equationIndices;
            
        end
      
     %% Update at each step the electrolyte volume fractions in the different regions (neg_elde, elyte, pos_elde)
        function state = updateElectrolyteVolumeFraction(model, state)
       

            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            am    = 'ActiveMaterial';
            sep   = 'Separator';

            elyte_cells = zeros(model.G.cells.num, 1);
            elyte_cells(model.(elyte).G.mappings.cellmap) = (1 : model.(elyte).G.cells.num)';

            
            % Initialisation of AD for the porosity of the elyte
            state.(elyte).volumeFraction = 0 * state.(elyte).c;

            % Define the porosity in the separator
            sep_cells = elyte_cells(model.(elyte).(sep).G.mappings.cellmap); 
            state.(elyte).volumeFraction = subsasgnAD(state.(elyte).volumeFraction,sep_cells, model.(elyte).(sep).porosity);


            % Define the volumeFraction in the electrodes
            eldes = {ne, pe};
            for ielde = 1 : numel(eldes)
                elde = eldes{ielde};
                if model.(elde).ActiveMaterial.isSwellingMaterial
                    state.(elyte).volumeFraction = subsasgnAD(state.(elyte).volumeFraction, elyte_cells(model.(elde).(am).G.mappings.cellmap), state.(elde).(am).porosity);
                else
                    state.(elyte).volumeFraction = subsasgnAD(state.(elyte).volumeFraction, elyte_cells(model.(elde).(am).G.mappings.cellmap), model.(elde).(am).porosity);
                end
            end

            if model.use_thermal
                model.(elyte).EffectiveThermalConductivity = model.(elyte).volumeFraction.*model.(elyte).thermalConductivity;
            end

        end

     %% Assign at each step the convective flux in the electrolyte region
        function state = updateConvFlux(model, state)

            elyte   = 'Electrolyte';
            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            am      = 'ActiveMaterial';
            itf     = 'Interface';
            eldes = {ne, pe};

            j = state.(elyte).j;
            state.(elyte).convFlux = 0 .* j;

            for ielde = 1 : numel(eldes)
                elde = eldes{ielde};
                if model.(elde).ActiveMaterial.isSwellingMaterial

                    itf_model = model.(elde).(am).(itf);
                    G = model.(elyte).G;
                    Gp = G.mappings.parentGrid;

                    c = state.(elde).(am).SolidDiffusion.cAverage;

                    molarVolumeLithiated = model.(elde).(am).updateMolarVolumeLithiated(c);
                    densitySi            = model.(elde).(am).Interface.density;
                    molarMassSi   = model.(elde).(am).molarMass;

                
                    F       = itf_model.constants.F;
                    n       = itf_model.n;
                    s       = -1;
                     molarVolumeSi = molarMassSi/densitySi;

                    a = state.(elde).(am).Interface.volumetricSurfaceArea;
                    j = state.(elyte).j;
                    j = j(model.(elde).G.mappings.cellmap);
                    c = state.(elyte).c;
                    c = c(model.(elde).G.mappings.cellmap);

                    elyte_cells = zeros(Gp.cells.num-1, 1);
                    elyte_cells(G.mappings.cellmap) = (1 : model.G.cells.num)';
                    elyte_cells_elde = elyte_cells(model.(elde).G.mappings.cellmap);

                    averageVelocity = (s./(n.*F)).*(molarVolumeLithiated - (4/15)*molarVolumeSi).*j;
                    Flux = c .* averageVelocity;


                    state.(elyte).convFlux(elyte_cells_elde) = Flux;
                end
            end
        end

      %% Same initialisation as for Battery but includes the porosity initialisation
        function initstate = setupInitialState(model)
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

            %% Setup initial state for electrodes
            
            eldes = {ne, pe};
            
            for ind = 1 : numel(eldes)
                
                elde = eldes{ind};
                
                elde_itf = bat.(elde).(am).(itf); 

                theta = SOC*(elde_itf.theta100 - elde_itf.theta0) + elde_itf.theta0;
                c     = theta*elde_itf.cmax;

                if bat.(elde).(am).isSwellingMaterial
                %need to rescale the concentration according to the initial size of the particle.

                    cmax     = elde_itf.cmax;
                    R_delith = bat.(elde).(am).(sd).rp;
    
                    %Calculating the initial radius of the particle
                    molarVolumeSi = 1.2e-05;
                    molarVolumeLi = bat.(elde).(am).constants.molarVolumeLi;
                    Q = (3.75*molarVolumeLi)/(molarVolumeSi);
                    radius = R_delith .* (1 + Q .* SOC)^(1/3);

                    c = cmax .* SOC .* (1+Q) .* (R_delith^3) ./ (radius^3);

                end

                nc    = elde_itf.G.cells.num;

                switch model.(elde).(am).diffusionModelType
                  case 'simple'
                    initstate.(elde).(am).(sd).cSurface = c*ones(nc, 1);
                    initstate.(elde).(am).c = c*ones(nc, 1);
                  case 'full'
                    initstate.(elde).(am).(sd).cSurface = c*ones(nc, 1);
                    N = model.(elde).(am).(sd).N;
                    np = model.(elde).(am).(sd).np; % Note : we have by construction np = nc
                    initstate.(elde).(am).(sd).c = c*ones(N*np, 1);
                  case 'interParticleOnly'
                    initstate.(elde).(am).c = c*ones(nc, 1);                    
                  otherwise
                    error('diffusionModelType not recognized')
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
            
            switch model.(ctrl).controlPolicy
              case 'CCCV'
                switch model.(ctrl).initialControl
                  case 'discharging'
                    initstate.(ctrl).ctrlType = 'CC_discharge1';
                    initstate.(ctrl).nextCtrlType = 'CC_discharge1';
                    initstate.(ctrl).I = model.(ctrl).Imax;
                  case 'charging'
                    initstate.(ctrl).ctrlType     = 'CC_charge1';
                    initstate.(ctrl).nextCtrlType = 'CC_charge1';
                    initstate.(ctrl).I = - model.(ctrl).Imax;
                  otherwise
                    error('initialControl not recognized');
                end
              case 'IEswitch'
                initstate.(ctrl).ctrlType = 'constantCurrent';
                switch model.(ctrl).initialControl
                  case 'discharging'
                    initstate.(ctrl).I = model.(ctrl).Imax;
                  case 'charging'
                    %initstate.(ctrl).I = model.(ctrl).Imax;
                    error('to implement (should be easy...)')
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
            for ind = 1 : numel(eldes)
                elde = eldes{ind};
                if model.(elde).(am).isSwellingMaterial
                    initstate.(elde).(am).porosity = model.(elde).(am).porosity;
                end   
            end
            
        end



      %% Assembly of the governing equation (same as for Battery but taking into account porosity variations)
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
            
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

            for ind = 1 : numel(electrodes)
                elde = electrodes{ind};

                %% Set the effective electrical conductivity if the electrode is a swelling material. added by Enguerran
                if battery.(elde).(am).isSwellingMaterial 
                    state.(elde).(am) = battery.(elde).(am).updateVolumeFraction(state.(elde).(am));
                end

                state.(elde).(am) = battery.(elde).(am).updatePhi(state.(elde).(am));            

               
             %% potential and concentration between interface and active material
                if (model.(elde).(am).use_particle_diffusion)
                    state.(elde).(am) = battery.(elde).(am).updateConcentrations(state.(elde).(am));
                else
                    state.(elde).(am).(itf).cElectrodeSurface = state.(elde).(am).c;
                end

            end
            

            %% Update the electrolyte volume fraction
            state = battery.updateElectrolyteVolumeFraction(state);
            

            %% Accumulation term in elyte

            state.(elyte) = battery.(elyte).updateAccumTerm(state.(elyte), state0.(elyte), dt);

            %% Update Electrolyte -> Electrodes coupling 
            
            state = battery.updateElectrodeCoupling(state);

            %% Update reaction rates in both electrodes

            for ind = 1 : numel(electrodes)
                elde = electrodes{ind};

                if battery.(elde).(am).isSwellingMaterial
                    state.(elde).(am).(sd)            = battery.(elde).(am).(sd).updateAverageConcentration(state.(elde).(am).(sd));
                    state.(elde).(am)                 = battery.(elde).(am).updateRadius(state.(elde).(am));
                    state.(elde).(am)                 = battery.(elde).(am).updateVolumetricSurfaceArea(state.(elde).(am));
                    state.(elde).(am)                 = battery.(elde).(am).updateReactionRateCoefficient(state.(elde).(am));
                else
                    state.(elde).(am).(sd)            = battery.(elde).(am).(sd).updateAverageConcentration(state.(elde).(am).(sd));
                    state.(elde).(am).(itf)           = battery.(elde).(am).(itf).updateReactionRateCoefficient(state.(elde).(am).(itf));
                    
                end
                
                state.(elde).(am).(itf) = battery.(elde).(am).(itf).updateOCP(state.(elde).(am).(itf));
                state.(elde).(am).(itf) = battery.(elde).(am).(itf).updateEta(state.(elde).(am).(itf));
                state.(elde).(am).(itf) = battery.(elde).(am).(itf).updateReactionRate(state.(elde).(am).(itf));
                state.(elde).(am) = battery.(elde).(am).updateRvol(state.(elde).(am));
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
                    state.(elde).(am) = battery.(elde).(am).updatejExternal(state.(elde).(am));
                else
                    state.(elde).(am) = battery.(elde).(am).updatejCoupling(state.(elde).(am));
                end
                state.(elde).(am) = battery.(elde).(am).updatejBcSource(state.(elde).(am));
                
            end
            
            %% elyte charge conservation

            state.(elyte) = battery.(elyte).updateCurrentBcSource(state.(elyte));
            state.(elyte) = battery.(elyte).updateConductivity(state.(elyte));
            state.(elyte) = battery.(elyte).updateDmuDcs(state.(elyte));
            state.(elyte) = battery.(elyte).updateCurrent(state.(elyte));
            state         = battery.updateConvFlux(state);
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


            %% Electrodes volume Conservation
            %added by Enguerran
           for ind = 1 : numel(electrodes)
                elde = electrodes{ind}; 
                if battery.(elde).(am).isSwellingMaterial
                    state.(elde).(am)                 = battery.(elde).(am).updatePorosityAccum(state.(elde).(am), state0.(elde).(am), dt);
                    state.(elde).(am)                 = battery.(elde).(am).updatePorositySource(state.(elde).(am));
                    state.(elde).(am)                 = battery.(elde).(am).updatePorosityFlux(state.(elde).(am));
                    state.(elde).(am)                 = battery.(elde).(am).updateVolumeConservation(state.(elde).(am));
                end
           end

            %% update solid diffustion equations
            for ind = 1 : numel(electrodes)
                elde = electrodes{ind};

                    switch model.(elde).(am).diffusionModelType
                      case 'simple'
                        state.(elde).(am).(sd) = battery.(elde).(am).(sd).updateDiffusionCoefficient(state.(elde).(am).(sd));
                        state.(elde).(am)      = battery.(elde).(am).assembleAccumTerm(state.(elde).(am), state0.(elde).(am), dt);
                        state.(elde).(am)      = battery.(elde).(am).updateMassSource(state.(elde).(am));
                        state.(elde).(am)      = battery.(elde).(am).updateMassFlux(state.(elde).(am));
                        state.(elde).(am)      = battery.(elde).(am).updateMassConservation(state.(elde).(am));
                        state.(elde).(am).(sd) = battery.(elde).(am).(sd).assembleSolidDiffusionEquation(state.(elde).(am).(sd));
                      case 'full'
                        state.(elde).(am).(sd) = battery.(elde).(am).(sd).updateDiffusionCoefficient(state.(elde).(am).(sd));
                        state.(elde).(am).(sd) = battery.(elde).(am).(sd).updateMassAccum(state.(elde).(am).(sd), state0.(elde).(am).(sd), dt);
                        state.(elde).(am).(sd) = battery.(elde).(am).(sd).updateMassSource(state.(elde).(am).(sd));
                        state.(elde).(am).(sd) = battery.(elde).(am).(sd).updateFlux(state.(elde).(am).(sd));
                        state.(elde).(am).(sd) = battery.(elde).(am).(sd).updateMassConservation(state.(elde).(am).(sd));
                        state.(elde).(am).(sd) = battery.(elde).(am).(sd).assembleSolidDiffusionEquation(state.(elde).(am).(sd));
                      case 'interParticleOnly'
                        state.(elde).(am) = battery.(elde).(am).assembleAccumTerm(state.(elde).(am), state0.(elde).(am), dt);
                        state.(elde).(am) = battery.(elde).(am).updateMassFlux(state.(elde).(am));
                        state.(elde).(am) = battery.(elde).(am).updateMassSource(state.(elde).(am));
                        state.(elde).(am) = battery.(elde).(am).updateMassConservation(state.(elde).(am));
                      otherwise
                        error('diffusionModelType not recognized')
                    end


     
                
            end

            if model.use_thermal
                
                %% update Face fluxes
                for ind = 1 : numel(electrodes)
                    elde = electrodes{ind};
                    state.(elde).(am) = battery.(elde).(am).updatejFaceBc(state.(elde).(am));
                    state.(elde).(am) = battery.(elde).(am).updateFaceCurrent(state.(elde).(am));
                    if model.include_current_collectors
                        state.(elde).(cc) = battery.(elde).(cc).updatejFaceBc(state.(elde).(cc));
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
                
                state.(thermal) = battery.(thermal).updateAccumTerm(state.(thermal), state0.(thermal), dt);
                
                %% Update energy conservation residual term
                
                state.(thermal) = model.(thermal).updateEnergyConservation(state.(thermal));
                
            end
            
            %% Setup relation between E and I at positive current collectror
            
            state = model.setupEIEquation(state);
            state = model.updateControl(state, drivingForces);
            state.(ctrl) = model.(ctrl).updateControlEquation(state.(ctrl));
            
            %% Set up the governing equations

            eqs = cell(1, numel(model.equationNames));
            
            %% We collect the governing equations
            % The governing equations are the mass and charge conservation equations for the electrolyte and the
            % electrodes and the solid diffusion model equations and the control equations. The equations are scaled to
            % a common value.

            ei = model.equationIndices;

            massConsScaling = model.con.F;
            V_scaling = 1000000000;
            M_scaling = 1;
            pescaling = 1;
            sc_ne_sd = 1;
            %Vscaling = 10000000000;
            
            % Equation name : 'elyte_massCons';
            eqs{ei.elyte_massCons} = M_scaling .*state.(elyte).massCons*massConsScaling;

            % Equation name : 'elyte_chargeCons';
            eqs{ei.elyte_chargeCons} = state.(elyte).chargeCons;
            
            % Equation name : 'ne_am_chargeCons';
            eqs{ei.ne_am_chargeCons} =  state.(ne).(am).chargeCons;

            % Equation name : 'pe_am_chargeCons';
            eqs{ei.pe_am_chargeCons} = pescaling .* state.(pe).(am).chargeCons;

            %Added by Enguerran
            if battery.(pe).(am).isSwellingMaterial
                % Equation name : 'pe_am_volumeCons';
                eqs{ei.pe_am_volumeCons} = state.(pe).(am).volumeCons;
            end
            if battery.(ne).(am).isSwellingMaterial
                % Equation name : 'ne_am_volumeCons';
                eqs{ei.ne_am_volumeCons} = V_scaling * state.(ne).(am).volumeCons;
            end
            
            switch model.(ne).(am).diffusionModelType
              case 'simple'
                eqs{ei.ne_am_massCons}       = state.(ne).(am).massCons*massConsScaling;
                eqs{ei.ne_am_sd_soliddiffeq} = state.(ne).(am).(sd).solidDiffusionEq.*massConsScaling.*battery.(ne).(am).(itf).G.cells.volumes/dt;
              case 'full'
                % Equation name : 'ne_am_sd_massCons';
                n    = model.(ne).(am).(itf).n; % number of electron transfer (equal to 1 for Lithium)
                F    = model.con.F;
               
                vol  = model.(ne).(am).operators.pv;
                rp   = model.(ne).(am).(sd).rp;
                vsf  = model.(ne).(am).Interface.volumetricSurfaceArea;              

                surfp = 4.*pi.*rp.^2;
                
                scalingcoef = (vsf.*vol(1).*n.*F)./surfp;
                eqs{ei.ne_am_sd_soliddiffeq} = scalingcoef.*state.(ne).(am).(sd).solidDiffusionEq;
                
                eqs{ei.ne_am_sd_massCons}    = M_scaling .* scalingcoef.*state.(ne).(am).(sd).massCons;
                
              case 'interParticleOnly'
                eqs{ei.ne_am_massCons} = state.(ne).(am).massCons*massConsScaling;
              otherwise
                error('diffusionModelType not recognized')                    
            end
            
            
            switch model.(pe).(am).diffusionModelType
              case 'simple'
                eqs{ei.pe_am_massCons}       = state.(pe).(am).massCons*massConsScaling;
                eqs{ei.pe_am_sd_soliddiffeq} = state.(pe).(am).(sd).solidDiffusionEq.*massConsScaling.*battery.(pe).(am).(itf).G.cells.volumes/dt;
              case 'full'
                % Equation name : 'pe_am_sd_massCons';
                n    = model.(pe).(am).(itf).n; % number of electron transfer (equal to 1 for Lithium)
                F    = model.con.F;

                vol  = model.(pe).(am).operators.pv;
                rp   = model.(pe).(am).(sd).rp;
                vsf  = model.(pe).(am).(itf).volumetricSurfaceArea;

                surfp = 4.*pi.*rp.^2;
                
                scalingcoef = (vsf.*vol(1).*n.*F)./surfp;
                eqs{ei.pe_am_sd_massCons} = pescaling .*M_scaling .* scalingcoef.*state.(pe).(am).(sd).massCons;

                eqs{ei.pe_am_sd_soliddiffeq} = pescaling .*scalingcoef.*state.(pe).(am).(sd).solidDiffusionEq;
              case 'interParticleOnly'
                eqs{ei.pe_am_massCons} = state.(pe).(am).massCons*massConsScaling;                
              otherwise
                error('diffusionModelType not recognized')
            end
            
            % Equation name : 'ne_cc_chargeCons';
            if model.(ne).include_current_collectors
                eqs{ei.ne_cc_chargeCons} =state.(ne).(cc).chargeCons;
            end
            
            % Equation name : 'pe_cc_chargeCons';
            if model.(pe).include_current_collectors
                eqs{ei.pe_cc_chargeCons} = state.(pe).(cc).chargeCons;
            end

            % Equation name : 'energyCons';
            if model.use_thermal
                eqs{ei.energyCons} = state.(thermal).energyCons;
            end
            
            % Equation name : 'EIeq';
            eqs{ei.EIeq} = - state.(ctrl).EIequation;
            
            % Equation name : 'controlEq'                                    
            eqs{ei.controlEq} = state.(ctrl).controlEquation;
            
            eqs{ei.elyte_massCons} = eqs{ei.elyte_massCons} - model.Electrolyte.sp.t(1)*eqs{ei.elyte_chargeCons};
            
            names = model.equationNames;
            types = model.equationTypes;
            
            %% The equations are reordered in a way that is consitent with the linear iterative solver 
            % (the order of the equation does not matter if we do not use an iterative solver)
            ctrltype = state.Control.ctrlType;
            switch ctrltype
              case {'constantCurrent', 'CC_discharge1', 'CC_discharge2', 'CC_charge1', 'charge', 'discharge'}
                types{ei.EIeq} = 'cell';
              case {'constantVoltage', 'CV_charge2'}
                eqs([ei.EIeq, ei.controlEq]) = eqs([ei.controlEq, ei.EIeq]);
              otherwise 
                error('control type not recognized')
            end

            primaryVars = model.getPrimaryVariables();

            
            %% Setup LinearizedProblem that can be processed by MRST Newton API
            
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
            
        end
        
        

      %% cap concentrations and porosity to reasonnable values
        function [state, report] = updateState(model, state, problem, dx, drivingForces)

            [state, report] = updateState@BaseModel(model, state, problem, dx, drivingForces);
            
            % cap concentrations
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
                if model.(elde).(am).use_interparticle_diffusion
                    state.(elde).(am).c = max(cmin, state.(elde).(am).c);
                    cmax = model.(elde).(am).(itf).cmax;
                    state.(elde).(am).c = min(cmax, state.(elde).(am).c);
                else
                    state.(elde).(am).(sd).c = max(cmin, state.(elde).(am).(sd).c);
                    state.(elde).(am).(sd).cSurface = max(cmin, state.(elde).(am).(sd).cSurface);
                    cmax = model.(elde).(am).(itf).cmax;
                    state.(elde).(am).(sd).c = min(cmax, state.(elde).(am).(sd).c);
                end
            end

            %cap porosity

            if model.(ne).(am).isSwellingMaterial
                state.(ne).(am).porosity = min(1, state.(ne).(am).porosity);
                state.(ne).(am).porosity = max(0, state.(ne).(am).porosity);
            end

            ctrl = 'Control';            
            state.(ctrl) = model.(ctrl).updateControlState(state.(ctrl));
            
            report = [];
            
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
