classdef SiliconGraphiteBattery < Battery

    methods
        
        function model = SiliconGraphiteBattery(paramobj)
            
            model = model@Battery(paramobj);

            assert(model.use_thermal == false, 'thermal model not covered yet');
            assert(model.include_current_collectors == false, 'current collectors  not included yet');
            
        end


        function model = setupSelectedModel(model);
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
            ctrl    = 'Control';

            gr = 'Graphite';
            si = 'Silicon';

            varEqTypes ={{elyte, 'c'}                 , 'elyte_massCons'         , 'cell'; ...  
                         {elyte, 'phi'}               , 'elyte_chargeCons'       , 'cell'; ...    
                         {ne, am, 'phi'}              , 'ne_am_chargeCons'       , 'cell'; ...    
                         {pe, am, 'phi'}              , 'pe_am_chargeCons'       , 'cell'; ...    
                         {ctrl, 'E'}                  , 'EIeq'                   , 'ctrl'; ...  
                         {ctrl, 'I'}                  , 'controlEq'              , 'ctrl'; ...
                         {ne, am, gr, sd, 'c'}        , 'ne_am_gr_sd_massCons'   , 'cell'; ...
                         {ne, am, gr, sd, 'cSurface'} , 'ne_am_gr_sd_soliddiffeq', 'cell'; ...
                         {ne, am, si, sd, 'c'}        , 'ne_am_si_sd_massCons'   , 'cell'; ...
                         {ne, am, si, sd, 'cSurface'} , 'ne_am_si_sd_soliddiffeq', 'cell'; ...
                         {pe, am, gr, sd, 'c'}        , 'pe_am_gr_sd_massCons'   , 'cell'; ...
                         {pe, am, gr, sd, 'cSurface'} , 'pe_am_gr_sd_soliddiffeq', 'cell'; ...
                         {pe, am, si, sd, 'c'}        , 'pe_am_si_sd_massCons'   , 'cell'; ...
                         {pe, am, si, sd, 'cSurface'} , 'pe_am_si_sd_soliddiffeq', 'cell'; ...
                         };

            addedVariableNames = {};
            addedVariableNames{end + 1} = {ctrl, 'ctrlType'};
            
            if strcmp(model.(ctrl).controlPolicy, 'CCCV')
                addedVariableNames{end + 1} = {ctrl, 'nextCtrlType'};
            end

            primaryVariableNames = varEqTypes(:, 1);
            equationNames        = varEqTypes(:, 2);
            equationTypes        = varEqTypes(:, 3);

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
        
        function model = registerVarAndPropfuncNames(model)
            
        %% Declaration of the Dynamical Variables and Function of the model
            % (setup of varnameList and propertyFunctionList)
            
            model = registerVarAndPropfuncNames@Battery(model);
            
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
            ctrl    = 'Control';

            gr = 'Graphite';
            si = 'Silicon';

            %%  removing some unused variables
            varnames = {{ne, am, itf, 'phiElectrolyte'},
                        {ne, am, itf, 'cElectrolyte'}};

            model = model.removeVarNames(varnames);
            
            %% Coupling functions
            
            % Dispatch electrolyte concentration and potential in the electrodes
            fn = @Battery.updateElectrodeCoupling;
            inputnames = {{elyte, 'c'}, ...
                          {elyte, 'phi'}};
            model = model.registerPropFunction({{ne, am, gr, itf, 'phiElectrolyte'}, fn, inputnames});
            model = model.registerPropFunction({{ne, am, si, itf, 'phiElectrolyte'}, fn, inputnames});
            model = model.registerPropFunction({{ne, am, gr, itf, 'cElectrolyte'}  , fn, inputnames});
            model = model.registerPropFunction({{ne, am, si, itf, 'cElectrolyte'}  , fn, inputnames});
            model = model.registerPropFunction({{pe, am, itf, 'phiElectrolyte'}, fn, inputnames});
            model = model.registerPropFunction({{pe, am, itf, 'cElectrolyte'}  , fn, inputnames});
            
            % Functions that update the source terms in the electolyte
            fn = @Battery.updateElectrolyteCoupling;
            inputnames = {{ne, am, 'eSource'}, ...
                          {pe, am, 'eSource'}};
            model = model.registerPropFunction({{elyte, 'massSource'}, fn, inputnames});
            model = model.registerPropFunction({{elyte, 'eSource'}, fn, inputnames});
            
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

                theta = SOC*(elde_itf.theta100 - elde_itf.theta0) + elde_itf.theta0;
                c     = theta*elde_itf.cmax;
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

            if model.(ne).include_current_collector
                OCP = initstate.(ne).(am).(itf).OCP;
                OCP = OCP(1) .* ones(bat.(ne).(cc).G.cells.num, 1);
                initstate.(ne).(cc).phi = OCP - ref;
            end
            
            if model.(pe).include_current_collector
                OCP = initstate.(pe).(am).(itf).OCP;
                OCP = OCP(1) .* ones(bat.(pe).(cc).G.cells.num, 1);
                initstate.(pe).(cc).phi = OCP - ref;
            end
            
            initstate.(ctrl).E = OCP(1) - ref;
            initstate.(ctrl).I = - model.(ctrl).Imax;
            
            switch model.(ctrl).controlPolicy
              case 'CCCV'
                initstate.(ctrl).ctrlType = 'CC_charge1';
                initstate.(ctrl).nextCtrlType = 'CC_charge1';
              case 'IEswitch'
                initstate.(ctrl).ctrlType = 'constantCurrent';
              otherwise
                error('control policy not recognized');
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
               disp('No AD initatlization in equation old style')
               state0 = model.initStateAD(state0);
            else
                assert(opts.ResOnly);
            end
            
           
            names = model.equationNames;
            types = model.equationTypes;
            
            %% The equations are reordered in a way that is consitent with the linear iterative solver 
            % (the order of the equation does not matter if we do not use an iterative solver)
            ctrltype = state.Control.ctrlType;
            switch ctrltype
              case {'constantCurrent', 'CC_discharge1', 'CC_discharge2', 'CC_charge1'}
                types{ei.EIeq} = 'cell';   
              case {'constantVoltage', 'CV_charge2'}
                % no changes
              otherwise 
                error('control type not recognized')
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
            
            F = battery.con.F;
            
            couplingterms = battery.couplingTerms;

            elyte_e_source = zeros(battery.(elyte).G.cells.num, 1);
            elyte_c_source = zeros(battery.(elyte).G.cells.num, 1);
            
            % setup AD 
            phi = state.(elyte).phi;
            if isa(phi, 'ADI')
                adsample = getSampleAD(phi);
                adbackend = model.AutoDiffBackend;
                elyte_e_source = adbackend.convertToAD(elyte_e_source, adsample);
            end
            
            coupnames = model.couplingNames;

            ne_e_source = state.(ne).(am).eSource;
            if isa(ne_e_source, 'ADI') & ~isa(elyte_e_source, 'ADI')
                adsample = getSampleAD(ne_e_source);
                adbackend = model.AutoDiffBackend;
                elyte_e_source = adbackend.convertToAD(elyte_e_source, adsample);
            end
            
            coupterm = getCoupTerm(couplingterms, 'NegativeElectrode-Electrolyte', coupnames);
            elytecells = coupterm.couplingcells(:, 2);
            elyte_e_source(elytecells) = - ne_e_source;
            
            pe_e_source = state.(pe).(am).Rvol;
            if isa(pe_e_source, 'ADI') & ~isa(elyte_e_source, 'ADI')
                adsample = getSampleAD(pe_e_source);
                adbackend = model.AutoDiffBackend;
                elyte_e_source = adbackend.convertToAD(elyte_e_source, adsample);
            end
            
            coupterm = getCoupTerm(couplingterms, 'PositiveElectrode-Electrolyte', coupnames);
            elytecells = coupterm.couplingcells(:, 2);
            elyte_e_source(elytecells) = - pe_e_source;
            
            elyte_c_source = elyte_e_source./(battery.(elyte).sp.z(1)*F); 
            
            state.Electrolyte.massSource = elyte_c_source; 
            state.Electrolyte.eSource = elyte_e_source;
            
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

              case 'None'
                % nothing done here. This case is only used for addVariables function
              otherwise
                error('control type not recognized');
            end
                
            
        end
        
        function state = updateElectrodeCoupling(model, state)
        % Setup the electrode coupling by updating the potential and concentration of the electrolyte in the active
        % component of the electrodes. There, those quantities are considered as input and used to compute the reaction
        % rate.

            bat = model;
            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            am    = 'ActiveMaterial';
            itf   = 'Interface';

            gr = 'Graphite';
            si = 'Silicon';
            
            phi_elyte = state.(elyte).phi;
            c_elyte = state.(elyte).cs{1};
            
            elyte_cells = zeros(model.G.cells.num, 1);
            elyte_cells(bat.(elyte).G.mappings.cellmap) = (1 : bat.(elyte).G.cells.num)';

            state.(ne).(am).(si).(itf).phiElectrolyte = phi_elyte(elyte_cells(bat.(ne).(am).(si).G.mappings.cellmap));
            state.(ne).(am).(si).(itf).cElectrolyte   = c_elyte(elyte_cells(bat.(ne).(am).(si).G.mappings.cellmap));
            state.(ne).(am).(gr).(itf).phiElectrolyte = phi_elyte(elyte_cells(bat.(ne).(am).(gr).G.mappings.cellmap));
            state.(ne).(am).(gr).(itf).cElectrolyte   = c_elyte(elyte_cells(bat.(ne).(am).(gr).G.mappings.cellmap));s
            
            state.(pe).(am).(itf).phiElectrolyte = phi_elyte(elyte_cells(bat.(pe).(am).G.mappings.cellmap));
            state.(pe).(am).(itf).cElectrolyte = c_elyte(elyte_cells(bat.(pe).(am).G.mappings.cellmap));
            
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
                if strcmp(model.(elde).(am).diffusionModelType, 'simple') | ~model.use_particle_diffusion
                    state.(elde).(am).c = max(cmin, state.(elde).(am).c);
                    cmax = model.(elde).(am).(itf).cmax;
                    state.(elde).(am).c = min(cmax, state.(elde).(am).c);
                else
                    state.(elde).(am).(sd).c = max(cmin, state.(elde).(am).(sd).c);
                    cmax = model.(elde).(am).(itf).cmax;
                    state.(elde).(am).(sd).c = min(cmax, state.(elde).(am).(sd).c);
                end
            end
            
            ctrl = 'Control';            
            state.(ctrl) = model.(ctrl).updateControlState(state.(ctrl));
            
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

        function model = setupCapping(model)
            model.cmin = 0;
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
