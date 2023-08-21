classdef CompositeBattery < Battery

    properties
        
        scenario

        primaryVarNames
        funcCallList
        
    end
    
    methods
        
        function model = CompositeBattery(paramobj)
            
            model = model@Battery(paramobj);

            fdnames = {'scenario'};
            model = dispatchParams(model, paramobj, fdnames);

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
            itf     = 'Interface';
            sd      = 'SolidDiffusion';
            thermal = 'ThermalModel';
            ctrl    = 'Control';

            gr = 'FirstMaterial';
            si = 'SecondMaterial';

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
                         {pe, am, sd, 'c'}            , 'pe_am_sd_massCons'      , 'cell'; ...
                         {pe, am, sd, 'cSurface'}     , 'pe_am_sd_soliddiffeq'   , 'cell' ...
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

            submodelnames = model.getSubModelNames();

            isthermal = ismember(submodelnames, 'ThermalModel');
            if any(isthermal)
                submodelnames = submodelnames(isthermal);
                model.subModelNameList = submodelnames;
            end
            
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

            gr = 'FirstMaterial';
            si = 'SecondMaterial';

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
            inputnames = {{ne, am, si, 'Rvol'}, ...
                          {ne, am, gr, 'Rvol'}, ...
                          {pe, am, 'Rvol'}};
            model = model.registerPropFunction({{elyte, 'massSource'}, fn, inputnames});
            model = model.registerPropFunction({{elyte, 'eSource'}, fn, inputnames});
            
        end

        function model = validateModel(model, varargin)

            model = validateModel@BaseModel(model, varargin{:});

            if isempty(model.funcCallList)
                model = model.setupComputationalGraph();
                cgt = model.computationalGraph;
                model.primaryVarNames = cgt.getPrimaryVariableNames();
                model.funcCallList = cgt.getOrderedFunctionCallList();
            end
            
        end

        
        function initstate = setupInitialState(model)
        % Setup the values of the primary variables at initial state

            nc       = model.G.getNumberOfCells();
            T        = model.initT;
            scenario = model.scenario;
            
            elyte   = 'Electrolyte';
            ne      = 'NegativeElectrode';
            pe      = 'PositiveElectrode';
            am      = 'ActiveMaterial';
            itf     = 'Interface';
            sd      = 'SolidDiffusion';
            cc      = 'CurrentCollector';
            thermal = 'ThermalModel';
            ctrl    = 'Control';

            gr = 'FirstMaterial';
            si = 'SecondMaterial';

            initstate.(thermal).T = T*ones(nc, 1);

            %% Synchronize temperatures
            initstate = model.updateTemperature(initstate);

            %% setup initial state for negative electrode

            mats = {gr, si};

            for imat = 1 : numel(mats)
                mat = mats{imat};
                % set primary variables
                N = model.(ne).(am).(mat).(sd).N;
                np = model.(pe).(am).(sd).np; % Note : we have by construction np = nc
                switch scenario
                  case 'discharge'
                    cElectrodeInit = (model.(ne).(am).(mat).(itf).theta100)*(model.(ne).(am).(mat).(itf).cmax);
                  case {'charge', 'first-charge'}
                    cElectrodeInit = (model.(ne).(am).(mat).(itf).theta0)*(model.(ne).(am).(mat).(itf).cmax);
                  otherwise
                    error('initCase not recognized')
                end
                
                initstate.(ne).(am).(mat).(sd).c        = cElectrodeInit*ones(N*np, 1);
                initstate.(ne).(am).(mat).(sd).cSurface = cElectrodeInit*ones(np, 1);
            end

            initstate.(ne).(am).(gr) = model.(ne).(am).(gr).updateConcentrations(initstate.(ne).(am).(gr));
            initstate.(ne).(am).(gr).(itf) = model.(ne).(am).(gr).(itf).updateOCP(initstate.(ne).(am).(gr).(itf));
            
            ref = initstate.(ne).(am).(gr).(itf).OCP(1);

            nc = model.(ne).(am).G.getNumberOfCells();
            
            initstate.(ne).(am).phi = zeros(nc, 1);
            

            %% setup initial state for positive electrode
            
            pe_itf = model.(pe).(am).(itf); 
            N = model.(pe).(am).(sd).N;
            np = model.(pe).(am).(sd).np; % Note : we have by construction np = nc

            switch scenario
              case 'discharge'
                cElectrodeInit = (model.(pe).(am).(itf).theta100)*(model.(pe).(am).(itf).cmax);
              case {'charge', 'first-charge'}
                cElectrodeInit = (model.(pe).(am).(itf).theta0)*(model.(pe).(am).(itf).cmax);
              otherwise
                error('initCase not recognized')
            end
            initstate.(pe).(am).(sd).cSurface = cElectrodeInit*ones(np, 1);
            np = model.(pe).(am).(sd).np; % Note : we have by construction np = nc
            initstate.(pe).(am).(sd).c = cElectrodeInit*ones(N*np, 1);

            initstate.(pe).(am) = model.(pe).(am).updateConcentrations(initstate.(pe).(am));
            initstate.(pe).(am).(itf) = pe_itf.updateOCP(initstate.(pe).(am).(itf));

            initstate.(pe).(am).phi = initstate.(pe).(am).(itf).OCP - ref;
            
            %% Setup initial Electrolyte state

            initstate.(elyte).phi = zeros(model.(elyte).G.getNumberOfCells(), 1) - ref;
            initstate.(elyte).c   = 1000*ones(model.(elyte).G.getNumberOfCells(), 1);

            initstate.(ctrl).E = initstate.(pe).(am).phi(1) - initstate.(ne).(am).phi(1);
            switch scenario
              case 'discharge'
                initstate.(ctrl).I = - model.(ctrl).Imax;
              case {'charge', 'first-charge'}
                initstate.(ctrl).I = model.(ctrl).Imax;
              otherwise
                error('scenario not recognized')
            end 
            
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

            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            am    = 'ActiveMaterial';
            sd    = 'SolidDiffusion';
            itf   = 'Interface';
            ctrl  = 'Control';
            gr    = 'FirstMaterial';
            si    = 'SecondMaterial';

            funcCallList = model.funcCallList;
            
            for ifunc = 1 : numel(funcCallList)
                eval(funcCallList{ifunc});
            end

            eqs = cell(1, numel(model.equationNames));
            ei = model.equationIndices;

            massConsScaling = model.con.F;
            
            eqs{ei.elyte_massCons}   = state.(elyte).massCons*massConsScaling;
            eqs{ei.elyte_chargeCons} = state.(elyte).chargeCons;
            eqs{ei.ne_am_chargeCons} = state.(ne).(am).chargeCons;
            eqs{ei.pe_am_chargeCons} = state.(pe).(am).chargeCons;
            eqs{ei.EIeq}             = state.(ctrl).EIequation;
            eqs{ei.controlEq}        = state.(ctrl).controlEquation;

            n    = model.(ne).(am).(gr).(itf).n; % number of electron transfer (equal to 1 for Lithium)
            F    = model.con.F;
            vol  = model.(ne).(am).(gr).G.getVolumes();
            rp   = model.(ne).(am).(gr).(sd).rp;
            vsf  = model.(ne).(am).(gr).(sd).volumetricSurfaceArea;
            surfp = 4*pi*rp^2;
            
            scalingcoef = (vsf*vol(1)*n*F)/surfp;
            
            eqs{ei.ne_am_gr_sd_massCons}    = scalingcoef*state.(ne).(am).(gr).(sd).massCons;
            eqs{ei.ne_am_gr_sd_soliddiffeq} = scalingcoef*state.(ne).(am).(gr).(sd).solidDiffusionEq;

            n    = model.(ne).(am).(si).(itf).n; % number of electron transfer (equal to 1 for Lithium)
            F    = model.con.F;
            vol  = model.(ne).(am).(si).G.getVolumes();
            rp   = model.(ne).(am).(si).(sd).rp;
            vsf  = model.(ne).(am).(si).(sd).volumetricSurfaceArea;
            surfp = 4*pi*rp^2;
            
            scalingcoef = (vsf*vol(1)*n*F)/surfp;
                        
            eqs{ei.ne_am_si_sd_massCons}    = scalingcoef*state.(ne).(am).(si).(sd).massCons;        
            eqs{ei.ne_am_si_sd_soliddiffeq} = scalingcoef*state.(ne).(am).(si).(sd).solidDiffusionEq;

            n    = model.(pe).(am).(itf).n; % number of electron transfer (equal to 1 for Lithium)
            F    = model.con.F;
            vol  = model.(pe).(am).G.getVolumes();
            rp   = model.(pe).(am).(sd).rp;
            vsf  = model.(pe).(am).(sd).volumetricSurfaceArea;
            surfp = 4*pi*rp^2;
            
            scalingcoef = (vsf*vol(1)*n*F)/surfp;

            eqs{ei.pe_am_sd_massCons}    = scalingcoef*state.(pe).(am).(sd).massCons;        
            eqs{ei.pe_am_sd_soliddiffeq} = scalingcoef*state.(pe).(am).(sd).solidDiffusionEq;

            names = model.equationNames;
            types = model.equationTypes;
            
            primaryVars = model.getPrimaryVariableNames();
            
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

            gr = 'FirstMaterial';
            si = 'SecondMaterial';
            
            % (here we assume that the ThermalModel has the "parent" grid)
            state.(elyte).T   = state.(thermal).T(model.(elyte).G.mappings.cellmap);
            state.(ne).(am).T = state.(thermal).T(model.(ne).(am).G.mappings.cellmap);
            state.(pe).(am).T = state.(thermal).T(model.(pe).(am).G.mappings.cellmap);
            
            % Update temperature in the active materials of the electrodes.
            state.(ne).(am)      = model.(ne).(am).dispatchTemperature(state.(ne).(am));
            state.(ne).(am).(si) = model.(ne).(am).(si).dispatchTemperature(state.(ne).(am).(si));
            state.(ne).(am).(gr) = model.(ne).(am).(gr).dispatchTemperature(state.(ne).(am).(gr));
            state.(pe).(am)      = model.(pe).(am).dispatchTemperature(state.(pe).(am));
            
        end
        
        function state = updateElectrolyteCoupling(model, state)
        % Assemble the electrolyte coupling by adding the ion sources from the electrodes
            
            battery = model;

            elyte = 'Electrolyte';
            ne    = 'NegativeElectrode';
            pe    = 'PositiveElectrode';
            am    = 'ActiveMaterial';
            itf   = 'Interface';
            
            gr = 'FirstMaterial';
            si = 'SecondMaterial';

            vols = battery.(elyte).G.getVolumes();
            F = battery.con.F;
            
            couplingterms = battery.couplingTerms;

            elyte_c_source = zeros(battery.(elyte).G.getNumberOfCells(), 1);
            elyte_e_source = zeros(battery.(elyte).G.getNumberOfCells(), 1);
            
            % setup AD 
            phi = state.(elyte).phi;
            if isa(phi, 'ADI')
                adsample = getSampleAD(phi);
                adbackend = model.AutoDiffBackend;
                elyte_c_source = adbackend.convertToAD(elyte_c_source, adsample);
            end
            
            coupnames = model.couplingNames;
            
            ne_si_Rvol = state.(ne).(am).(si).Rvol;
            ne_gr_Rvol = state.(ne).(am).(gr).Rvol;
            if (isa(ne_si_Rvol, 'ADI') | isa(ne_gr_Rvol, 'ADI')) & ~isa(elyte_c_source, 'ADI')
                adsample = getSampleAD(ne_Rvol);
                adbackend = model.AutoDiffBackend;
                elyte_c_source = adbackend.convertToAD(elyte_c_source, adsample);
            end
            
            coupterm = getCoupTerm(couplingterms, 'NegativeElectrode-Electrolyte', coupnames);
            elytecells = coupterm.couplingcells(:, 2);
            elyte_c_source(elytecells) = (ne_si_Rvol + ne_gr_Rvol).*vols(elytecells);
            
            pe_Rvol = state.(pe).(am).Rvol;
            if isa(pe_Rvol, 'ADI') & ~isa(elyte_c_source, 'ADI')
                adsample = getSampleAD(pe_Rvol);
                adbackend = model.AutoDiffBackend;
                elyte_c_source = adbackend.convertToAD(elyte_c_source, adsample);
            end
            
            coupterm = getCoupTerm(couplingterms, 'PositiveElectrode-Electrolyte', coupnames);
            elytecells = coupterm.couplingcells(:, 2);
            elyte_c_source(elytecells) = pe_Rvol.*vols(elytecells);
            
            elyte_e_source = elyte_c_source.*battery.(elyte).sp.z(1)*F; 
            
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

            gr = 'FirstMaterial';
            si = 'SecondMaterial';
            
            phi_elyte = state.(elyte).phi;
            c_elyte = state.(elyte).c;
            
            elyte_cells = zeros(model.G.getNumberOfCells(), 1);
            elyte_cells(bat.(elyte).G.mappings.cellmap) = (1 : bat.(elyte).G.getNumberOfCells())';

            state.(ne).(am).(si).(itf).phiElectrolyte = phi_elyte(elyte_cells(bat.(ne).(am).(si).G.mappings.cellmap));
            state.(ne).(am).(si).(itf).cElectrolyte   = c_elyte(elyte_cells(bat.(ne).(am).(si).G.mappings.cellmap));
            state.(ne).(am).(gr).(itf).phiElectrolyte = phi_elyte(elyte_cells(bat.(ne).(am).(gr).G.mappings.cellmap));
            state.(ne).(am).(gr).(itf).cElectrolyte   = c_elyte(elyte_cells(bat.(ne).(am).(gr).G.mappings.cellmap));
            
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
