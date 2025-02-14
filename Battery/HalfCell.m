classdef HalfCell < BaseModel
% 
% Used when at least one of the two electrodes is a swelling Material.
% It functions as the Battery class but implements a new equation necessary because of the particle swelling :
% the volume Conservation equation 
%
    properties

        Electrolyte
        Separator
        ActiveMaterial

        Control
        SOC
        use_thermal

        externalCoupling
        
        initT
        
        primaryVarNames
        funcCallList
        
    end
    
    methods
        
        function model = HalfCell(inputparams)
            
            model = model@BaseModel()

            model.Electrolyte    = Electrolyte(inputparams.Electrolyte);
            model.ActiveMaterial = ActiveMaterial(inputparams.ActiveMaterial);
            
            model.Control        = IEswitchControlModel(inputparams.Control);

            model.initT = 298;
            
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
            model.Control = ControlModel([]);
            model.Control.controlPolicy = 'None';
            drivingForces = model.getValidDrivingForces();

            % We call getEquations to update state
            
            [~, state] = getEquations(model, state0, state, dt, drivingForces, 'ResOnly', true);

            % we set to empty the fields we know that are not meaningfull

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


            switch model.(am).diffusionModelType
              case 'full'
                state.(am).(sd).massAccum = [];
                state.(am).(sd).massCons = [];
              case {'simple', 'interParticleOnly'}
                state.(am).massAccum = [];
                state.(am).massCons = [];                    
              otherwise
                error('diffusion model type not recognized');
            end
           
           
        
            switch model.(am).diffusionModelType
              case 'full'
                state.(elde).(sd) = model.(am).(sd).updateAverageConcentration(state.(elde).(am).(sd));
                state.(elde).(am) = model.(am).updateSOC(state.(elde).(am));
                state.(elde).(am) = model.(am).updateAverageConcentration(state.(elde).(am));
              case {'simple', 'interParticleOnly'}
                % do nothing
              otherwise
                error('diffusion model type not recognized');
            end
            

            if model.use_thermal
                state = model.updateThermalIrreversibleReactionSourceTerms(state);
                state = model.updateThermalReversibleReactionSourceTerms(state);
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
        
        function [problem, state] = getEquations(model, state0, state, dt, drivingForces, varargin)
        % Assembly of the governing equation
            
            opts = struct('ResOnly', false, 'iteration', 0, 'reverseMode', false); 
            opts = merge_options(opts, varargin{:});
            
            time = state0.time + dt;
            %if(not(opts.ResOnly) && not(opts.reverseMode))
            %    state = model.initStateAD(state);
            %    state.Electrolyte.c   = model.Electrolyte.cInit;
            %    state.Electrolyte.phi = model.Electrolyte.phiInit;
            %elseif(opts.reverseMode)
            %   disp('No AD initatlization in equation old style')
            %   state0 = model.initStateAD(state0);
            %else
            %    assert(opts.ResOnly);
            %end


            %% We call the assembly equations ordered from the graph
            %state = model.updateElectrolyte();

            funcCallList = model.funcCallList;

            for ifunc = 1 : numel(funcCallList)
                eval(funcCallList{ifunc});
            end

            eqs = {};
            names = {};
            
            eqs{end + 1} = state.(am).chargeCons
            names{end + 1} = 'am_chargeCons';

            eqs{end + 1} = state.(am).(sd).solidDiffusionEq;
            names{end + 1} = 'am_sd_solidDiffusionEq';
            
            eqs{end + 1} = state.(am).(sd).massCons;
            names{end + 1} = 'am_sd_massCons';
            
            eqs{end + 1} = state.(ctrl).controlEquation;
            names{end + 1} = 'ctrl_controlEquation';
            
            eqs{end + 1} = state.(ctrl).EIequation;
            names{end + 1} = 'ctrl_EIequation';
            
            types = repmat({'cell'}, 1, numel(names));
            
            primaryVars = model.getPrimaryVariables();
            
            %% Setup LinearizedProblem that can be processed by MRST Newton API
            
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);

        end

        function primaryvarnames = getPrimaryVariableNames(model)

            primaryvarnames = model.primaryVarNames;
            
        end
     
        function model = validateModel(model, varargin)

            model.ActiveMaterial = model.ActiveMaterial.setupDependentProperties();
            
            model = model.setupElectrolyteModel();

            if isempty(model.computationalGraph)
                model = model.setupComputationalGraph();
            end
            
            cgt = model.computationalGraph;
            model.primaryVarNames = cgt.getPrimaryVariableNames();
            model.funcCallList = cgt.getOrderedFunctionCallList();
            
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
            ctrl    = 'Control';

            if ~model.use_thermal
                % we register the temperature variable, as it is not done by the ThermalModel which is empty in this case
                model = model.registerVarName({thermal, 'T'});
            end
            
            varnames = {{ne, am, 'E'}, ...
                        {ne, am, 'I'}};
            model = model.removeVarNames(varnames);

            
            %% Temperature dispatch functions
            fn = @Battery.updateTemperature;
            
            inputnames = {{thermal, 'T'}};
            model = model.registerPropFunction({{elyte, 'T'}  , fn, inputnames});
            model = model.registerPropFunction({{am, 'T'} , fn, inputnames});

                  
            %% Coupling functions
            
            % Dispatch electrolyte concentration and potential in the electrodes
            fn = @Battery.updateElectrodeCoupling;
            inputnames = {{elyte, 'c'}, ...
                          {elyte, 'phi'}};
            model = model.registerPropFunction({{am, itf, 'phiElectrolyte'}, fn, inputnames});
            model = model.registerPropFunction({{am, itf, 'cElectrolyte'}  , fn, inputnames});

            % Functions that update the source terms in the electolyte
            fn = @Battery.updateElectrolyteCoupling;
            inputnames = {{am, 'Rvol'}};
            model = model.registerPropFunction({{elyte, 'massSource'}, fn, inputnames});
            model = model.registerPropFunction({{elyte, 'eSource'}, fn, inputnames});
            
            % Function that assemble the control equation
            fn = @Battery.setupEIEquation;
            inputnames = {{ctrl, 'E'}, ...
                          {ctrl, 'I'}}; 
            
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
                              {am, 'jFace'}       , ...
                              {elyte, 'conductivity'} , ...
                              {am, 'conductivity'}};

                               
                model = model.registerPropFunction({{thermal, 'jHeatOhmSource'}, fn, inputnames});
                
                %% Function that updates the Thermal Chemical Terms
                fn = @Battery.updateThermalChemicalSourceTerms;
                inputnames = {{elyte, 'diffFlux'}, ...
                              {elyte, 'D'}       , ...
                              VarName({elyte}, 'dmudcs', 2)};
                model = model.registerPropFunction({{thermal, 'jHeatChemicalSource'}, fn, inputnames});
                
                %% Function that updates Thermal Reaction Terms
                fn = @Battery.updateThermalReactionSourceTerms;
                inputnames = {{thermal, 'T'}       , ...
                              {am, 'Rvol'}     , ...
                              {am, itf, 'eta'} , ...
                              {am, itf, 'dUdT'}};
                model = model.registerPropFunction({{thermal, 'jHeatReactionSource'}, fn, inputnames});

            else
                model = model.removeVarName({elyte, 'diffFlux'});
                model = model.removeVarName({am, itf, 'dUdT'});
            end
            
            %% Functions that setup external  coupling for negative electrode
            

            fn = @Battery.setupExternalCouplingNegativeElectrode;

            inputnames = {{am, 'phi'}, ...
                          {am, 'conductivity'}};
            model = model.registerPropFunction({{am, 'jExternal'}, fn, inputnames});
            if model.use_thermal
                model = model.registerPropFunction({{am, 'jFaceExternal'}, fn, inputnames});
            end
                             

            %% Declare the "static" variables
            varnames = {};
            if ~model.use_thermal
                varnames{end + 1} = {thermal, 'T'};
            end
            model = model.registerStaticVarNames(varnames);

            %% Declare the variables used for thermal model
            eldes = {ne, pe};

            if ~model.use_thermal

                    varnames = {{am, 'jFace'}, ...
                                {am, 'jFaceCoupling'}, ...
                                {am, 'jFaceBc'}};
                    model = model.removeVarNames(varnames);

            end
            
        end

        function state = updateTemperature(model, state)

            am    = 'ActiveMaterial';
            thermal = 'ThermalModel';
            
            state.(am).T = state.(thermal).T(model.(am).G.mappings.cellmap);
            state.(am) = model.(am).dispatchTemperature(state.(am));
            
        end

        function state = setupEIEquation(model, state)
            
            ctrl = 'Control';
            am = 'ActiveMaterial';
            
            I = state.(ctrl).I;
            E = state.(ctrl).E;
            
            phi = state.(am).phi;
            
            coupterm = model.(am).externalCouplingTerm;
            faces    = coupterm.couplingfaces;
            cond_pcc = model.(am).EffectiveElectricalConductivity;
            [trans_pcc, cells] = model.(am).operators.harmFaceBC(cond_pcc, faces);
            
            state.Control.EIequation = sum(trans_pcc.*(state.(am).phi(cells) - E)) - I;

        end

        function state = updateControl(model, state, drivingForces)
            
            ctrl = "Control";
            
            E    = state.(ctrl).E;
            I    = state.(ctrl).I;
            time = state.time;
            
            [ctrlVal, ctrltype] = drivingForces.src(time, value(I), value(E));
            
            state.(ctrl).ctrlVal  = ctrlVal;
            state.(ctrl).ctrlType = ctrltype;
            
        end

        function initstate = setupInitialState(model)
        % Setup the values of the primary variables at initial state

            nc = model.ActiveMaterial.G.cells.num;

            SOC = model.SOC;
            T   = model.initT;
            
            halfCell = model;
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

            %% Setup initial state for Electrode

             elde_itf = halfCell.(am).(itf); 

             theta = SOC*(elde_itf.theta100 - elde_itf.theta0) + elde_itf.theta0;
             c     = theta*elde_itf.cmax;
             nc    = elde_itf.G.cells.num;
             switch model.(am).diffusionModelType
               case 'simple'
                 initstate.(am).(sd).cSurface = c*ones(nc, 1);
                 initstate.(am).c = c*ones(nc, 1);
               case 'full'
                 initstate.(am).(sd).cSurface = c*ones(nc, 1);
                 N = model.(am).(sd).N;
                 np = model.(am).(sd).np; % Note : we have by construction np = nc
                 initstate.(am).(sd).c = c*ones(N*np, 1);
               case 'interParticleOnly'
                 initstate.(am).c = c*ones(nc, 1);                    
               otherwise
                 error('diffusionModelType not recognized')
             end
             
             initstate.(am) = model.(am).updateConcentrations(initstate.(am));
             initstate.(am).(itf) = elde_itf.updateOCP(initstate.(am).(itf));
             OCP = initstate.(am).(itf).OCP;
                 
             % The value in the first cell is used as reference.
             ref = OCP(1);
             initstate.(am).phi = OCP - ref;

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
            
            ne_Rvol = state.(ne).(am).Rvol;
            if isa(ne_Rvol, 'ADI') & ~isa(elyte_c_source, 'ADI')
                adsample = getSampleAD(ne_Rvol);
                adbackend = model.AutoDiffBackend;
                elyte_c_source = adbackend.convertToAD(elyte_c_source, adsample);
            end
            
            coupterm = getCoupTerm(couplingterms, 'NegativeElectrode-Electrolyte', coupnames);
            elytecells = coupterm.couplingcells(:, 2);
            elyte_c_source(elytecells) = ne_Rvol.*vols(elytecells);
            
            pe_Rvol = state.(pe).(am).Rvol;
            if isa(pe_Rvol, 'ADI') & ~isa(elyte_c_source, 'ADI')
                adsample = getSampleAD(pe_Rvol);
                adbackend = model.AutoDiffBackend;
                elyte_c_source = adbackend.convertToAD(elyte_c_source, adsample);
            end
            
            coupterm = getCoupTerm(couplingterms, 'PositiveElectrode-Electrolyte', coupnames);
            elytecells = coupterm.couplingcells(:, 2);
            elyte_c_source(elytecells) = - ne_Rvol.*vols(elytecells);
            
            elyte_e_source = elyte_c_source.*battery.(elyte).sp.z(1)*F; 
            
            state.Electrolyte.massSource = elyte_c_source; 
            state.Electrolyte.eSource = elyte_e_source;
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
            am    = 'ActiveMaterial';
            itf   = 'Interface';
            cc    = 'CurrentCollector';

            phi_elyte = state.(elyte).phi;
            c_elyte = state.(elyte).c;
            
            elyte_cells = zeros(model.G.cells.num, 1);
            elyte_cells(bat.(elyte).G.mappings.cellmap) = (1 : bat.(elyte).G.cells.num)';

           state.(am).(itf).phiElectrolyte = phi_elyte(elyte_cells(bat.(am).G.mappings.cellmap));
           state.(am).(itf).cElectrolyte = c_elyte(elyte_cells(bat.(am).G.mappings.cellmap));

            
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
