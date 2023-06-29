classdef BatterySwelling < Battery
% 
% Used when at least one of the two electrodes is a swelling Material.
% It functions as the Battery class but implements a new equation necessary because of the particle swelling :
% the volume Conservation equation 
%
    properties
        
        primaryVarNames
        funcCallList
        
    end
    
    methods
        
        function model = BatterySwelling(paramobj)
            
            model = model@Battery(paramobj)
            
        end

        function model = registerVarAndPropfuncNames(model)
            model = registerVarAndPropfuncNames@Battery(model);

            elyte   = 'Electrolyte';
            ne      = 'NegativeElectrode';
            am      = 'ActiveMaterial';
            sd      = 'SolidDiffusion';
            itf     = 'Interface';

            varnames = {{elyte, 'volumeFraction'} ,...
                        {elyte, 'convFlux'}};

            model = model.registerVarNames(varnames);

            fn = @BatterySwelling.updateElectrolyteVolumeFraction;
            inputNames = {{ne,am, 'porosity'}};
            model = model.registerPropFunction({{elyte,'volumeFraction'}, fn, inputNames});

            fn = @BatterySwelling.updateConvFlux;
            inputNames = {{elyte, 'j'}, {elyte, 'c'}, {ne,am,sd, 'cAverage'}, {ne,am,itf, 'volumetricSurfaceArea'}};
            model = model.registerPropFunction({{elyte,'convFlux'}, fn, inputNames});

            fn = @BatterySwelling.updateAccumTerm;
            fn = {fn, @(propfunction) PropFunction.accumFuncCallSetupFn(propfunction)};
            model = model.registerPropFunction({{elyte, 'massAccum'}, fn, {{elyte,'c'},{elyte, 'volumeFraction'},{ne, am, 'porosity'}}});

        end
        
        function model = setupElectrolyte(model, paramobj)
            
            model.Electrolyte = ElectrolyteSwelling(paramobj.Electrolyte);
            
        end

        
        %% Definition of the accumulation term (dc/dt)
        function state = updateAccumTerm(model, state, state0, dt)

            elyte   = 'Electrolyte';
            ne      = 'NegativeElectrode';
            am      = 'ActiveMaterial';

            c = state.(elyte).c;
            vf = state.(elyte).volumeFraction;
            c0 = state0.(elyte).c;

            elyte_cells = zeros(model.G.cells.num, 1);
            elyte_cells(model.(elyte).G.mappings.cellmap) = (1 : model.(elyte).G.cells.num)';
            ne_cells = elyte_cells(model.(ne).(am).G.mappings.cellmap);
            vf0 = vf;
            porosity0 = state0.(ne).(am).porosity;
            vf0 = subsasgnAD(vf0, ne_cells, porosity0);
            

            cdotcc  = (vf .* c - vf0 .* c0)/dt;
            vol = model.(elyte).G.cells.volumes;

            state.(elyte).massAccum  = vol.*cdotcc;
            
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
            state.(elyte).volumeFraction = subsasgnAD(state.(elyte).volumeFraction, sep_cells, model.(elyte).(sep).porosity);


            % Define the volumeFraction in the electrodes
            eldes = {ne, pe};
            for ielde = 1 : numel(eldes)
                elde = eldes{ielde};
                if isa(model.(elde).(am), 'SwellingMaterial')
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
                if isa(model.(elde).(am), 'SwellingMaterial')

                    itf_model = model.(elde).(am).(itf);
                    G = model.(elyte).G;
                    Gp = G.mappings.parentGrid;

                    c = state.(elde).(am).SolidDiffusion.cAverage;

                    molarVolumeLithiated = model.(elde).(am).computeMolarVolumeLithiated(c);
                    densitySi            = model.(elde).(am).(itf).density;
                    molarMassSi   = model.(elde).(am).molarMass;

                    theta0 = model.(elde).(am).(itf).theta0;
                    molarVolumeSi = model.(elde).(am).constants.molarVolumeSi;
                    molarVolumeLi = model.(elde).(am).constants.molarVolumeLi;
                    molarVolumeDelithiated = (4/15)*(molarVolumeSi + 3.75*theta0*molarVolumeLi);
                    molarVolumeDelithiated = model.(elde).(am).computeMolarVolumeLithiated(c);

                    
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

                    averageVelocity = (s./(n.*F)).*(molarVolumeLithiated - molarVolumeDelithiated).*j;
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

                    if isa(model.(elde).(am), 'SwellingMaterial')
                        %Calculating the initial radius of the particle if it
                        %is swelling
                        cmax     = elde_itf.cmax;
                        theta0   = elde_itf.theta0;
                        theta100 = elde_itf.theta100;
                        R_delith = bat.(elde).(am).(sd).rp;
                        
                        
                        molarVolumeSi = bat.(elde).(am).constants.molarVolumeSi;
                        molarVolumeLi = bat.(elde).(am).constants.molarVolumeLi;
                        Q = (3.75*molarVolumeLi)/(molarVolumeSi);
    
                        theta = theta0 + SOC .* (theta100-theta0);
    
                        radius = R_delith .* (1 + Q .* theta)^(1/3);
                        
                        c = cmax .* theta .* (1+Q) .* (R_delith^3) ./ (radius^3);
                    end


                    nc    = elde_itf.G.cells.num;

                    initstate.(elde).(am).(sd).cSurface = c*ones(nc, 1);
                    N = model.(elde).(am).(sd).N;
                    np = model.(elde).(am).(sd).np; % Note : we have by construction np = nc
                    initstate.(elde).(am).(sd).c = c*ones(N*np, 1);

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
                    initstate.(ctrl).I = model.(ctrl).Imax;
                    %error('to implement (should be easy...)')
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
            if isa(model.(ne).(am), 'SwellingMaterial')
                vf = model.(ne).(am).(itf).volumeFraction;
                ADstruc = model.(ne).(am).porosity ./ model.(ne).(am).porosity;
                initstate.(ne).(am).porosity = (1 - vf) .* ADstruc;
            end

            
        end

        function control = setupControl(model, paramobj)

            C = computeCellCapacity(model, 'isSwellingMaterial', true);

            switch paramobj.controlPolicy
              case "IEswitch"
                control = IEswitchControlModel(paramobj); 
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
            
            %% We call the assembly equations ordered from the graph
            
            funcCallList = model.funcCallList;

            for ifunc = 1 : numel(funcCallList)
                eval(funcCallList{ifunc});
            end
            
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
            

            massConsScaling = model.con.F;
            
            V_scaling = 1;
            M_scaling = 1;
            pescaling = 1;
            sc_ne_sd  = 1;
            %Vscaling = 10000000000;

            eqs = {};
            names = {};
            
            % Equation name : 'elyte_massCons';
            eqs{end + 1} = M_scaling .*state.(elyte).massCons*massConsScaling;
            names{end + 1} = 'ei.elyte_massCons';

            % Equation name : 'elyte_chargeCons';
            eqs{end + 1} = state.(elyte).chargeCons;
            names{end + 1} = 'ei.elyte_chargeCons';
            
            % Equation name : 'ne_am_chargeCons';
            eqs{end + 1} =  state.(ne).(am).chargeCons;
            names{end + 1} = 'ei.ne_am_chargeCons';

            % Equation name : 'pe_am_chargeCons';
            eqs{end + 1} = pescaling .* state.(pe).(am).chargeCons;
            names{end + 1} = 'ei.pe_am_chargeCons';

            %Added by Enguerran
            if isa(battery.(pe).(am), 'SwellingMaterial')
                % Equation name : 'pe_am_volumeCons';
                eqs{end + 1} = state.(pe).(am).volumeCons;
                names{end + 1} = 'ei.pe_am_volumeCons';
            end
            if isa(battery.(ne).(am), 'SwellingMaterial')
                % Equation name : 'ne_am_volumeCons';
                eqs{end + 1} = V_scaling * state.(ne).(am).volumeCons;
                names{end + 1} = 'ei.ne_am_volumeCons';
            end

            % Equation name : 'ne_am_sd_massCons';
            n    = model.(ne).(am).(itf).n; % number of electron transfer (equal to 1 for Lithium)
            F    = model.con.F;
            
            vol  = model.(ne).(am).operators.pv;
            rp   = model.(ne).(am).(sd).rp;
            vsf  = model.(ne).(am).Interface.volumetricSurfaceArea;

            surfp = 4.*pi.*rp.^2;
            
            scalingcoef = (vsf.*vol(1).*n.*F)./surfp;
            eqs{end + 1} = scalingcoef.*state.(ne).(am).(sd).solidDiffusionEq;
            names{end + 1} = 'ei.ne_am_sd_soliddiffeq';
            
            eqs{end + 1}    = M_scaling .* scalingcoef.*state.(ne).(am).(sd).massCons;
            names{end + 1} = 'ei.ne_am_sd_massCons';
            

            % Equation name : 'pe_am_sd_massCons';
            n    = model.(pe).(am).(itf).n; % number of electron transfer (equal to 1 for Lithium)
            F    = model.con.F;

            vol  = model.(pe).(am).operators.pv;
            rp   = model.(pe).(am).(sd).rp;
            vsf  = model.(pe).(am).(itf).volumetricSurfaceArea;

            surfp = 4.*pi.*rp.^2;
            
            scalingcoef = (vsf.*vol(1).*n.*F)./surfp;
            eqs{end + 1} = pescaling .*M_scaling .* scalingcoef.*state.(pe).(am).(sd).massCons;
            names{end + 1} = 'ei.pe_am_sd_massCons';

            eqs{end + 1} = pescaling .*scalingcoef.*state.(pe).(am).(sd).solidDiffusionEq;
            names{end + 1} = 'ei.pe_am_sd_soliddiffeq';

            % Equation name : 'ne_cc_chargeCons';
            if model.(ne).include_current_collectors
                eqs{end + 1} =state.(ne).(cc).chargeCons;
                names{end + 1} = 'ei.ne_cc_chargeCons';
            end
            
            % Equation name : 'pe_cc_chargeCons';
            if model.(pe).include_current_collectors
                eqs{end + 1} = state.(pe).(cc).chargeCons;
                names{end + 1} = 'ei.pe_cc_chargeCons';
            end

            % Equation name : 'energyCons';
            if model.use_thermal
                eqs{end + 1} = state.(thermal).energyCons;
                names{end + 1} = 'ei.energyCons';
            end
            
            % Equation name : 'EIeq';
            eqs{end + 1} = - state.(ctrl).EIequation;
            names{end + 1} = 'ei.EIeq';
            
            % Equation name : 'controlEq'                                    
            eqs{end + 1} = state.(ctrl).controlEquation;
            names{end + 1} = 'ei.controlEq';
            types = repmat({'cell'}, 1, numel(names));
            
            primaryVars = model.getPrimaryVariables();


           if state.(ctrl).E > 4.1
               theta = model.NegativeElectrode.ActiveMaterial.computeTheta(state.NegativeElectrode.ActiveMaterial.SolidDiffusion.cAverage);
               theta0 = model.NegativeElectrode.ActiveMaterial.Interface.theta0;
               theta100 = model.NegativeElectrode.ActiveMaterial.Interface.theta100;

               soc = (theta.val - theta0) ./ (theta100 - theta0);
               soc = soc(1);
           end

            
            %% Setup LinearizedProblem that can be processed by MRST Newton API
            
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
            
        end
        

        %% cap concentrations and porosity to reasonnable values
        function [state, report] = updateState(model, state, problem, dx, drivingForces)

            [state, report] = updateState@Battery(model, state, problem, dx, drivingForces);
            
            ne      = 'NegativeElectrode';
            am      = 'ActiveMaterial';

            state.(ne).(am).porosity = min(1, state.(ne).(am).porosity);
            state.(ne).(am).porosity = max(0, state.(ne).(am).porosity);
            
        end
        
        function primaryvarnames = getPrimaryVariableNames(model)

            primaryvarnames = model.primaryVarNames;
            
        end
        

        function model = validateModel(model, varargin)

            model.PositiveElectrode.ActiveMaterial = model.PositiveElectrode.ActiveMaterial.setupDependentProperties();
            model.NegativeElectrode.ActiveMaterial = model.NegativeElectrode.ActiveMaterial.setupDependentProperties();
            model.Electrolyte.Separator = model.Electrolyte.Separator.setupDependentProperties();
            
            model = model.setupElectrolyteModel();

            if isempty(model.computationalGraph)
                model = model.setupComputationalGraph();
            end
            
            cgt = model.computationalGraph;
            model.primaryVarNames = cgt.getPrimaryVariableNames();
            model.funcCallList = cgt.getOrderedFunctionCallList();
            
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
