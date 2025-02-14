classdef BatterySwelling < GenericBattery
% 
    methods
        
        function model = BatterySwelling(inputparams)
            
            model = model@GenericBattery(inputparams)
            
        end

        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@GenericBattery(model);

            elyte   = 'Electrolyte';
            ne      = 'NegativeElectrode';
            co      = 'Coating';
            am      = 'ActiveMaterial';
            sd      = 'SolidDiffusion';
            itf     = 'Interface';

            fn = @BatterySwelling.updateElectrolyteVolumeFraction;
            inputNames = {{ne, co, am, 'porosity'}};
            model = model.registerPropFunction({{elyte,'volumeFraction'}, fn, inputNames});

            fn = @BatterySwelling.updateConvFlux;
            inputNames = {{elyte, 'j'}, {elyte, 'c'}, {ne, co, am, sd, 'cAverage'}, {ne, co, am, itf, 'volumetricSurfaceArea'}};
            model = model.registerPropFunction({{elyte,'convFlux'}, fn, inputNames});

            fn = @BatterySwelling.updateAccumTerm;
            fn = {fn, @(propfunction) PropFunction.accumFuncCallSetupFn(propfunction)};
            inputNames = {{elyte,'c'}, {elyte, 'volumeFraction'}, {ne, co, am, 'porosity'}};
            model = model.registerPropFunction({{elyte, 'massAccum'}, fn, inputNames});

        end
        
        function model = setupElectrolyteModel(model, inputparams)
            
            model.Electrolyte = ElectrolyteSwelling(inputparams.Electrolyte);
            
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
            porosity0 = state0.(ne).(co).(am).porosity;
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
                    state.(elyte).volumeFraction = subsasgnAD(state.(elyte).volumeFraction, elyte_cells(model.(elde).(am).G.mappings.cellmap), state.(elde).(co).(am).porosity);
                else
                    state.(elyte).volumeFraction = subsasgnAD(state.(elyte).volumeFraction, elyte_cells(model.(elde).(am).G.mappings.cellmap), model.(elde).(co).(am).porosity);
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

                    G              = model.(elyte).G;
                    Gp             = G.mappings.parentGrid;

                    cmax           = model.(elde).(am).(itf).cmax;
                    theta0         = model.(elde).(am).(itf).theta0;
                    densitySi      = model.(elde).(am).(itf).density;
                    molarMassSi    = model.(elde).(am).molarMass;
                    F              = model.(elde).(am).(itf).constants.F;
                    s = -1;
                    n = 1;

                    c = state.(elde).(am).SolidDiffusion.cAverage;

                    theta = c./cmax;
                    
                    molarVolumeLithiated = model.(elde).(am).computeMolarVolumeLithiated(theta);
                    molarVolumeDelithiated = model.(elde).(am).computeMolarVolumeLithiated(theta0);

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
        function initstate = setupInitialState(model, jsonstruct)
            
            if isa(model.(ne).(co).(am), 'SwellingMaterial')
                vf = model.(ne).(co).(am).(itf).volumeFraction;
                ADstruc = model.(ne).(co).(am).porosity ./ model.(ne).(co).(am).porosity;
                initstate.(ne).(co).(am).porosity = (1 - vf) .* ADstruc;
            end

            
        end

        function control = setupControl(model, inputparams)

            % C = computeCellCapacity(model, 'isSwellingMaterial', true);
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


        %% Assembly of the governing equation (same as for Battery but taking into account porosity variations)
        function [problem, state] = setupScalingFix(model, state0, state, dt, drivingForces, varargin)
            

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
            eqs{end + 1} =  state.(ne).(co).(am).chargeCons;
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
            if isa(battery.(ne).(co).(am), 'SwellingMaterial')
                % Equation name : 'ne_am_volumeCons';
                eqs{end + 1} = V_scaling * state.(ne).(co).(am).volumeCons;
                names{end + 1} = 'ei.ne_am_volumeCons';
            end

            % Equation name : 'ne_am_sd_massCons';
            n    = model.(ne).(co).(am).(itf).n; % number of electron transfer (equal to 1 for Lithium)
            F    = model.con.F;
            
            vol  = model.(ne).(co).(am).operators.pv;
            rp   = model.(ne).(co).(am).(sd).rp;
            vsf  = model.(ne).(co).(am).Interface.volumetricSurfaceArea;

            surfp = 4.*pi.*rp.^2;
            
            scalingcoef = (vsf.*vol(1).*n.*F)./surfp;
            eqs{end + 1} = scalingcoef.*state.(ne).(co).(am).(sd).solidDiffusionEq;
            names{end + 1} = 'ei.ne_am_sd_soliddiffeq';
            
            eqs{end + 1}    = M_scaling .* scalingcoef.*state.(ne).(co).(am).(sd).massCons;
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

            
            %% Setup LinearizedProblem that can be processed by MRST Newton API
            
            problem = LinearizedProblem(eqs, types, names, primaryVars, state, dt);
            
        end

        %% cap concentrations and porosity to reasonnable values
        function [state, report] = updateState(model, state, problem, dx, drivingForces)

            [state, report] = updateState@Battery(model, state, problem, dx, drivingForces);
            
            ne = 'NegativeElectrode';
            am = 'ActiveMaterial';
            co = 'Coating';

            state.(ne).(co).(am).porosity = min(1, state.(ne).(co).(am).porosity);
            state.(ne).(co).(am).porosity = max(0, state.(ne).(co).(am).porosity);
            
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
