classdef HalfCell < BaseModel
% 
% Used when at least one of the two electrodes is a swelling Material.
% It functions as the Battery class but implements a new equation necessary because of the particle swelling :
% the volume Conservation equation 
%
    properties

        Electrolyte
        ActiveMaterial
        Control

        externalCoupling
        
        initT
        
        primaryVarNames
        funcCallList
        
    end
    
    methods
        
        function model = HalfCell(paramobj)
            
            model = model@BaseModel()

            model.Electrolyte    = SingleCellElectrolyte(paramobj.Electrolyte);
            model.ActiveMaterial = ActiveMaterial(paramobj.ActiveMaterial);
            model.Control        = IEswitchControlModel(paramobj.Control);
            
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


            %% We call the assembly equations ordered from the graph
            
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
        
        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@BaseModel(model);

            elyte = 'Electrolyte';
            itf   = 'Interface';
            am    = 'ActiveMaterial';
            sd    = 'SolidDiffusion';
            ctrl  = 'Control';

            model = model.removeVarName({am, itf, 'dUdT'});
            model = model.removeVarName({elyte, 'massCons'});
            
            fn = @HalfCell.updateElectrolyte;
            inputnames = {{elyte, 'c'}, {elyte, 'phi'}};
            model = model.registerPropFunction({{am, itf, 'cElectrolyte'}, fn, inputnames});
            model = model.registerPropFunction({{am, itf, 'phiElectrolyte'}, fn, inputnames});

            fn = @HalfCell.updateTemperature;
            model = model.registerPropFunction({{am, 'T'}, fn, {}});
            
            fn = @HalfCell.setupEIEquation;
            inputnames = {{ctrl, 'E'}, {ctrl, 'I'}, {am, 'phi'}};
            model = model.registerPropFunction({{ctrl, 'EIequation'}, fn, inputnames});

            inputnames = {};
            fn = @Battery.updateControl;
            fn = {fn, @(propfunction) PropFunction.drivingForceFuncCallSetupFn(propfunction)};
            model = model.registerPropFunction({{ctrl, 'ctrlVal'}, fn, inputnames});            
            model = model.registerPropFunction({{ctrl, 'ctrlType'}, fn, inputnames});
            
            
        end

        function model = updateElectrolyte(model, state)

            elyte = 'Electrolyte';
            am    = 'ActiveMaterial';
            itf   = 'Interface';

            state.(am).(itf).cElectrolyte = state.(elyte).c;
            state.(am).(itf).phiElectrolyte = state.(elyte).phi;
           
        end
        
        function model = updateTemperature(model, state)

            am    = 'ActiveMaterial';
            
            state.(am).T = model.initT;
            
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
