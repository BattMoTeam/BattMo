classdef CCcontrolModel < ControlModel

    
    methods

        function model = CCcontrolModel(inputparams)
            
            model = model@ControlModel(inputparams);
            
        end
        
        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@ControlModel(model);
            
            varnames = {};
            % Control type (string)
            % - 'constantCurrent'
            % - 'constantVoltage'
            varnames{end + 1} = 'ctrlType';            
            % control value that can be either a voltage or a current
            varnames{end + 1} = 'ctrlVal';            

            model = model.registerVarNames(varnames);
            
            fn = @CCDischargeControlModel.updateControlEquation;
            model = model.registerPropFunction({'controlEquation', fn, {'ctrlType', 'ctrlVal', 'E', 'I'}});
            
        end

        function state = updateControlEquation(model, state)
            
            E        = state.E;
            I        = state.I;            
            ctrlVal  = state.ctrlVal;
            ctrlType = state.ctrlType;

            switch ctrlType
              case 'constantCurrent'
                ctrleq = I - ctrlVal;
              case 'constantVoltage'
                %% TODO : fix hard-coded scaling
                ctrleq = (E - ctrlVal)*1e5;
              otherwise
                error('ctrlType not recognized');
            end
            
            state.controlEquation = ctrleq;
                        
        end
        
        function cleanState = addStaticVariables(model, cleanState, state)

            cleanState.ctrlType = state.ctrlType;
            
        end
        
        function step = setupScheduleStep(model, timeSteppingParams)
            
        % Setup and a return the step structure that is part of the schedule which is used as input for
        % :mrst:`simulateScheduleAD`. For some control type, there is a natural construction for this structure. This is
        % why we include this method here, for convenience. It can be overloaded by derived classes. The
        % timeSteppingParams structure by default is given by the data described in :battmofile:`Utilities/JsonSchemas/TimeStepping.schema.json`


            if (nargin > 1)
                params = timeSteppingParams;
            else
                params = [];
            end

            params = model.parseTimeSteppingStruct(params);
            
            if isa(model, 'CCChargeControlModel')
                rate = model.CRate;
            elseif isa(model, 'CCDischargeControlModel')
                rate = model.DRate;
            else
                error('model not recognized')
            end

            if ~isempty(params.totalTime)
                totalTime = params.totalTime;
            else
                totalTime = 1.4*(1*hour/rate);
            end

            if ~isempty(params.timeStepDuration)
                dt = params.timeStepDuration
            else
                n  = params.numberOfTimeSteps;
                dt = totalTime/n;
            end

            if params.useRampup
                n = params.numberOfRampupSteps;
            else
                n = 0;
            end

            dts = rampupTimesteps(totalTime, dt, n);
            
            step = struct('val', dts, 'control', ones(numel(dts), 1));

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
