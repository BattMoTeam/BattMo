classdef TimeControlModel < ControlModel


    properties

        value % function setup (see Utilities/JsonSchemas/Function.schema.json) that returns the value for the control,
              % either current or voltage, depending on the returned value of the type function
        type % function setup (see Utilities/JsonSchemas/Function.schema.json) that returns the value for the control
             % type,
             % - 1 for current
             % - 2 for voltage

        %% helper properties

        getValue % function object (see Utilities/FunctionInterface/Function.m) which is setup from value property
        getType % function object (see Utilities/FunctionInterface/Function.m) which is setup from type property

    end

    methods

        function model = TimeControlModel(inputparams)

            model = model@ControlModel(inputparams);

            fdnames = {'value', ...
                       'type'};

            model = dispatchParams(model, inputparams, fdnames);

            model.getValue = setupFunction(model.value);
            model.getType  = setupFunction(model.type);

        end

        function model = registerVarAndPropfuncNames(model)

            model = registerVarAndPropfuncNames@ControlModel(model);

            varnames = {};
            % Control type
            % - 1 : Current control
            % - 2 : Voltage control
            varnames{end + 1} = 'ctrlType';
            % control value that can be either a voltage or a current
            varnames{end + 1} = 'ctrlVal';

            model = model.registerVarNames(varnames);

            fn = @TimeControlModel.updateControlEquation;
            model = model.registerPropFunction({'controlEquation', fn, {'ctrlType', 'ctrlVal', 'E', 'I'}});

        end


        function [ctrlVal, ctrlType] = computeInput(model, t)

            ctrlVal  = model.getValue(t);
            ctrlType = model.getType(t);

        end

        function state = updateControlEquation(model, state)

            E        = state.E;
            I        = state.I;
            ctrlVal  = state.ctrlVal;
            ctrlType = state.ctrlType;

            switch ctrlType

              case 1

                ctrleq = I - ctrlVal;

              case 2

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

        function func = setupControlFunction(model)

            func = [];

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

            totalTime = getStructField(params, 'totalTime');

            if isa(totalTime, 'UnAssigned')
                if isa(model.getValue, 'TabulatedFunction1D')
                    totalTime = model.getValue.dataX(end);
                else
                    error('total time is not given and the input function is not a tabulated function');
                end
            end

            if isAssigned(params, {'timeStepDuration'})
                dt = params.timeStepDuration;
            else
                assert(isAssigned(params, 'numberOfTimeSteps'), 'No timeStepDuration and numberOfTimeSteps are given');
                n  = params.numberOfTimeSteps;
                dt = totalTime/n;
            end

            dts = rampupTimesteps(totalTime, dt, 0);

            step = struct('val', dts, 'control', ones(numel(dts), 1));

        end

    end


end

%{
Copyright 2021-2026 SINTEF Industry, Sustainable Energy Technology
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
