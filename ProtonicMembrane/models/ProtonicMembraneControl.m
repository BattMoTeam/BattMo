classdef ProtonicMembraneControl < BaseModel
    
    properties

        controlType
        useCurrentDensity % If true, then we use the value of currentDensity. Otherwise, a value for I the total current should
                          % be given.
        currentDensity
        I
        
        area % if current density is used, we need the area. This value will typically be assigned by the grid constructor
        
    end
    
    
    methods
        
        function model = ProtonicMembraneControl(inputparams)

            model = model@BaseModel();

            fdnames = {'controlType'      , ...
                       'useCurrentDensity', ...
                       'currentDensity'   , ...
                       'I'                , ...
                       'area'};
            model = dispatchParams(model, inputparams, fdnames);

            if model.useCurrentDensity
                model.I = model.currentDensity*model.area;
            else
                model.currentDensity = model.I/model.area;
            end
            
        end

        function model = registerVarAndPropfuncNames(model)

            varnames = {};

            % Potential
            varnames{end + 1} = 'U';
            % Current
            varnames{end + 1} = 'I';            
            % Value of the variable that is controled
            varnames{end + 1} = 'ctrlVal';
            % Control Equation
            varnames{end + 1} = 'controlEquation';

            model = model.registerVarNames(varnames);
            
            fn = @ProtonicMembraneControl.setupControl;
            switch model.controlType
              case 'current'
                inputnames = {'I', 'ctrlVal'};
              case 'voltage'
                inputnames = {'U', 'ctrlVal'};
              otherwise
                error('controlType not recognized');
            end
            model = model.registerPropFunction({'controlEquation', fn, inputnames});

        end

        function state = setupControl(model, state)

            switch model.controlType
              case 'current'
                state.controlEquation = state.I - state.ctrlVal;
              case 'voltage'
                state.controlEquation = state.U - state.ctrlVal;
              otherwise
                error('control type not recognized');
            end
            
        end

        function step = setupScheduleStep(model, jsonstruct)

            N = jsonstruct.TimeStepping.numberOfTimeSteps;

            dt = 1/N;
            
            step.val     = dt*ones(N, 1);
            step.control = ones(numel(step.val), 1);
            
        end

        function control = setupScheduleControl(model, jsonstruct)

            useSwitch      = jsonstruct.TimeStepping.useSwitch;
            fractionSwitch = jsonstruct.TimeStepping.fractionSwitch;
            orderSwitch    = jsonstruct.TimeStepping.orderSwitch;
            
            I = model.I;
            
            control.src = @(time, I) pmControlFunc(time, I, fractionSwitch, 1, 'order', orderSwitch);
            
        end
        
        function schedule = setupSchedule(model, jsonstruct)
        % Convenience function to setup schedule from main jsonstruct with property TimeStepping

            step    = model.setupScheduleStep(jsonstruct);
            control = model.setupScheduleControl(jsonstruct);
            
            schedule = struct('step', step, ...
                              'control', control);
            
        end

    end
    
end
