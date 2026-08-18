classdef TimeControlModelInputParams < ControlModelInputParams

    properties

        value % function setup (see Utilities/JsonSchemas/Function.schema.json) that returns the value for the control,
              % either current or voltage, depending on the returned value of the type function
        type % function setup (see Utilities/JsonSchemas/Function.schema.json) that returns the value for the control
             % type,
             % - 1 for current
             % - 2 for voltage

        % Advanced parameters

        tolerance = 1e-4 % tolerance to skip timesteps (in second)

    end

    methods
        
        function inputparams = TimeControlModelInputParams(jsonstruct);

            inputparams =  inputparams@ControlModelInputParams(jsonstruct);
            
        end
        
    end
    
end

