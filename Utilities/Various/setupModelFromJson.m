function [model, inputparams, gridGenerator, jsonstruct] = setupModelFromJson(jsonstruct)
% We setup a model from a json structure
    
    % We convert all the numerical value to SI unit.
    jsonstruct = resolveUnitInputJson(jsonstruct);
    
    inputparams = BatteryInputParams(jsonstruct);

    % Setup the geometry
    [inputparams, gridGenerator] = setupBatteryGridFromJson(inputparams, jsonstruct);

    % Check the consistency of the parameters.
    inputparams = inputparams.validateInputParams();

    model = Battery(inputparams);
    
end
