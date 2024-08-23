classdef RestControlModelInputParams < ControlModelInputParams
    
    methods

        function inputparams = RestControlModelInputParams(jsonstruct)

            inputparams = inputparams@ControlModelInputParams(jsonstruct);

            inputparams.controlPolicy = 'rest';
            
        end

    end
end
