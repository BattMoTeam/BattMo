classdef ActiveElectroChemicalComponentInputParams < ElectroChemicalComponentInputParams
% for the moment, not much here
    
    methods

        function paramobj = setup(parmobj, params)
        % params struct should contain valid fields for ElectroChemicalComponentInputParams
            paramobj = setup@ElectroChemicalComponentInputParams(paramobj, params);
        end

    end
    
end
