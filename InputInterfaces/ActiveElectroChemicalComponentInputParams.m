classdef ActiveElectroChemicalComponentInputParams < ElectroChemicalComponentInputParams
 
    properties
        %% parameters for the electrode components
        % shortcut used here
        % am : ActiveMaterial parameters (class ActiveMaterialInputParams)
        am 
    end

    methods

        function paramobj = ActiveElectroChemicalComponentInputParams()
            paramobj = paramobj@ElectroChemicalComponentInputParams();
            paramobj.am = ActiveMaterialInputParams();
        end
        
        
        function paramobj = setup(parmobj, params)
        % params struct should contain valid fields for ElectroChemicalComponentInputParams
            paramobj = setup@ElectroChemicalComponentInputParams(paramobj, params);
        end

    end
    
end
