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
            paramobj.EffectiveElectronicConductivity = 'not used';
        end
        
    end
    
end
