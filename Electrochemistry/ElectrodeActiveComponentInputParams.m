classdef ElectrodeActiveComponentInputParams < ElectroChemicalComponentInputParams
 
    properties
        %% parameters for the electrode components
        % shortcut used here
        % am : ActiveMaterial parameters (class ActiveMaterialInputParams)
        am
        amName
    end

    methods

        function paramobj = ElectrodeActiveComponentInputParams()
            paramobj = paramobj@ElectroChemicalComponentInputParams();
            paramobj.am = ActiveMaterialInputParams();
            paramobj.amName = char();
            paramobj.EffectiveElectronicConductivity = 'not used';
        end
        
    end
    
end
