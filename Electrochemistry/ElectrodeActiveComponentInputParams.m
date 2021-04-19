classdef ElectrodeActiveComponentInputParams < ThermoElectroChemicalComponentInputParams
 
    properties
        %% parameters for the electrode components
        % shortcut used here
        % am : ActiveMaterial parameters (class ActiveMaterialInputParams)
        am
        amName
    end

    methods

        function paramobj = ElectrodeActiveComponentInputParams()
            paramobj = paramobj@ThermoElectroChemicalComponentInputParams();
            paramobj.am = ActiveMaterialInputParams();
            paramobj.amName = char();
            paramobj.EffectiveElectronicConductivity = 'not used';
        end
        
    end
    
end
