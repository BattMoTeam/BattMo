classdef ElectrodeActiveComponentInputParams < ElectroChemicalComponentInputParams
 
    properties
        
        %% parameters for the electrode components
        ActiveMaterial
        amName
                
        % Interdiffusion coefficient parameter (diffusion between the particles)
        InterDiffusionCoefficient
        
        thermalConductivity
        heatCapacity
        
    end

    methods

        function paramobj = ElectrodeActiveComponentInputParams()
            paramobj = paramobj@ElectroChemicalComponentInputParams();
            paramobj.ActiveMaterial = ActiveMaterialInputParams();
            paramobj.amName = char();
            paramobj.EffectiveElectricalConductivity = 'not used';
        end
        
    end
    
end
