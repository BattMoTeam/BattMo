classdef ElectrodeActiveComponentInputParams < ElectroChemicalComponentInputParams
%
% Input class for :class:`ElectrodeActiveComponent <Electrochemistry.ElectrodeActiveComponent>`
% 
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

        function paramobj = ElectrodeActiveComponentInputParams(jsonstruct)
            paramobj = paramobj@ElectroChemicalComponentInputParams(jsonstruct);

            pick = @(fd) pickField(jsonstruct, fd);
            paramobj.ActiveMaterial = ActiveMaterialInputParams(pick('ActiveMaterial'));
        end
        
    end
    
end
