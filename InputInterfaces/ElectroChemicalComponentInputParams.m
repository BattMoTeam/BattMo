classdef ElectroChemicalComponentInputParams < ElectronicComponentInputParams
    
    properties
       
        ionName
        ionFluxName 
        ionSourceName
        ionMassConsName
        ionAccumName        

    end
    
    methods
        
        function params = ElectroChemicalComponentInputParams(params)
            params = params@ElectronicComponentInputParams(params);
        end
        
    end
    
end
