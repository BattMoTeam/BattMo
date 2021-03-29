classdef ElectroChemicalComponentInputParams < ElectronicComponentInputParams
    
    properties
       
        ionName
        ionFluxName 
        ionSourceName
        ionMassConsName
        ionAccumName        

    end
    
    methods
        
        function paramobj = ElectroChemicalComponentInputParams(params)
        % params struct should contain valid fields for ElectronicComponentInputParams and the fields
        % - ionName
        % - ionFluxName 
        % - ionSourceName
        % - ionMassConsName
        % - ionAccumName        
            
            paramobj = paramobj@ElectronicComponentInputParams(params);
            paramobj.ionName         = params.ionName;
            paramobj.ionFluxName     = params.ionFluxName;
            paramobj.ionSourceName   = params.ionSourceName;
            paramobj.ionMassConsName = params.ionMassConsName;
            paramobj.ionAccumName    = params.ionAccumName;
            
        end
        
    end
    
end
