classdef ElectroChemicalComponentInputParams < ElectronicComponentInputParams
    
    properties
       
        ionName
        ionFluxName 
        ionSourceName
        ionMassConsName
        ionAccumName        

    end
    
    methods
        
        function paramobj = setup(parmobj, params)
        % params struct should contain valid fields for ElectronicComponentInputParams and the fields
        % - ionName
        % - ionFluxName 
        % - ionSourceName
        % - ionMassConsName
        % - ionAccumName        
            
            paramobj.ionName         = params.ionName;
            paramobj.ionFluxName     = params.ionFluxName;
            paramobj.ionSourceName   = params.ionSourceName;
            paramobj.ionMassConsName = params.ionMassConsName;
            paramobj.ionAccumName    = params.ionAccumName;
            
            paramobj = setup@ElectronicComponentInputParams(paramobj, params);
            
        end
        
    end
    
end
