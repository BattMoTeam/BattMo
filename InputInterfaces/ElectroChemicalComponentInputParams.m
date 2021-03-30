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
            
            fdnames = {'ionName', ...
                       'ionFluxName', ...
                       'ionSourceName', ...
                       'ionMassConsName', ...
                       'ionAccumName'};
            paramobj = dispatchParams(paramobj, params, fdnames);

            paramobj = setup@ElectronicComponentInputParams(paramobj, params);
            
        end
        
    end
    
end
