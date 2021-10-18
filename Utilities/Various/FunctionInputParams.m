classdef FunctionInputParams < InputParams
    
    properties
        
        functionname
        filepath
        argumentlist
        name
        
    end
    
    methods

        function paramobj = FunctionInputParams(jsonstruct)
            paramobj = paramobj@InputParams(jsonstruct);
        end
        
    end
    
end
