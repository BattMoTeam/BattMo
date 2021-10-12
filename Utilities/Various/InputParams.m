classdef InputParams 
    
    methods

        function paramobj = InputParams(jsonstruct)
            paramobj = assignJsonParams(paramobj, jsonstruct);
        end
        
    end
    
end
