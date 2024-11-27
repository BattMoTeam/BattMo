classdef ConstantFunction
% The function is a constant
    
    properties
        
        constant
        
    end

    methods

        function fn = Function(jsonstruct)

            fdnames = {'constant'};

            fn = dispatchParams(fn, jsonstruct, fdnames);

        end
        
        function y = eval(fn)

            y = fn.constant;
            
        end

        
    end
    
    
end
