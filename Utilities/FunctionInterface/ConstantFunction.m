classdef ConstantFunction < Function
% The function is a constant
    
    properties
        
        value
        
    end

    methods

        function fn = ConstantFunction(jsonstruct)

            fn = fn@Function(jsonstruct);
            
            fdnames = {'value'};

            fn = dispatchParams(fn, jsonstruct, fdnames);

        end
        
        function y = eval(fn, varargin)

            y = fn.value;
            
        end

        
    end
    
    
end
