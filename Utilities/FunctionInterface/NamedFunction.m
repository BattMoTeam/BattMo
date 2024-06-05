classdef NamedFunction
% The function is implemented and available in the path. It then called using its name.
    properties
        
        functionName

        %% helper

        functionHandler
    end

    methods


        function fn = Function(jsonstruct)

            fdnames = {'functionName'};

            fn = dispatchParams(fn, jsonstruct, fdnames);

            fn.functionHandler = str2func(fn.functionName);
            
        end
        
        function y = eval(fn, varargin)

            y = fn.functionHandler(varargin{:});
            
        end
        
    end
    
    
end
